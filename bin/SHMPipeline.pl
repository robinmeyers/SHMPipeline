#!/usr/bin/perl

##
## This program analyzes deletions and insertions
## in NGS amplicon sequence data
##
## run with "--help" for usage information
##
## Robin Meyers <robin.meyers@childrens.harvard.edu>, 05mar2013

use strict;
use warnings;
use Getopt::Long;
use Carp;
use IO::File;
use Text::CSV;
use Bio::SeqIO;
use Bio::DB::Sam;
use Bio::Seq::Quality;
use Switch;
use List::MoreUtils qw(pairwise);
use threads;
use threads::shared;
use IPC::System::Simple qw(capture);
use Interpolation 'arg:@->$' => \&argument;
use Time::HiRes qw(gettimeofday tv_interval);

use Cwd qw(abs_path);
use FindBin;
use lib abs_path("$FindBin::Bin/../lib");

require "SHMHelper.pl";
require "pslHelper.pl";


# Flush output after every write
select( (select(STDOUT), $| = 1 )[0] );



# Forward declarations
sub parse_command_line;
sub read_in_meta_file;
sub check_existance_of_files;
sub initialize_stats_hash;
sub process_illumina_experiment ($);
sub align_to_reference ($);
sub merge_alignments ($);
sub write_summary_stats;
sub create_summary;


# Global flags and arguments, 
# Set by command line arguments
my $meta_file;
my $indir;
my $outdir;
my $refdir;


my $user_bowtie_opt = "";
my $max_threads = 2;
my $expt_threads = 4;
my $phred;
my $min_qual = 20;
my $min_sw_score = 100;
my $dup_threshold = 0.9;
my $shm_threshold = 0;
my $ow;

# Global variabless
my %meta_hash;


my $exptfile;
my $shmexptfile;
my $delshmexptfile;
my $clonefile;
my $shmclonefile;
my $delshmclonefile;
my $bxexptfile;
my $mutfile;
#
# Start of Program
#

parse_command_line;

my $default_pe_bowtie_opt = "--local -D 20 -R 3 -N 1 -L 12 -i C,6 --rfg 8,1 --rdg 8,1 --score-min C,100 --no-discordant --no-mixed -p $expt_threads --reorder -t";
my $default_merge_bowtie_opt = "-D 20 -R 3 -N 1 -L 12 -i C,6 --np 0 --rfg 8,1 --rdg 8,1 --score-min C,-200 -p $expt_threads --reorder -t";

my $t0 = [gettimeofday];

my $pe_bt2_opt = manage_program_options($default_pe_bowtie_opt,$user_bowtie_opt);
my $merge_bt2_opt = manage_program_options($default_merge_bowtie_opt,$user_bowtie_opt);

read_in_meta_file;

check_existance_of_files;

#initialize_stats_hash;

my $vizdir = $outdir."/viz";
mkdir $vizdir;

my @threads = ();

foreach my $expt_id (sort keys %meta_hash) {

    while (1) {

    # joins any threads if possible
        foreach my $thr (@threads) {
            $thr->join() if $thr->is_joinable();
        }

        my @running = threads->list(threads::running);
        
        # if there are open threads, create a new one, push it onto list, and exit while loop
        if (scalar @running < $max_threads) {
            my $thr = threads->create( sub {
                        my $t0_expt = [gettimeofday];
                        print "\nStarting $expt_id\n";
                        
                        process_illumina_experiment( $meta_hash{$expt_id} );
                        my $t1 = tv_interval($t0_expt);
                        printf("\nFinished %s in %.2f seconds.\n", $expt_id,$t1);
                    });
            push(@threads,$thr);
            sleep(1);
            last;
        }
        sleep(1);
    } 
}

# waits for all threads to finish
while( scalar threads->list(threads::all) > 0) {
  for my $thr (@threads) {
      $thr->join() if $thr->is_joinable;
  }
  sleep(1);
}

# create_summary;

my $t1 = tv_interval($t0);

printf("\nFinished all processes in %.2f seconds.\n", $t1);


#
# End of program
#




sub process_illumina_experiment ($) {

  my $expt = shift;
  my $expt_id = $expt->{experiment};
  my $exptdir = "$outdir/$expt_id";

  mkdir $exptdir unless -d $outdir or croak "Error: cannot create experiment directory";


  $expt->{logfh} = IO::File->new(">".$expt->{log});
  *STDOUT = $expt->{logfh};
  *STDERR = $expt->{logfh};

  align_to_reference($expt);

  merge_alignments($expt);

  realign_to_reference($expt);

  parse_alignments($expt);


  print Capture(join(" ","Rscript $FindBin::Bin/../R/removeDupClones.R",
                         $expt->{mutfile},$expt->{readfile},$expt->{reference},
                         "cores=$expt_threads",
                         "j_thresh=$dup_threshold",
                         "tstart=".$expt->{start},
                         "tend=".$expt->{end},
                         "2>&1"));


  print Capture(join(" ","Rscript $FindBin::Bin/../R/SHMProfile.R",
                         $expt->{mutfile},$expt->{readfile},$expt->{reference},"$exptdir/$expt_id",
                         "tstart=".$expt->{start},
                         "tend=".$expt->{end},
                         "2>&1"));



  print Capture(join(" ","Rscript $FindBin::Bin/../R/SHMPlot.R",
                         "$exptdir/${expt_id}_profile.txt",
                         "$vizdir/${expt_id}_mutrate.pdf",
                         "tstart=".$expt->{start},
                         "tend=".$expt->{end},
                         "2>&1"));

  print Capture(join(" ","Rscript $FindBin::Bin/../R/SHMViz.R",
                         $expt->{mutfile},$expt->{clonefile},$expt->{reference},
                         "$vizdir/${expt_id}_viz.pdf",
                         "tstart=".$expt->{start},
                         "tend=".$expt->{end},
                         "2>&1"));



}

sub parse_alignments ($) {
  use Switch;
  my $expt = shift;

  my $samobj = Bio::DB::Sam->new(-bam => $expt->{bam},
                                     -fasta => $expt->{reference},
                                     -expand_flags => 1);

  my $align_stream = $samobj->get_seq_stream;

  my $mutfh = IO::File->new(">".$expt->{mutfile});
  my $readfh = IO::File->new(">".$expt->{readfile});

  

  # my @bases = qw(A C G T N);
  # my @bx_header = ("ID");
  # foreach my $b1 (@bases) {
  #   next if $b1 eq "N";
  #   foreach my $b2 (@bases) {
  #     next if $b2 eq "N";
  #     next if $b1 eq $b2;
  #     push(@bx_header,$b1."->".$b2);
  #   }
  # }
  # $bxfh->print(join("\t",@bx_header)."\n");




  while (my $Aln = $align_stream->next_seq) {
    next if $Aln->unmapped;


    my $read = {};
    $read->{bps} = 0;
    # my $clone = {};
    # $clone->{bps} = 0;
    # $clone->{sub} = 0;
    # $clone->{del} = 0;
    # $clone->{delbp} = 0;
    # $clone->{largedel} = 0;
    # $clone->{ins} = 0;
    # $clone->{insbp} = 0;

    
    # foreach my $base1 (@bases) {
    #   $read->{$base1} = 0;
    #   foreach my $base2 (@bases) {
    #     if ($base1 eq $base2) {
    #       $read->{bx}->{$base1}->{$base2} = "-";
    #     } else {
    #       $read->{bx}->{$base1}->{$base2} = 0;
    #     }
    #   }
    # }

    my @Rseq = split("",uc($samobj->seq($Aln->seq_id)));
    my @Qseq = split("",uc($Aln->query->seq->seq));
    my @Qual = $Aln->qscore;
    my @Cigar = ();
    foreach my $i (@{$Aln->cigar_array}) {
      next if $i->[0] eq "S";
      push(@Cigar,split("",($i->[0] x $i->[1])));
    }
    my $Rpos = $Aln->start;
    my $Qpos = $Aln->query->start;

    my @coords;
    my $readstart;
    my $readend;

    my $refstart = $expt->{start};
    my $refend = $expt->{end};

    while ($Rpos <= $Aln->end) {
      my $c = shift(@Cigar);
      croak "Error: cigar error" unless defined $c;
      switch ($c) {
        case "M" {

          if ($Rpos < $refstart || $Rpos > $refend || 
            ($Qual[$Qpos-1] < $min_qual && $Qseq[$Qpos-1] ne $Rseq[$Rpos-1])) {
            if (defined $readstart) {
              $readend = $Rpos-1;
            }
            
          } else {

            $readstart = $Rpos unless defined $readstart;
            $read->{bps}++;
            if ($Qseq[$Qpos-1] ne $Rseq[$Rpos-1] && $Qseq[$Qpos-1] ne 'N') {
              $mutfh->print(join("\t",$expt->{experiment},
                                      $Aln->qname,
                                      $Rpos,
                                      "sub",
                                      $Rseq[$Rpos-1],
                                      $Qseq[$Qpos-1])."\n");
              
              # $read->{bx}->{$Rseq[$Rpos-1]}->{$Qseq[$Qpos-1]}++;
            }
            
          }

          $Rpos++;
          $Qpos++;

        }
        case "D" {

          if (defined $readstart) {
            $readend = $Rpos-1;
          }

          my $del_size = 0;
          my $del_start = $Rpos;
          do {
            $del_size++;
            $Rpos++;
            $c = shift(@Cigar);
          } while ($c eq 'D');
          unshift(@Cigar,$c);


          unless ($del_start < $refstart || $Rpos > $refend) {

            $mutfh->print(join("\t",$expt->{experiment},
                                    $Aln->qname,
                                    $del_start,
                                    "del",
                                    "",
                                    "",
                                    $del_size,
                                    $Rpos-1)."\n");
          }
        }

        case "I" {
          my $ins_bases = "";
          do {
            $ins_bases .= $Qseq[$Qpos-1];
            $Qpos++;
            $c = shift(@Cigar);
          } while ($c eq "I");
          unshift(@Cigar,$c);
          
          unless ($Rpos < $refstart || $Rpos > $refend) {
            $mutfh->print(join("\t",$expt->{experiment},
                                    $Aln->qname,
                                    $Rpos,
                                    "ins",
                                    "",
                                    "",
                                    length($ins_bases),
                                    "",
                                    $ins_bases)."\n");
          }
        }
      }

      if (defined $readstart && defined $readend) {
        push(@coords,"$readstart-$readend");
        # $readfh->print(join("\t",$expt->{expt},
        #                          $Aln->qname,
        #                          $readstart,
        #                          $readend)."\n");
        undef $readstart;
        undef $readend;
      }

    }
    if (defined $readstart) {
      $readend = $refend;
      push(@coords,"$readstart-$readend");
    }

    $readfh->print(join("\t",$expt->{experiment},
                             $Aln->qname,
                             $read->{bps},
                             join(",",@coords))."\n");

    # my @bx_row = ($Aln->qname);
    # foreach my $b1 (@bases) {
    #   next if $b1 eq "N";
    #   foreach my $b2 (@bases) {
    #     next if $b2 eq "N";
    #     next if $b1 eq $b2;
    #     push(@bx_row,$read->{bx}->{$b1}->{$b2});
    #   }
    # }
    # $bxfh->print(join("\t",@bx_row)."\n");

  }

}

sub merge_alignments ($) {

  use Switch;

  my $expt = shift;

  my $samobj = Bio::DB::Sam->new(-bam => $expt->{pe_bam},
                                     -fasta => $expt->{reference},
                                     -expand_flags => 1);

  my $align_stream = $samobj->get_seq_stream(-type => 'read_pair');

  my $merged_fq = $expt->{exptdir}."/".$expt->{experiment}."_merged.fq";
  $expt->{merged_fq} = $merged_fq;

  my $merged_fq_fh = Bio::SeqIO->new(-file => ">$merged_fq", -format => 'fastq');

  while (my $pair = $align_stream->next_seq) {
    my ($Aln1,$Aln2) = $pair->get_SeqFeatures;

    next if $Aln1->unmapped || $Aln2->unmapped;
    
    my @Rseq = split("",$samobj->seq($Aln1->seq_id));

    my @Qseq1 = split("",$Aln1->query->seq->seq);
    my @Qseq2 = split("",$Aln2->query->seq->seq);
    my @Qual1 = $Aln1->qscore;
    my @Qual2 = $Aln2->qscore;
    my $Start1 = $Aln1->start;
    my $End1 = $Aln1->end;
    my $Start2 = $Aln2->start;
    my $End2 = $Aln2->end;

    if ($Start2 < $Start1) {
      carp "Warning: R2 start is less than R1 start";
      next;
    }

    my @Cigar1 = ();
    foreach my $i (@{$Aln1->cigar_array}) {
      next if $i->[0] eq "S";
      push(@Cigar1,split("",($i->[0] x $i->[1])));
    }
    my @Cigar2 = ();
    foreach my $i (@{$Aln2->cigar_array}) {
      next if $i->[0] eq "S";
      push(@Cigar2,split("",($i->[0] x $i->[1])));
    }

    
    my $Rpos = $Aln1->start;
    my $Qpos1 = $Aln1->query->start;
    my $Qpos2 = $Aln2->query->start;

    my @merged_seq = ();
    my @merged_qual = ();

    # if ($Qpos1 > 1) {
    #   push(@merged_seq,@Qseq1[0..($Qpos1-2)]);
    #   push(@merged_qual,@Qual1[0..($Qpos1-2)]);
    # }

    while ($Rpos <= $Aln1->end && $Rpos < $Aln2->start) {
      
      my $c1 = shift(@Cigar1);
      switch ($c1) {
        case 'M' {
          push(@merged_seq,$Qseq1[$Qpos1-1]);
          push(@merged_qual,$Qual1[$Qpos1-1]);
          $Qpos1++;
          $Rpos++;
        }
        case 'D' {
          $Rpos++;
        }
        case 'I' {
          push(@merged_seq,$Qseq1[$Qpos1-1]);
          push(@merged_qual,$Qual1[$Qpos1-1]);
          $Qpos1++;
        }
      }
    }

    while ($Rpos > $Aln1->end && $Rpos < $Aln2->start) {
      push(@merged_seq,"N");
      push(@merged_qual,2);
      $Rpos++;
    }

    while ($Rpos <= $Aln1->end && $Rpos >= $Aln2->start) {
      my $c1 = shift @Cigar1;
      my $c2 = shift @Cigar2;
      my $q1 = $Qseq1[$Qpos1-1];
      my $q2 = $Qseq2[$Qpos2-1];
      my $r = $Rseq[$Rpos-1];

      switch ($c1) {
        case "M" {
          switch ($c2) {
            case "M" {
              if ($q1 eq $q2) {
                push(@merged_seq,$q1);
                push(@merged_qual,max($Qual1[$Qpos1-1],$Qual2[$Qpos2-1]));
              } elsif ($q1 eq $r) {
                push(@merged_seq,$q1);
                push(@merged_qual,$Qual1[$Qpos1-1])
              } elsif ($q2 eq $r) {
                push(@merged_seq,$q2);
                push(@merged_qual,$Qual2[$Qpos2-1])
              } else {
                push(@merged_seq,"N");
                push(@merged_qual,2)
              }
              $Qpos1++;
              $Qpos2++;
              $Rpos++;
            }
            case "D" {
              if ($q1 eq $r) {
                push(@merged_seq,$q1);
                push(@merged_qual,$Qual1[$Qpos1-1])
              }
              $Qpos1++;
              $Rpos++;
            }
            case "I" {
              $Qpos2++;
              unshift(@Cigar1,$c1);
            }
          }
        }
        case "I" {
          switch ($c2) {
            case "M" {
              $Qpos1++;
              unshift(@Cigar2,$c2);
            }
            case "D" {
              $Qpos1++;
              unshift(@Cigar2,$c2);
            }
            case "I" {
              if ($Qual1[$Qpos1-1] > $Qual2[$Qpos2-1]) {
                push(@merged_seq,$q1);
                push(@merged_qual,$Qual1[$Qpos1-1])
              } else {
                push(@merged_seq,$q2);
                push(@merged_qual,$Qual2[$Qpos2-1])
              }
              $Qpos1++;
              $Qpos2++;
            }
          }
        }
        case "D" {
          switch ($c2) {
            case "M" {
              if ($q2 eq $r) {
                push(@merged_seq,$q2);
                push(@merged_qual,$Qual2[$Qpos2-1])
              }
              $Qpos2++;
              $Rpos++;
            }
            case "D" {
              $Rpos++;
            }
            case "I" {
              $Qpos2++;
              unshift(@Cigar1,$c1);
            }
          }
        }
      }
    }

    while ($Rpos <= $Aln2->end) {
      my $c2 = shift(@Cigar2);
      switch ($c2) {
        case "M" {
          push(@merged_seq,$Qseq2[$Qpos2-1]);
          push(@merged_qual,$Qual2[$Qpos2-1]);
          $Qpos2++;
          $Rpos++;
        }
        case "D" {
          $Rpos++;
        }
        case "I" {
          push(@merged_seq,$Qseq2[$Qpos2-1]);
          push(@merged_qual,$Qual2[$Qpos2-1]);
          $Qpos2++;
        }
      }
    }

    # if ($Qpos2 < @Qseq2) {
    #   push(@merged_seq,@Qseq2[$Qpos2..$#Qseq2]);
    #   push(@merged_qual,@Qual1[$Qpos2..$#Qseq2]);
    # }

    my $fq = Bio::Seq::Quality->new(-id => $Aln1->qname, -seq => join("",@merged_seq), -qual => \@merged_qual);

    $merged_fq_fh->write_seq($fq);

  }

}

sub align_to_reference ($) {

  my $expt = shift;

  my $ref_fa = $expt->{reference};
  my $R1_fa = $expt->{R1};
  my $R2_fa = $expt->{R2};
  my $pe_sam = $expt->{exptdir}."/tmp/".$expt->{experiment}."_pe.sam";
  $expt->{pe_sam} = $pe_sam;
  (my $pe_bam = $pe_sam) =~ s/\.sam$/.bam/;
  $expt->{pe_bam} = $pe_bam;
  my $log = $expt->{log};

  mkdir $expt->{exptdir}."/tmp" unless -d $expt->{exptdir}."/tmp";
  print "\nRunning PE Bowtie2 alignment for ".$expt->{experiment}." against reference sequence\n";

  print Capture("bowtie2-build -q $ref_fa $ref_fa 2>&1") unless -r "$ref_fa.1.bt2";

  my $bt2_cmd = "bowtie2 $pe_bt2_opt -x $ref_fa -1 $R1_fa -2 $R2_fa -S $pe_sam 2>&1";

  print Capture($bt2_cmd);

  print Capture("samtools view -bS -o $pe_bam $pe_sam 2>&1");
}

sub realign_to_reference ($) {

  my $expt = shift;

  my $ref_fa = $expt->{reference};
  my $merged_fq = $expt->{merged_fq};
  my $sam = $expt->{exptdir}."/".$expt->{experiment}.".sam";
  $expt->{sam} = $sam;
  (my $bam = $sam) =~ s/\.sam$/.bam/;
  $expt->{bam} = $bam;
  my $log = $expt->{log};

  mkdir $expt->{exptdir}."/tmp" unless -d $expt->{exptdir}."/tmp";
  print "\nRunning merged Bowtie2 alignment for ".$expt->{experiment}." against reference sequence\n";

  my $bt2_cmd = "bowtie2 $merge_bt2_opt -x $ref_fa -U $merged_fq -S $sam 2>&1";

  print Capture($bt2_cmd);

  print Capture("samtools view -bS -o $bam $sam 2>&1");
}

sub check_existance_of_files {
	print "\nSearching for sequence files...\n";


	foreach my $expt_id (sort keys %meta_hash) {
		my $input = $indir."/".$expt_id;

    if (-r "${input}_R1.fastq" && -r "${input}_R2.fastq") {
      $meta_hash{$expt_id}->{R1} = "${input}_R1.fastq";
      $meta_hash{$expt_id}->{R2} = "${input}_R2.fastq";
      next;

    } else {
      croak "Error: No valid experiment $expt_id in $indir" ;  
    }
	}
	

	print "\nSearching for reference files...\n";

	foreach my $expt_id (sort keys %meta_hash) {
		my $reffile = $refdir."/".$meta_hash{$expt_id}->{reference};
		croak "Error: Could not locate reference file $reffile in $indir" unless (-r $reffile);
		$meta_hash{$expt_id}->{reference} = $reffile;
	}
	
	print "Done.\n";

	print "\nChecking on output files...\n";
	foreach my $expt_id (sort keys %meta_hash) {
    my $exptdir = $outdir."/$expt_id/";
    $meta_hash{$expt_id}->{exptdir} = $outdir."/$expt_id";
    mkdir $meta_hash{$expt_id}->{exptdir};
		my $mutfile = $exptdir."/".$expt_id."_muts.txt";
    my $clonefile = $exptdir."/".$expt_id."_clones.txt";
    my $readfile = $exptdir."/".$expt_id."_reads.txt";
    my $bxfile = $exptdir."/".$expt_id."_basex.txt";
    my $logfile = $exptdir."/".$expt_id.".log";
		croak "Error: Output file $mutfile already exists and overwrite switch --ow not set" if (-r $mutfile && ! defined $ow);
		System("touch $mutfile",1);
		croak "Error: Cannot write to $mutfile" unless (-w $mutfile);
		$meta_hash{$expt_id}->{mutfile} = $mutfile;
    $meta_hash{$expt_id}->{clonefile} = $clonefile;
    $meta_hash{$expt_id}->{readfile} = $readfile;
    $meta_hash{$expt_id}->{bxfile} = $bxfile;
    $meta_hash{$expt_id}->{log} = $logfile;

    System("echo \"Running experiment $expt_id\" > $logfile",1);

	}
	$exptfile = "$outdir/Expts.txt";
  $clonefile = "$outdir/Clones.txt";
  $bxexptfile = "$outdir/BaseX.txt";
  $mutfile = "$outdir/Mutations.txt";

  croak "Error: Output file $exptfile already exists and overwrite switch --ow not set" if (-r $exptfile && ! defined $ow);
  System("touch $exptfile",1);
  croak "Error: Cannot write to $exptfile" unless (-w $exptfile);
  print "Done.\n";

}

sub create_summary {

  my %stats;
  my $exptsfh = IO::File->new(">$exptfile");
  my $clonesfh = IO::File->new(">$clonefile");
  $exptsfh->print(join("\t",qw(Expt Allele Clones Bp Subs Del DelBp Ins InsBp RefA RefC RefG RefT RefN))."\n");
  $clonesfh->print(join("\t",qw(Expt Allele Clone Bp Subs Del DelBp LargeDel Ins InsBp RefA RefC RefG RefT RefN Coords))."\n");



  my %bx;
  my $bxexptsfh = IO::File->new(">$bxexptfile");
  my @bx_combs = ();
  my @bases = qw(A C G T);
  foreach my $b1 (@bases) {
    foreach my $b2 (@bases) {
      next if $b1 eq $b2;
      push(@bx_combs,$b1."->".$b2);
    }
  }
  $bxexptsfh->print(join("\t","Expt","Allele",@bx_combs)."\n");

  foreach my $expt (sort keys %meta_hash) {
    my $clonefh = IO::File->new("<".$meta_hash{$expt}->{clonefile});
    my $csv = Text::CSV->new({sep_char => "\t"});
    my $header = $csv->getline($clonefh);
    $csv->column_names(@$header);

    my $bxfh = IO::File->new("<".$meta_hash{$expt}->{bxfile});
    my $bx_csv = Text::CSV->new({sep_char => "\t"});
    my $bx_header = $bx_csv->getline($bxfh);
    $bx_csv->column_names(@$bx_header);

    @{$bx{$expt}}{@bx_combs} = (0) x @bx_combs;

    while (my $clone = $bx_csv->getline_hr($bxfh)) {
      $bx{$expt}->{$clone->{ID}} = $clone;
    }
    $bxfh->close;


    $stats{$expt}->{Clones} = 0;
    $stats{$expt}->{Bp} = 0;
    $stats{$expt}->{Subs} = 0;
    $stats{$expt}->{Del} = 0;
    $stats{$expt}->{DelBp} = 0;
    $stats{$expt}->{Ins} = 0;
    $stats{$expt}->{InsBp} = 0;
    $stats{$expt}->{RefA} = 0;
    $stats{$expt}->{RefC} = 0;
    $stats{$expt}->{RefG} = 0;
    $stats{$expt}->{RefT} = 0;
    $stats{$expt}->{RefN} = 0;


    while (my $clone = $csv->getline_hr($clonefh)) {

      next unless $clone->{Bp} > 0;
      next if $clone->{LargeDel} > 10 && $clone->{Subs} < 2;
      $stats{$expt}->{Clones}++;
      $stats{$expt}->{Bp} += $clone->{Bp};
      $stats{$expt}->{Subs} += $clone->{Subs};
      $stats{$expt}->{Del} += $clone->{Dels};
      $stats{$expt}->{DelBp} += $clone->{DelBp};
      $stats{$expt}->{Ins} += $clone->{Ins};
      $stats{$expt}->{InsBp} += $clone->{InsBp};
      $stats{$expt}->{RefA} += $clone->{RefA};
      $stats{$expt}->{RefC} += $clone->{RefC};
      $stats{$expt}->{RefG} += $clone->{RefG};
      $stats{$expt}->{RefT} += $clone->{RefT};
      $stats{$expt}->{RefN} += $clone->{RefN};

      my @tmp1 = @{$bx{$expt}}{@bx_combs};
      my @tmp2 = @{$bx{$expt}->{$clone->{ID}}}{@bx_combs};
      # print "$expt ".$clone->{ID}.":\n";
      # print "@tmp1\n";
      # print "@tmp2\n";
      our ($a,$b);
      @{$bx{$expt}}{@bx_combs} = pairwise { $a + $b } @tmp1, @tmp2;



      $clonesfh->print(join("\t",$expt,$meta_hash{$expt}->{allele},
                                $clone->{ID},
                                $clone->{Bp},
                                $clone->{Subs},
                                $clone->{Dels},
                                $clone->{DelBp},
                                $clone->{LargeDel},
                                $clone->{Ins},
                                $clone->{InsBp},
                                $clone->{RefA},
                                $clone->{RefC},
                                $clone->{RefG},
                                $clone->{RefT},
                                $clone->{RefN},
                                $clone->{Coords})."\n");
    }



    $clonefh->close;

    $exptsfh->print(join("\t",$expt,$meta_hash{$expt}->{allele},$stats{$expt}->{Clones},
                                  $stats{$expt}->{Bp},
                                  $stats{$expt}->{Subs},
                                  $stats{$expt}->{Del},                                  
                                  $stats{$expt}->{DelBp},
                                  $stats{$expt}->{Ins},
                                  $stats{$expt}->{InsBp},                                  
                                  $stats{$expt}->{RefA},
                                  $stats{$expt}->{RefC},
                                  $stats{$expt}->{RefG},
                                  $stats{$expt}->{RefT},
                                  $stats{$expt}->{RefN})."\n");


    $bxexptsfh->print(join("\t",$expt,$meta_hash{$expt}->{allele},@{$bx{$expt}}{@bx_combs})."\n");

    
  }

  $exptsfh->close;
  $clonesfh->close;
  $bxexptsfh->close;
  System("cat $outdir/*/*_muts.txt > $mutfile");

  mkdir "$outdir/group_viz";
  my $viz_cmd = "Rscript $FindBin::Bin/../R/mutationVizGrouped.R $meta_file $mutfile $clonefile $refdir $outdir/group_viz/";
  System("$viz_cmd > $outdir/group_viz/R.out 2>&1");

  my $group_cmd = "Rscript $FindBin::Bin/../R/mutationStats.R $exptfile $clonefile $meta_file $outdir/Group.txt grouping=\"genotype,allele,tissue,pna\"";
  System("$group_cmd >> $outdir/group_viz/R.out 2>&1");
  my $shm_group_cmd = "Rscript $FindBin::Bin/../R/mutationStats.R $shmexptfile $shmclonefile $meta_file $outdir/GroupSHM.txt grouping=\"genotype,allele,tissue,pna\"";
  System("$shm_group_cmd >> $outdir/group_viz/R.out 2>&1");
  my $del_shm_group_cmd = "Rscript $FindBin::Bin/../R/mutationStats.R $delshmexptfile $delshmclonefile $meta_file $outdir/GroupSHMDel.txt grouping=\"genotype,allele,tissue,pna\"";
  System("$del_shm_group_cmd >> $outdir/group_viz/R.out 2>&1");

  my $bx_group_cmd = "Rscript $FindBin::Bin/../R/baseExGroup.R $bxexptfile $meta_file $outdir/GroupBXTabular.txt";
  System("$bx_group_cmd");

  
  my $infh = IO::File->new("<$outdir/GroupBXTabular.txt");
  my $outfh = IO::File->new(">$outdir/GroupBX.txt");

  my $csv = Text::CSV->new({sep_char => "\t"});
  my $header = $csv->getline($infh);
  $csv->column_names(@$header);


  while (my $group = $csv->getline_hr($infh)) {
    $outfh->print(join(" ",$group->{genotype},$group->{allele},$group->{tissue},$group->{pna})."\n");
    $outfh->print(join("\t","Base",@bases)."\n");
    foreach my $b1 (@bases) {
      my @row = ($b1);
      foreach my $b2 (@bases) {
        if ($b1 eq $b2) {
          push(@row, "-");
        } else {
          push(@row,$group->{$b1."..".$b2})
        }
      }
      $outfh->print(join("\t",@row)."\n");
    }
    $outfh->print("\n");

  }

  my $mouse_cmd = "Rscript $FindBin::Bin/../R/mutationStats.R $exptfile $clonefile $meta_file $outdir/Mouse.txt grouping=\"genotype,allele,mouse,tissue,pna\"";
  System("$mouse_cmd >> $outdir/group_viz/R.out 2>&1");
  my $shm_mouse_cmd = "Rscript $FindBin::Bin/../R/mutationStats.R $shmexptfile $shmclonefile $meta_file $outdir/MouseSHM.txt grouping=\"genotype,allele,mouse,tissue,pna\"";
  System("$shm_mouse_cmd >> $outdir/group_viz/R.out 2>&1");
  my $del_shm_mouse_cmd = "Rscript $FindBin::Bin/../R/mutationStats.R $delshmexptfile $delshmclonefile $meta_file $outdir/MouseSHMDel.txt grouping=\"genotype,allele,mouse,tissue,pna\"";
  System("$del_shm_mouse_cmd >> $outdir/group_viz/R.out 2>&1");
}

sub read_in_meta_file {
	System("perl -pi -e 's/\\r/\\n/g' $meta_file");

	print "\nReading in meta file...\n";

	my $meta = IO::File->new("<$meta_file");
	my $csv = Text::CSV->new({sep_char => "\t"});
	my $header = $csv->getline($meta);
	$csv->column_names(@$header);

	while (my $row = $csv->getline_hr($meta)) {
    next unless $row->{experiment} =~ /\S/;
		$meta_hash{$row->{experiment}} = $row;
	}

}

sub parse_command_line {
	my $help;

	usage() if (scalar @ARGV == 0);

	my $result = GetOptions ( 
														"meta=s" => \$meta_file ,
														"in=s" => \$indir ,
														"out=s" => \$outdir ,
                            "ref=s" => \$refdir ,
                            "minqual=i" => \$min_qual ,
														"threads=i" => \$max_threads ,
                            "phred" => \$phred ,
														"ow" => \$ow ,
														"help" => \$help

				            			);
	
	usage() if ($help);


  #Check options

  croak "Error: must specifiy --meta" unless (defined $meta_file);
  croak "Error: cannot find $meta_file" unless (-r $meta_file);
  croak "Error: must specify --in" unless (defined $indir);
  croak "Error: cannot find $indir" unless (-d $indir);
  croak "Error: must specify --out" unless (defined $outdir);
	unless (-d $outdir) {
  	mkdir $outdir or croak "Error: output directory $outdir does not exist and cannot be created";
  }
	exit unless $result;
}


sub usage()
{
print<<EOF;
SHMPipeline.pl, by Robin Meyers, 12aug2013


Usage: $0 --metafile FILE (--in FILE | --indir DIR) --outdir DIR
		[--bcmismatch N] [--blastopt "-opt val"] [--threads N] [--ow]

Arguments (defaults in parentheses):

$arg{"--meta","Tab-delimited file containing experiment information"}
$arg{"--in","Input directory"}
$arg{"--out","Output directory"}
$arg{"--ref","Reference directory"}
$arg{"--threads","Number of processing core threads to run on",$max_threads}
$arg{"--ow","Overwrite output files if they already exist"}
$arg{"--help","This helpful help screen."}


EOF

exit 1;
}
