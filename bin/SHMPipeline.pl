#!/usr/bin/env perl

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
use Time::HiRes qw(gettimeofday tv_interval usleep);

use Cwd qw(abs_path);
use FindBin;
use lib abs_path("$FindBin::Bin/../lib");

require "SHMHelper.pl";


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
my $min_qual = 20;
my $dup_threshold = 0.9;
my $bt2_rfg = "5,3";
my $bt2_rdg = "8,1";
my $bt2_mp = "6,2";
my $bt2_dpad = 1000;
my $fragment;
my $dovetail;
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


my $default_se_bowtie_opt = "--very-sensitive-local -N 1 --mp $bt2_mp --rfg $bt2_rfg --rdg $bt2_rdg --dpad $bt2_dpad -p $expt_threads --reorder -t";
my $default_pe_bowtie_opt = "--very-sensitive-local -N 1 --mp $bt2_mp --rfg $bt2_rfg --rdg $bt2_rdg --dpad $bt2_dpad -X 2000 --no-discordant --no-mixed -p $expt_threads --reorder -t";
my $default_merge_bowtie_opt = "--very-sensitive -N 1 --score-min L,0,-1 --np 0 --mp $bt2_mp --rfg $bt2_rfg --rdg $bt2_rdg --dpad $bt2_dpad -p $expt_threads --reorder -t";

my $t0 = [gettimeofday];

my $se_bt2_opt = manage_program_options($default_se_bowtie_opt,$user_bowtie_opt);
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
            usleep(100000);
            last;
        }
        usleep(100000);
    } 
}

# waits for all threads to finish
while( scalar threads->list(threads::all) > 0) {
  for my $thr (@threads) {
      $thr->join() if $thr->is_joinable;
  }
  usleep(100000);
}

create_summary;

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


  print Capture(join(" ","Rscript $FindBin::Bin/../R/SHMDedup.R",
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
      $readend = $Rpos-1;
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


  my $merged_fq = $expt->{exptdir}."/".$expt->{experiment}."_merged.fq";
  $expt->{merged_fq} = $merged_fq;

  my $merged_fq_fh = Bio::SeqIO->new(-file => ">$merged_fq", -format => 'fastq');

  if (defined $expt->{R1} && defined $expt->{R2}) {

    my $align_stream = $samobj->get_seq_stream(-type => 'read_pair');
    while (my $pair = $align_stream->next_seq) {
      my ($Aln1,$Aln2) = $pair->get_SeqFeatures;

      ($Aln2,$Aln1) = ($Aln1,$Aln2) if $Aln1->reversed;


      # Assert that both pairs must map
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


      # Assert that the alignments extend to expected start and end of reference
      next if !defined $fragment && ($Start1 > $expt->{start} || $End2 < $expt->{end});




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

      if ($Start2 < $Start1) {
        carp "Warning: R2 start is less than R1 start";
        my $tmp_Start2 = $Start2;
        while ($tmp_Start2 < $Start1) {
          my $c2 = shift(@Cigar2);
          switch ($c2) {
            case 'M' {
              $Qpos2++;
              $tmp_Start2++;
            }
            case 'D' {
              $tmp_Start2++;
            }
            case 'I' {
              $Qpos2++;
            }
          }
        }
      }

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

      while ($Rpos <= $Aln1->end && $Rpos >= $Aln2->start && $Rpos <= $Aln2->end) {
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
  
  } else {
    
    my $align_stream = $samobj->get_seq_stream;

    while (my $Aln = $align_stream->next_seq) {

      next if $Aln->unmapped;
      next if !defined $fragment && ($Aln->start > $expt->{start} || $Aln->end < $expt->{end});

      my $Qseq = $Aln->query->seq->seq;
      my @Qual = $Aln->qscore;

      $Qseq = substr($Qseq,$Aln->query->start - 1, $Aln->query->end - $Aln->query->start + 1);
      @Qual = @Qual[($Aln->query->start - 1)..($Aln->query->end - 1)];
      
      if ($Aln->reversed) {
        $Qseq = reverseComplement($Qseq);
        @Qual =  reverse @Qual;
      }

      my $fq = Bio::Seq::Quality->new(-id => $Aln->qname, -seq => $Qseq, -qual => \@Qual);

      $merged_fq_fh->write_seq($fq);


    }
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

  my $bt2_cmd;

  if (defined $R1_fa && defined $R2_fa) {
    $pe_bt2_opt .= " --dovetail" if defined $dovetail;
    $bt2_cmd = "bowtie2 $pe_bt2_opt -x $ref_fa -1 $R1_fa -2 $R2_fa -S $pe_sam 2>&1";
  } elsif (defined $R1_fa) {
    $bt2_cmd = "bowtie2 $se_bt2_opt -x $ref_fa -U $R1_fa -S $pe_sam 2>&1";
  } elsif (defined $R2_fa) {
    $bt2_cmd = "bowtie2 $se_bt2_opt -x $ref_fa -U $R2_fa -S $pe_sam 2>&1";
  }

  print Capture($bt2_cmd);

  print Capture("samtools view -bSh -o $pe_bam $pe_sam 2>&1");
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

  print Capture("samtools view -bSh -o $bam $sam 2>&1");
}

sub check_existance_of_files {
	print "\nSearching for sequence files...\n";


	foreach my $expt_id (sort keys %meta_hash) {
		my $input = $indir."/".$expt_id;

    for my $ext (qw(fastq fq fastq.gz fq.gz)) {
      if (-r "${input}_R1.$ext" || -r "${input}_R2.$ext") {
        $meta_hash{$expt_id}->{R1} = "${input}_R1.$ext" if -r "${input}_R1.$ext";
        $meta_hash{$expt_id}->{R2} = "${input}_R2.$ext" if -r "${input}_R2.$ext";
        last;
      } elsif (-r "$input.$ext") {
        $meta_hash{$expt_id}->{R1} = "$input.$ext";
        last;
      }
    }


    croak "Error: No valid experiment $expt_id in $indir"
      unless defined $meta_hash{$expt_id}->{R1} || defined $meta_hash{$expt_id}->{R2};  
    
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


  croak "Error: Output file $exptfile already exists and overwrite switch --ow not set" if (-r $exptfile && ! defined $ow);
  System("touch $exptfile",1);
  croak "Error: Cannot write to $exptfile" unless (-w $exptfile);
  print "Done.\n";

}

sub create_summary {


  System(join(" ","Rscript $FindBin::Bin/../R/SHMAggregate.R",
                   $meta_file,$outdir));

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
                            "fragment" => \$fragment,
                            "dovetail" => \$dovetail,
                            "bt2-mp=s" => \$bt2_mp ,
                            "bt2-rfg=s" => \$bt2_rfg ,
                            "bt2-rdg=s" => \$bt2_rdg ,
                            "bt2-dpad=i" => \$bt2_dpad ,
                            "min-qual=i" => \$min_qual ,
                            "dup-thresh=f" =>\$dup_threshold ,
														"threads=i" => \$max_threads ,
                            "expt-threads=i" => \$expt_threads ,
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
$arg{"--bt2-mp","Maximum and minimum mismatch penalty, separated by a comma",$bt2_mp}
$arg{"--bt2-rfg","Reference gap open and extend penalty, separated by a comma",$bt2_rfg}
$arg{"--bt2-rfg","Read gap open and extend penalty, separated by a comma",$bt2_rdg}
$arg{"--bt2-dpad","Dynamic programming pad to allow gaps",$bt2_dpad}
$arg{"--min-qual","Minimum quality score for a base to count as a mutation",$min_qual}
$arg{"--dup-thresh","Minimum similarity fraction (0-1) for two reads to be called dups",$dup_threshold}
$arg{"--threads","Number of experiments to run simultaneously",$max_threads}
$arg{"--expt-threads","Number of cores to use for each experiment ",$expt_threads}
$arg{"--ow","Overwrite output files if they already exist"}
$arg{"--help","This helpful help screen."}


EOF

exit 1;
}
