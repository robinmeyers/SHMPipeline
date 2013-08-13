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
use threads;
use threads::shared;
use Interpolation 'arg:@->$' => \&argument;
use Time::HiRes qw(gettimeofday tv_interval);

use Cwd qw(abs_path);
use FindBin;
use lib abs_path("$FindBin::Bin/../lib");

require "mutationVDJSwitchHelper.pl";
require "pslHelper.pl";


# Flush output after every write
select( (select(STDOUT), $| = 1 )[0] );



# Forward declarations
sub parse_command_line;
sub read_in_meta_file;
sub check_existance_of_files;
sub initialize_stats_hash;
sub process_experiment ($$);
sub write_summary_stats;
sub create_summary;


# Global flags and arguments, 
# Set by command line arguments
my $meta_file;
my $indir;
my $outdir;
my $refdir;

my $userscreenopt = "";
my $userphredopt = "";
my $userphrapopt = "";
my $userswopt = "";
my $max_threads = 2;
my $phred;
my $min_qual = 15;
my $min_sw_score = 100;
my $shm_threshold = 3;
my $ow;

# Global variabless
my %meta_hash;
my $defaultphredopt = "-trim_alt \"\" -trim_cutoff 0.05 -trim_fasta";
my $defaultscreenopt = "-minmatch 12 -minscore 20";
my $defaultphrapopt = "-minscore 100";
my $defaultswopt = "-gapopen 13 -gapextend 0.1 -supper1 -supper2 -datafile /usr/local/emboss/data/EDNACUST";
my $sw_outfmt = "-aformat markx10";
my $exptfile;
my $shmexptfile;
my $clonefile;
my $shmclonefile;
my $mutfile;
#
# Start of Program
#

parse_command_line;

my $t0 = [gettimeofday];


my $phredopt = manage_program_options($defaultphredopt,$userphredopt);
my $screenopt = manage_program_options($defaultscreenopt,$userscreenopt);
my $phrapopt = manage_program_options($defaultphrapopt,$userphrapopt);
my $swopt = manage_program_options($defaultswopt,$userswopt) . " $sw_outfmt";


read_in_meta_file;

check_existance_of_files;

#initialize_stats_hash;


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
                        
                        process_experiment($expt_id, $meta_hash{$expt_id} );
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

#write_summary_stats;

create_summary;

my $t1 = tv_interval($t0);

printf("\nFinished all processes in %.2f seconds.\n", $t1);


#
# End of program
#


sub create_summary {

  my %stats;
  my $exptfh = IO::File->new(">$exptfile");
  my $shmexptfh = IO::File->new(">$shmexptfile");
  my $clonestatsfh = IO::File->new(">$clonefile");
  my $shmclonestatsfh = IO::File->new(">$shmclonefile");
  $exptfh->print(join("\t",qw(Expt Clones Bp Subs Del DelBp Ins InsBp RefA RefC RefG RefT RefN))."\n");
  $shmexptfh->print(join("\t",qw(Expt Clones Bp Subs Del DelBp Ins InsBp RefA RefC RefG RefT RefN))."\n");
  $clonestatsfh->print(join("\t",qw(Expt Clone Bp Subs Del DelBp Ins InsBp RefA RefC RefG RefT RefN Coords))."\n");
  $shmclonestatsfh->print(join("\t",qw(Expt Clone Bp Subs Del DelBp Ins InsBp RefA RefC RefG RefT RefN Coords))."\n");


  foreach my $expt (sort keys %meta_hash) {
    my $clonefh = IO::File->new("<".$meta_hash{$expt}->{clonefile});
    my $csv = Text::CSV->new({sep_char => "\t"});
    my $header = $csv->getline($clonefh);
    $csv->column_names(@$header);

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

    $stats{$expt}->{SHMClones} = 0;
    $stats{$expt}->{SHMBp} = 0;
    $stats{$expt}->{SHMSubs} = 0;
    $stats{$expt}->{SHMDel} = 0;
    $stats{$expt}->{SHMDelBp} = 0;
    $stats{$expt}->{SHMIns} = 0;
    $stats{$expt}->{SHMInsBp} = 0;
    $stats{$expt}->{SHMRefA} = 0;
    $stats{$expt}->{SHMRefC} = 0;
    $stats{$expt}->{SHMRefG} = 0;
    $stats{$expt}->{SHMRefT} = 0;
    $stats{$expt}->{SHMRefN} = 0;

    while (my $clone = $csv->getline_hr($clonefh)) {

      next unless $clone->{Bp} > 0;
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

      $clonestatsfh->print(join("\t",$expt,
                                $clone->{ID},
                                $clone->{Bp},
                                $clone->{Subs},
                                $clone->{Dels},
                                $clone->{DelBp},
                                $clone->{Ins},
                                $clone->{InsBp},
                                $clone->{RefA},
                                $clone->{RefC},
                                $clone->{RefG},
                                $clone->{RefT},
                                $clone->{RefN},
                                $clone->{Coords})."\n");

      next unless $clone->{Subs} > $shm_threshold;
      $stats{$expt}->{SHMClones}++;
      $stats{$expt}->{SHMBp} += $clone->{Bp};
      $stats{$expt}->{SHMSubs} += $clone->{Subs};
      $stats{$expt}->{SHMDel} += $clone->{Dels};
      $stats{$expt}->{SHMDelBp} += $clone->{DelBp};
      $stats{$expt}->{SHMIns} += $clone->{Ins};
      $stats{$expt}->{SHMInsBp} += $clone->{InsBp};
      $stats{$expt}->{SHMRefA} += $clone->{RefA};
      $stats{$expt}->{SHMRefC} += $clone->{RefC};
      $stats{$expt}->{SHMRefG} += $clone->{RefG};
      $stats{$expt}->{SHMRefT} += $clone->{RefT};
      $stats{$expt}->{SHMRefN} += $clone->{RefN};

      $shmclonestatsfh->print(join("\t",$expt,
                                $clone->{ID},
                                $clone->{Bp},
                                $clone->{Subs},
                                $clone->{Dels},
                                $clone->{DelBp},
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

    my $mutrate = $stats{$expt}->{Bp} - $stats{$expt}->{DelBp} > 0 ? $stats{$expt}->{Subs}/($stats{$expt}->{Bp} - $stats{$expt}->{DelBp}) : "";
    my $shmmutrate = $stats{$expt}->{SHMBp} - $stats{$expt}->{SHMDelBp} > 0 ? $stats{$expt}->{SHMSubs}/($stats{$expt}->{SHMBp} - $stats{$expt}->{SHMDelBp}) : "";

    $exptfh->print(join("\t",$expt,$stats{$expt}->{Clones},
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

    $shmexptfh->print(join("\t",$expt,$stats{$expt}->{SHMClones},
                                  $stats{$expt}->{SHMBp},
                                  $stats{$expt}->{SHMSubs},
                                  $stats{$expt}->{SHMDel},                                  
                                  $stats{$expt}->{SHMDelBp},
                                  $stats{$expt}->{SHMIns},
                                  $stats{$expt}->{SHMInsBp},                                  
                                  $stats{$expt}->{SHMRefA},
                                  $stats{$expt}->{SHMRefC},
                                  $stats{$expt}->{SHMRefG},
                                  $stats{$expt}->{SHMRefT},
                                  $stats{$expt}->{SHMRefN})."\n");

    
  }

  $exptfh->close;
  $shmexptfh->close;
  $clonestatsfh->close;
  $shmclonestatsfh->close;
  System("cat $outdir/*/*_muts.txt > $mutfile");

  mkdir "$outdir/group_viz";
  my $viz_cmd = "Rscript $FindBin::Bin/../R/mutationVizGrouped.R $meta_file $mutfile $clonefile $refdir $outdir/group_viz/";
  System("$viz_cmd > $outdir/group_viz/R.out 2>&1");
}

sub process_experiment ($$) {

	my $expt_id = shift;
	my $expt_hash = shift;


  if ($phred || defined $expt_hash->{raw}) {
    my $phred_cmd = join(" ","phred -id",$expt_hash->{raw},"-sa",$expt_hash->{phred},"-qa",$expt_hash->{phred}.".qual",$phredopt,">>",$expt_hash->{exptdir}."/phred.out 2>&1");
    System($phred_cmd);
  }

  my $screen_cmd = join(" ","cross_match",$expt_hash->{phred},$expt_hash->{vector},$screenopt,"-screen",">>",$expt_hash->{exptdir}."/screen.out 2>&1");
  System($screen_cmd);

  $expt_hash->{seqdir} = $expt_hash->{exptdir}."/sequences";

  mkdir $expt_hash->{seqdir} unless -d $expt_hash->{seqdir};


  $expt_hash->{phred} =~ s/\/\//\//g;
  my ($phr_path,$phr_name,$phr_ext ) = parseFilename($expt_hash->{phred});
  System(join(" ","mv",$expt_hash->{phred}.".screen",$expt_hash->{seqdir}."/".$phr_name.$phr_ext),1);
  $expt_hash->{phred} = $expt_hash->{seqdir}."/".$phr_name.$phr_ext;

  $expt_hash->{qual} =~ s/\/\//\//g;
  my ($qual_path,$qual_name,$qual_ext) = parseFilename($expt_hash->{qual});
  System(join(" ","cp",$expt_hash->{qual},$expt_hash->{seqdir}."/".$qual_name.$qual_ext),1);
  $expt_hash->{qual} = $expt_hash->{seqdir}."/".$qual_name.$qual_ext;


  $expt_hash->{phrapdir} = $expt_hash->{exptdir}."/phrap";
  mkdir $expt_hash->{phrapdir} unless -d $expt_hash->{phrapdir};

  phrap_sequence_pairs($expt_id,$expt_hash,$phrapopt);

  $expt_hash->{water} = $expt_hash->{exptdir} . "/${expt_id}.water";
  $expt_hash->{waternice} = $expt_hash->{exptdir} . "/${expt_id}_nice.water";

  (my $niceswopt = $swopt) =~ s/-aformat markx10//;
  System(join(" ","water","-asequence",$expt_hash->{reference},"-bsequence",$expt_hash->{contig},"-outfile",$expt_hash->{water},$swopt));
  System(join(" ","water","-asequence",$expt_hash->{reference},"-bsequence",$expt_hash->{contig},"-outfile",$expt_hash->{waternice},$niceswopt));

  parse_smith_water_file($expt_id,$expt_hash,$min_sw_score,$min_qual);

  System(join(" ","Rscript $FindBin::Bin/../R/removeDupClones.R",$expt_hash->{mutfile},$expt_hash->{clonefile}));

  my $vizdir = $outdir."/viz";
  mkdir $vizdir;
  System(join(" ","Rscript $FindBin::Bin/../R/mutationVizSuite.R",$expt_hash->{mutfile},$expt_hash->{clonefile},$expt_hash->{reference},"$vizdir/$expt_id",">",$expt_hash->{exptdir}."/R.out 2>&1"));


}

# sub initialize_stats_hash {
#  	foreach my $expt_id (sort keys %meta_hash) {
#  		my %temp_hash :shared;
#     $temp_hash{bases} = 0;
#  		$temp_hash{subs} = 0;
#  		$temp_hash{ins} = 0;
#  		$temp_hash{del} = 0;
#     $temp_hash{insbp} = 0;
#     $temp_hash{delbp} = 0;
#  		$stats{$expt_id} = \%temp_hash;
#  	}
# }


sub check_existance_of_files {
	print "\nSearching for sequence files...\n";


	foreach my $expt_id (sort keys %meta_hash) {
		my $input = $indir."/".$expt_id;
		my @exts = qw(.fa .fasta);
		foreach my $ext (@exts) {
			if (-r $input.$ext && -r $input.$ext.".qual") {
				$meta_hash{$expt_id}->{phred} = $input.$ext;
        $meta_hash{$expt_id}->{qual} = $input.$ext.".qual";
				last;
			}
		}
    if (! defined $meta_hash{$expt_id}->{phred} || $phred) {
      if (-d $input && scalar listFilesInDir($input) > 0 ) {
        $meta_hash{$expt_id}->{raw} = $input;
        $meta_hash{$expt_id}->{phred} = $input.".fa";
        $meta_hash{$expt_id}->{qual} = $meta_hash{$expt_id}->{phred}.".qual";

      } else {
        croak "Error: No valid experiment $expt_id in $indir" ;
      }
    }
	}
	

	print "\nSearching for reference files...\n";

	foreach my $expt_id (sort keys %meta_hash) {
		my $reffile = $refdir."/".$meta_hash{$expt_id}->{reference};
		croak "Error: Could not locate reference file $reffile in $indir" unless (-r $reffile);
		$meta_hash{$expt_id}->{reference} = $reffile;
    my $vecfile = $refdir."/".$meta_hash{$expt_id}->{vector};
    croak "Error: Could not locate reference file $vecfile in $indir" unless (-r $vecfile);
    $meta_hash{$expt_id}->{vector} = $vecfile;
	}
	
	print "Done.\n";

	print "\nChecking on output files...\n";
	foreach my $expt_id (sort keys %meta_hash) {
    my $exptdir = $outdir."/$expt_id/";
    $meta_hash{$expt_id}->{exptdir} = $outdir."/$expt_id";
    mkdir $meta_hash{$expt_id}->{exptdir};
		my $file = $exptdir."/".$expt_id."_muts.txt";
    my $stats = $exptdir."/".$expt_id."_clones.txt";
    my $base_ex = $exptdir."/".$expt_id."_baseEx.txt";
		croak "Error: Output file $file already exists and overwrite switch --ow not set" if (-r $file && ! defined $ow);
		System("touch $file",1);
		croak "Error: Cannot write to $file" unless (-w $file);
		$meta_hash{$expt_id}->{mutfile} = $file;
    $meta_hash{$expt_id}->{clonefile} = $stats;
    $meta_hash{$expt_id}->{base_ex} = $base_ex;

	}
	$exptfile = "$outdir/ExptStats.txt";
  $shmexptfile = "$outdir/SHMExptStats.txt";
  $clonefile = "$outdir/CloneStats.txt";
  $shmclonefile = "$outdir/SHMCloneStats.txt";
  $mutfile = "$outdir/MutFile.txt";
	croak "Error: Output file $exptfile already exists and overwrite switch --ow not set" if (-r $exptfile && ! defined $ow);
	System("touch $exptfile",1);
	croak "Error: Cannot write to $exptfile" unless (-w $exptfile);
	print "Done.\n";

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
		my $expt_id = $row->{experiment};
		$meta_hash{$expt_id} = $row;
	}

}

sub write_summary_stats {
#	print "\nWriting summary statistics...";
#	my $sumfh = IO::File->new(">$summaryfile");
#	$sumfh->print(join("\t",qw(Experiment Deletions Insertions SingleBpDel SingleBpIns Complex Total))."\n");
#	foreach my $expt_id (sort keys %stats) {
#		$sumfh->print(join("\t",$expt_id,
#														$stats{$expt_id}->{dels},
#														$stats{$expt_id}->{ins},
#														$stats{$expt_id}->{singledel},
#														$stats{$expt_id}->{singleins},
#														$stats{$expt_id}->{complex},
#														$stats{$expt_id}->{total})."\n");
#	}
#	$sumfh->close;
#	print "Done\n";
#

}

sub parse_command_line {
	my $help;

	usage() if (scalar @ARGV == 0);

	my $result = GetOptions ( 
														"meta=s" => \$meta_file ,
														"in=s" => \$indir ,
														"out=s" => \$outdir ,
                            "ref=s" => \$refdir ,
                            "screenopt=s" => \$userscreenopt ,
														"phredopt=s" => \$userphredopt ,
                            "phrapopt=s" => \$userphrapopt ,
                            "swopt=s" => \$userswopt ,
                            "minscore=i" => \$min_sw_score ,                            
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
  croak "Error: changing -aformat is not allowed" if ($userswopt =~ /aformat/);
	exit unless $result;
}


sub usage()
{
print<<EOF;
OttDeBruin.pl, by Robin Meyers, 05mar2013

This program analyzes insertions and deletions in NGS amplicon sequence data.

Usage: $0 --metafile FILE (--in FILE | --indir DIR) --outdir DIR
		[--bcmismatch N] [--blastopt "-opt val"] [--threads N] [--ow]

Arguments (defaults in parentheses):

$arg{"--meta","Tab-delimited file containing experiment information"}
$arg{"--in","Input directory"}
$arg{"--out","Output directory"}
$arg{"--ref","Reference directory"}
$arg{"--screenopt","Specify cross_match options for vector screening",$defaultscreenopt}
$arg{"--threads","Number of processing core threads to run on",$max_threads}
$arg{"--phred","Force phred to run"}
$arg{"--minscore","Minimum score to accept alignment",$min_sw_score}
$arg{"--minqual","Minimum quality to accept mutation",$min_qual}
$arg{"--ow","Overwrite output files if they already exist"}
$arg{"--help","This helpful help screen."}

Meta file format
----------------
The meta file must be a simple tab-delimited text file with the following headers:

experiment mid reference bindstart bindend cutstart cutend

Here is the text from a sample meta file ODB_meta.txt:

experiment         mid         reference  bindstart  bindend  cutstart  cutend
TALEN2_SCID058-IT  AGCACTGTAG  TALEN2.fa  74         123      92        107
TALEN2_neg_ctrl    ATCAGACACG  TALEN2.fa  74         123      92        107
TALEN3_SCID058-IT  TCGTCGCTCG  TALEN3.fa  206        255      223       238
TALEN3_neg_ctrl    ACATACGCGT  TALEN3.fa  206        255      223       238
TALEN4_SCID058-IT  CATAGTAGTG  TALEN4.fa  237        287      255       270
TALEN4_neg_ctrl    ATACGACGTA  TALEN4.fa  237        287      255       270


Program Description
-------------------

Example Commands
----------------

1)
OttDeBruin.pl --metafile ./in/ODB_meta.txt --in ./in/1.TCA.454ReadsLisa1_14_12.fna --outdir ./out/

EOF

exit 1;
}
