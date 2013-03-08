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

require "OttDeBruinHelper.pl";

# Flush output after every write
select( (select(STDOUT), $| = 1 )[0] );



# Forward declarations
sub parse_command_line;
sub read_in_meta_file;
sub check_existance_of_files;
sub demultiplex_reads;
sub initialize_stats_hash;
sub process_experiment ($$);
sub write_summary_stats;


# Global flags and arguments, 
# Set by command line arguments
my $meta_file;
my $infile;
my $indir;
my $outdir;
my $summaryfile;
my $userblastopt = "";
my $bc_mismatches = 1;
my $max_threads = 2;
my $ow;

# Global variabless
my %meta_hash;
my %stats :shared;
my $defaultblastopt = "-task blastn-short -strand plus -gapopen 1 -gapextend 1 -reward 1 -penalty -4 -max_target_seqs 1";
my @blast_header = qw(qseqid sseqid pident length mismatch gapopen qstart qend sstart send qseq sseq);
my $blast_outfmt = "-outfmt \"6 ".join(" ",@blast_header)."\"";

#
# Start of Program
#

parse_command_line;

($indir, my $filename, my $inext) = parseFilename($infile) if defined $infile;

my $t0 = [gettimeofday];

my $blastopt = manage_blast_options($defaultblastopt,$userblastopt) . " $blast_outfmt";

read_in_meta_file;

check_existance_of_files;

demultiplex_reads if (defined $infile);

initialize_stats_hash;


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
                        $meta_hash{$expt_id}->{exptdir} = $outdir;
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

write_summary_stats;

my $t1 = tv_interval($t0);

printf("\nFinished all processes in %.2f seconds.\n", $t1);


#
# End of program
#

sub process_experiment ($$) {

	my $expt_id = shift;
	my $expt_hash = shift;

	blast_to_reference_file($expt_id,$expt_hash,$blastopt);

	parse_blast_output($expt_id,$expt_hash,$stats{$expt_id},\@blast_header);

}

sub initialize_stats_hash {
 	foreach my $expt_id (sort keys %meta_hash) {
 		my %temp_hash :shared;
 		$temp_hash{total} = 0;
 		$temp_hash{ins} = 0;
 		$temp_hash{dels} = 0;
 		$temp_hash{singledel} = 0;
 		$temp_hash {singleins} = 0;
 		$temp_hash{complex} = 0;
 		$stats{$expt_id} = \%temp_hash;
 	}
}


sub demultiplex_reads {

	# create barcode file
	my $bcfile = "$indir/bcfile.txt";
	my $bcfh = IO::File->new(">$bcfile") or croak "Error: cannot open $bcfile for writing";
	foreach my $expt_id (sort keys %meta_hash) {
		$bcfh->print(join("\t",$meta_hash{$expt_id}->{experiment},$meta_hash{$expt_id}->{mid})."\n");
		$meta_hash{$expt_id}->{raw} = "$indir/".$meta_hash{$expt_id}->{experiment}.".fa";

	}
	$bcfh->close;

	print "\nDemultiplexing reads...\n";

	my $cmd = "$FindBin::Bin/fasta_barcode_splitter.pl --bcfile $bcfile --fasta $infile --prefix $indir --suffix .fa --bol --mismatches $bc_mismatches";
	System($cmd);
	print "Done\n";
}

sub check_existance_of_files {
	print "\nSearching for sequence files...\n";

	unless (defined $infile) {
		foreach my $expt_id (sort keys %meta_hash) {
			my $file = $indir."/".$expt_id;
			my @exts = qw(.fa .fasta);
			foreach my $ext (@exts) {
				if (-r $file.$ext) {
					$meta_hash{$expt_id}->{raw} = $file.$ext;
					last;
				}
			}
			croak "Error: Could not locate sequence file $file in $indir - make sure it has .fa or .fasta extension" unless (defined $meta_hash{$expt_id}->{raw});
		}
	} else {
		croak "Error: could not locate input file $infile" unless (-r $infile);
	}

	print "\nSearching for reference files...\n";

	foreach my $expt_id (sort keys %meta_hash) {
		my $file = $indir."/".$meta_hash{$expt_id}->{reference};
		croak "Error: Could not locate reference file $file in $indir" unless (-r $file);
		$meta_hash{$expt_id}->{reference} = $file;
	}
	
	print "Done.\n";

	print "\nChecking on output files...\n";
	foreach my $expt_id (sort keys %meta_hash) {
		my $file = $outdir."/".$expt_id.".txt";
		croak "Error: Output file $file already exists and overwrite switch --ow not set" if (-r $file && ! defined $ow);
		System("touch $file");
		croak "Error: Cannot write to $file" unless (-w $file);
		$meta_hash{$expt_id}->{output} = $file;
	}
	$summaryfile = "$outdir/SummaryStats.txt";
	croak "Error: Output file $summaryfile already exists and overwrite switch --ow not set" if (-r $summaryfile && ! defined $ow);
	System("touch $summaryfile");
	croak "Error: Cannot write to $summaryfile" unless (-w $summaryfile);
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
		my $expt_id = $row->{experiment};
		$meta_hash{$expt_id} = $row;
	}

}

sub write_summary_stats {
	print "\nWriting summary statistics...";
	my $sumfh = IO::File->new(">$summaryfile");
	$sumfh->print(join("\t",qw(Experiment Deletions Insertions SingleBpDel SingleBpIns Complex Total))."\n");
	foreach my $expt_id (sort keys %stats) {
		$sumfh->print(join("\t",$expt_id,
														$stats{$expt_id}->{dels},
														$stats{$expt_id}->{ins},
														$stats{$expt_id}->{singledel},
														$stats{$expt_id}->{singleins},
														$stats{$expt_id}->{complex},
														$stats{$expt_id}->{total})."\n");
	}
	$sumfh->close;
	print "Done\n";


}

sub parse_command_line {
	my $help;

	usage() if (scalar @ARGV == 0);

	my $result = GetOptions ( 
														"metafile=s" => \$meta_file ,
														"in=s" => \$infile ,
														"indir=s" => \$indir,
														"outdir=s" => \$outdir ,
														"bcmismatch=i" => \$bc_mismatches ,
														"blastopt=s" => \$userblastopt ,
														"threads=i" => \$max_threads ,
														"ow" => \$ow ,
														"help" => \$help

				            			);
	
	usage() if ($help);


  #Check options

  croak "Error: must specifiy --metafile" unless (defined $meta_file);
  croak "Error: cannot find $meta_file" unless (-r $meta_file);
  croak "Error: must specifiy one and only one of --in and --indir"
  	unless ((defined $infile && ! defined $indir) || (defined $indir && ! defined $infile));
  croak "Error: cannot read $infile" if (defined $infile && ! -r $infile);
  croak "Error: cannot find $indir" if (defined $indir && ! -d $indir);
  croak "Error: must specify --outdir" unless (defined $outdir);
	unless (-d $outdir) {
  	mkdir $outdir or croak "Error: output directory $outdir does not exist and cannot be created";
  }
  croak "Error: changing -outfmt is not allowed" if ($userblastopt =~ /outfmt/);
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

$arg{"--metafile","Tab-delimited file containing experiment information"}
$arg{"--in","Fasta file with multiplexed reads from sequencer"}
$arg{"--indir","Directory containing demultiplexed fasta files"}
$arg{"--outdir","Output directory"}
$arg{"--bcmismatch","Number of mismatches allowed in matching barcode",$bc_mismatches}
$arg{"--blastopt","Specify Blast options - refer to 'blastn -help'",$defaultblastopt}
$arg{"--threads","Number of processing core threads to run on",$max_threads}
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
The multiplexed sequence file will be split up by MID, into experiment.fa files.
If --indir is given instead of --in, the program will look for experiment.fa files
in this directory, raising an error if any are not found.
The program will then look for the reference file (e.g. TALEN2.fa)
in the same folder as the sequence file(s).

The program then uses the blastn binary from the NCBI Blast+ toolkit to align
each sequence to the reference sequence.

Finally, the program parses this output, printing any insertion and deletion info
between the cut start and cut end coordinates of the reference sequence,
as well as the aligned sequences between the bind start and bind end coordinates.
This data goes into the experiment.txt file in the --outdir folder.

For example for the TALEN2_SCID058-IT experiment above, the output file will contain
aligned sequences between base 74 (inclusive) and 123 (exclusive), but will only count
indels that intersect the region between base 92 (inclusive) and 107 (exclusive).

The program prints out a summary file as well. The logic behind classifying each read
should be updated as desired in the parse_blast_output subroutine in lib/OttDeBruinHelper.pl

Example Commands
----------------

1)
OttDeBruin.pl --metafile ./in/ODB_meta.txt --in ./in/1.TCA.454ReadsLisa1_14_12.fna --outdir ./out/

Given this file tree
./ (working directory)
|-in/
|---ODB_meta.txt
|---1.TCA.454ReadsLisa1_14_12.fna
|---TALEN2.fa
|---TALEN3.fa
|---TALEN4.fa

This command will separate the input file into the six experiments above,
and run the default analysis.

This will be the resulting file tree
./ (working directory)
|-in/
|---ODB_meta.txt
|---1.TCA.454ReadsLisa1_14_12.fna
|---TALEN2.fa
|---TALEN3.fa
|---TALEN4.fa
|---TALEN2_SCID058-IT.fa
|---TALEN2_neg_ctrl.fa
|---TALEN3_SCID058-IT.fa
|---TALEN3_neg_ctrl.fa
|---TALEN4_SCID058-IT.fa
|---TALEN4_neg_ctrl.fa
|-out/
|---SummaryStats.txt
|---TALEN2_SCID058-IT.txt
|---TALEN2_neg_ctrl.txt
|---TALEN3_SCID058-IT.txt
|---TALEN3_neg_ctrl.txt
|---TALEN4_SCID058-IT.txt
|---TALEN4_neg_ctrl.txt
|---blasts/
|-----TALEN2_SCID058-IT_blast.txt
|-----TALEN2_neg_ctrl_blast.txt
|-----TALEN3_SCID058-IT_blast.txt
|-----TALEN3_neg_ctrl_blast.txt
|-----TALEN4_SCID058-IT_blast.txt
|-----TALEN4_neg_ctrl_blast.txt

2)
OttDeBruin.pl --metafile ./in/ODB_meta.txt --indir ./in/ --outdir ./out/
	--blastopt "-penalty -3" --threads 4 --ow

Given this file tree
./ (working directory)
|-in/
|---ODB_meta.txt
|---1.TCA.454ReadsLisa1_14_12.fna
|---TALEN2.fa
|---TALEN3.fa
|---TALEN4.fa
|---TALEN2_SCID058-IT.fa
|---TALEN2_neg_ctrl.fa
|---TALEN3_SCID058-IT.fa
|---TALEN3_neg_ctrl.fa
|---TALEN4_SCID058-IT.fa
|---TALEN4_neg_ctrl.fa

This command skips the demultiplexing steps and goes straight to analyzing
the reads. In contrast to the first example, this commmand will overwrite
any existing output files, changes the blast options to use a mismatch penalty
of -3, and runs on 4 threads instead of the default 2. The resulting file tree
will be the same as in example 1.
EOF

exit 1;
}
