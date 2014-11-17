#!/usr/bin/env perl

##
## This program analyzes somatic hypermutation
## from Sanger sequencing data
##
## run with "--help" for usage information
##
## Robin Meyers <robin.meyers@childrens.harvard.edu>, 05mar2013

use strict;
use warnings;
use Getopt::Long;
use Carp;
use File::Spec;
use IO::File;
use Text::CSV;
use List::MoreUtils qw(pairwise);
use Bio::SeqIO;
use Bio::Seq::Quality;
use Data::Dumper;
use threads;
use threads::shared;
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
sub process_experiment ($);



# Global flags and arguments, 
# Set by command line arguments
my $meta_file;
my $indir;


# Global variabless
my %meta;
my $phredopt = "-trim_alt \"\" -trim_cutoff 0.05 -trim_fasta";
my $defaultphredopt = "-trim_alt \"\" -trim_cutoff 0.05 -trim_fasta";
my $userphredopt;
my $max_threads = 4;

#
# Start of Program
#

parse_command_line;

my $t0 = [gettimeofday];

read_in_meta_file;


my @threads = ();


# Standard multithreading loop - each experiment is a thread
foreach my $expt_id (sort keys %meta) {

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
                        
                        process_experiment($meta{$expt_id} );
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

my $t1 = tv_interval($t0);

printf("\nFinished all processes in %.2f seconds.\n", $t1);


#
# End of program
#

sub process_experiment ($) {
  my $expt = shift;

  my $expt_path = File::Spec->catfile($indir,$expt->{experiment});

  my $logfile = $expt_path.".out";

  $expt->{logfh} = IO::File->new(">".$logfile);

  local *STDOUT = $expt->{logfh};
  local *STDERR = $expt->{logfh};
  select( (select(STDOUT), $| = 1 )[0] );
  select( (select(STDERR), $| = 1 )[0] );


  my $phred_cmd = join(" ", "phred -id", $expt_path,
                            "-sa", $expt_path.".fa",
                            "-qa", $expt_path.".fa.qual",
                            $phredopt,
                            ">>",$logfile,"2>&1");

  System($phred_cmd);

  my $fafh = Bio::SeqIO->new(-file => $expt_path.".fa", -format => 'fasta');
  my $qualfh = Bio::SeqIO->new(-file => $expt_path.".fa.qual", -format => 'qual');


  my $forward = $expt->{forward} if $expt->{forward} =~ /\S/;
  my $reverse = $expt->{reverse} if $expt->{reverse} =~ /\S/;

  if (defined $forward || defined $reverse) {
    my $seqs = {};

    while (my $seq = $fafh->next_seq) {

      my $qual = $qualfh->next_seq;

      if (defined $forward && $seq->id =~ /(.+)_$forward\./) {
        my $fq = Bio::Seq::Quality->new( -id => $1, -seq => $seq->seq, -qual => $qual->qual);
        
        $seqs->{$1} = {} unless defined $seqs->{$1};
        $seqs->{$1}->{R1} = $fq;
      
      } elsif (defined $reverse && $seq->id =~ /(.+)_$reverse\./) {
        my $fq = Bio::Seq::Quality->new( -id => $1, -seq => $seq->seq, -qual => $qual->qual);
        
        $seqs->{$1} = {} unless defined $seqs->{$1};
        $seqs->{$1}->{R2} = $fq;
      
      } else {
        next;
      }

    }

    my $R1 = Bio::SeqIO->new(-file => ">".$expt_path."_R1.fq", -format => 'fastq') if defined $forward;
    my $R2 = Bio::SeqIO->new(-file => ">".$expt_path."_R2.fq", -format => 'fastq') if defined $reverse;

    foreach my $id (keys %$seqs) {
      $R1->write_seq($seqs->{$id}->{R1}) if defined $seqs->{$id}->{R1};
      $R2->write_seq($seqs->{$id}->{R2}) if defined $seqs->{$id}->{R2};
    }


  } else {

    my $R1 = Bio::SeqIO->new(-file => ">".$expt_path.".fq", -format => 'fastq');

    while (my $seq = $fafh->next_seq) {
      my $qual = $qualfh->next_seq;
      my $fq = Bio::SeqIO::Quality->new( -id => $seq->id, -seq => $seq->seq, -qual => $qual->qual);
      $R1->write_seq($fq);
    }
  }

  unlink $expt_path.".fa";
  unlink $expt_path.".fa.qual";

}


sub read_in_meta_file {
  System("perl -pi -e 's/\\r/\\n/g' $meta_file");

  print "\nReading in meta file...\n";

  my $metafh = IO::File->new("<$meta_file");
  my $csv = Text::CSV->new({sep_char => "\t"});
  my $header = $csv->getline($metafh);
  $csv->column_names(@$header);

  while (my $row = $csv->getline_hr($metafh)) {
    next unless $row->{experiment} =~ /\S/;
    my $expt_id = $row->{experiment};
    $meta{$expt_id} = $row;
  }

}


sub parse_command_line {
  my $help;

  usage() if (scalar @ARGV == 0);

  my $result = GetOptions ( 
                            "meta=s" => \$meta_file ,
                            "in=s" => \$indir ,
                            "phredopt=s" => \$userphredopt ,
                            "threads=i" => \$max_threads ,
                            
                            "help" => \$help

                          );
  
  usage() if ($help);


  #Check options

  croak "Error: must specifiy --meta" unless (defined $meta_file);
  croak "Error: cannot find $meta_file" unless (-r $meta_file);
  croak "Error: must specify --in" unless (defined $indir);
  croak "Error: cannot find $indir" unless (-d $indir);
  exit unless $result;
}


sub usage()
{
print<<EOF;
SHMSanger.pl, by Robin Meyers, 05mar2013

This program analyzes somatic hypermutation from Sanger sequencing data.

Usage: $0 --metafile FILE (--in FILE | --indir DIR) --outdir DIR
    [--bcmismatch N] [--blastopt "-opt val"] [--threads N] [--ow]

Arguments (defaults in parentheses):

$arg{"--meta","Tab-delimited file containing experiment information"}
$arg{"--in","Input directory"}
$arg{"--threads","Number of processing core threads to run on",$max_threads}
$arg{"--help","This helpful help screen."}


EOF

exit 1;
}
