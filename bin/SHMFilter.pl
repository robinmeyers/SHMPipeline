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
use Interpolation 'arg:@->$' => \&argument;

use Cwd qw(abs_path);
use FindBin;
use lib abs_path("$FindBin::Bin/../lib");

require "SHMHelper.pl";

# Flush output after every write
select( (select(STDOUT), $| = 1 )[0] );


# Forward declarations
sub parse_command_line;
sub align_sequences;


# Global flags and arguments, 
# Set by command line arguments
my $input_prefix;
my $output_prefix;
my $filter_seq;
my $infq;
my $inreads;
my $inmuts;
my $sam;
my $bam;
my $outfq;
my $outreads;
my $outmuts;
my $min_score = 40;
my $filtered_reads = 0;
my $total_reads = 0;
my $filtered_muts = 0;
my $total_muts = 0;

# Global variabless
my %badids;

#
# Start of Program
#

parse_command_line;

align_sequences;

my $samobj = Bio::DB::Sam->new(-bam => $bam,
                               -fasta => $filter_seq,
                               -expand_flags => 1);
my $align_stream = $samobj->get_seq_stream(-type => 'match');

while (my $aln = $align_stream->next_seq) {
  $badids{$aln->query->name} = 1;
}

my $infqfh = Bio::SeqIO->new(-file => $infq, -format => 'fastq');
my $outfqfh = Bio::SeqIO->new(-file => ">$outfq", -format => 'fastq');

while (my $fq = $infqfh->next_seq) {
  $outfqfh->write_seq($fq) unless exists $badids{$fq->id};
}

my $inreadfh = IO::File->new($inreads);
my $readcsv = Text::CSV->new({sep_char => "\t"});
my $outreadfh = IO::File->new(">$outreads");

while (my $read = $readcsv->getline($inreadfh)) {
  $outreadfh->print(join("\t",@$read)."\n") unless exists $badids{$read->[1]};
  $total_reads++;
}

my $inmutfh = IO::File->new($inmuts);
my $mutcsv = Text::CSV->new({sep_char => "\t"});
my $outmutfh = IO::File->new(">$outmuts");

while (my $mut = $mutcsv->getline($inmutfh)) {
  $outmutfh->print(join("\t",@$mut)."\n") unless exists $badids{$mut->[1]};
  $filtered_muts++ if exists $badids{$mut->[1]};
  $total_muts++;
}

$filtered_reads = scalar keys %badids;

print "Filtered $filtered_reads reads from $total_reads total\n";
print "Filtered $filtered_muts mutations from $total_muts total\n";

#
# End of program
#


sub align_sequences {

  $sam = "${output_prefix}.sam";
  $bam = "${output_prefix}.bam";

  System("bowtie2-build -q $filter_seq $filter_seq") unless -r "$filter_seq.1.bt2";

  System("bowtie2 --no-unal --very-sensitive-local -L 10 -i C,1 --min-score C,$min_score -x $filter_seq -U $infq -S $sam");

  System("samtools view -bSh -o $bam $sam");

}

sub parse_command_line {
  my $help;

  usage() if (scalar @ARGV == 0);

  my $result = GetOptions ( 
                            "in=s" => \$input_prefix ,
                            "out=s" => \$output_prefix ,
                            "seq=s" => \$filter_seq ,
                            "min-score=i" => \$min_score ,
                            "help" => \$help

                          );
  
  usage() if ($help);


  #Check options

  $infq = $input_prefix."_merged.fq";
  $inreads = $input_prefix."_reads.txt";
  $inmuts = $input_prefix."_muts.txt";

  $outfq = $output_prefix."_merged.fq";
  $outreads = $output_prefix."_reads.txt";
  $outmuts = $output_prefix."_muts.txt";

  croak "Error: must specify --in" unless (defined $input_prefix);
  croak "Error: must specify --out" unless (defined $output_prefix);
  croak "Error: cannot find sequence file to filter against" unless (-r $filter_seq);

  croak "Error: input prefix and output prefix may not be identical" if ($input_prefix eq $output_prefix);

  croak "Error: cannot find fastq file $infq" unless (-r $infq);
  croak "Error: cannot find reads file $inreads" unless (-r $inreads);
  croak "Error: cannot find muts file $inmuts" unless (-r $inmuts);

  exit unless $result;
}


sub usage()
{
print<<EOF;
SHMPipeline.pl, by Robin Meyers, 12aug2013


Usage: $0 --in PREFIX --out PREFIX --seq FASTA
    [--min-score N] [--help]

Arguments (defaults in parentheses):

$arg{"--in","Input prefix"}
$arg{"--out","Output prefix"}
$arg{"--seq","Fasta file to filter sequences against"}
$arg{"--min-score","Minimum alignment score to filter a read",$min_score}
$arg{"--help","This helpful help screen."}


EOF

exit 1;
}