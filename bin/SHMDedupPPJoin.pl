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
use File::Basename;
use File::Spec;
use IO::File;
use Text::CSV;
use Switch;
use Data::Dumper;
use List::Util qw(min max);
use List::Compare;
use POSIX qw(ceil floor);
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


# Global flags and arguments, 
# Set by command line arguments
my $readsfile;
my $mutsfile;
my $j_thresh = 0.8;
my $tstart;
my $tend;


# Global variables
my $total_tokens = {};
my $read_tokens = {};
my $read_size = {};
my $inv_index = {};
my $dup_pairs = {};


#
# Start of Program
#

parse_command_line;

my $mutfh = IO::File->new("<$mutsfile");
my $mutcsv = Text::CSV->new({sep_char => "\t"});
$mutcsv->column_names(qw(Expt Read Pos Type From To Size End Ins));


my $readfh = IO::File->new("<$readsfile");
my $readcsv = Text::CSV->new({sep_char => "\t"});
$readcsv->column_names(qw(Expt Read Bp Coords Dup));

# Read in mutations and fill read tokens and total tokens hashes
while (my $mut = $mutcsv->getline_hr($mutfh)) {
  my @tokens;

  if (defined $tstart) {
    next if $mut->{Pos} < $tstart;
  }

  if (defined $tend) {
    next if $mut->{Pos} > $tend;
  }

  switch ($mut->{Type}) {
    case "sub" {
      @tokens = ($mut->{Pos}.$mut->{To});
    }
    case "del" {
      if (defined $tend) {
        next if $mut->{End} > $tend;
      }
      @tokens = ($mut->{Pos}."<",$mut->{End}.">");
    }
    case "ins" {
      @tokens = ($mut->{Pos}."i".$mut->{Ins});
    }
    else { carp "Warning: skipped line in mut file - could not determine mut type"; }
  }

  next unless @tokens > 0;

  if (exists $read_tokens->{$mut->{Read}}) {
    push($read_tokens->{$mut->{Read}},@tokens);
  } else {
    $read_tokens->{$mut->{Read}} = \@tokens;
  }

  foreach my $token (@tokens) {
    # if (exists $total_tokens->{$token}) {
      $total_tokens->{$token}++;
    # } else {
      # $total_tokens->{$token} = 1;
    # }
  }

  # if (exists $read_size->{$muts->{Read}}) {
    $read_size->{$mut->{Read}} += scalar @tokens;
  # } else {
    # $read_size->{$muts->{Read}} = scalar @tokens;
  # }

}

$mutfh->close;

foreach my $read (keys $read_tokens) {

  @{$read_tokens->{$read}} = sort {$total_tokens->{$a} <=> $total_tokens->{$b}} @{$read_tokens->{$read}};

}


# Apply algorithm

my $count = 0;

my @sorted_reads = sort {$read_size->{$a} <=> $read_size->{$b}} keys $read_size;
foreach my $x_id (@sorted_reads) {

  # if ($count++ % 1000 == 0) { print "$count\n";}


  next unless $read_size->{$x_id} > 0;

  my @x = @{$read_tokens->{$x_id}};
  my $x_prefix = @x - ceil($j_thresh * @x) + 1;
  my $overlap_map = {};

  for my $i (1..$x_prefix) {
    my $w = $x[$i-1];

    foreach my $index_lookup ( @{$inv_index->{$w}} ) {
      my ($y_id, $j) = @$index_lookup;
      my @y = @{$read_tokens->{$y_id}};
      next unless @y >= $j_thresh * @x;

      $overlap_map->{$y_id} = 0 unless exists $overlap_map->{$y_id};
      my $alpha = ceil($j_thresh/(1+$j_thresh)*(@x + @y));
      my $ubound = 1 + min(@x-$i,@y-$j);

      if ($overlap_map->{$y_id} + $ubound >= $alpha) {
        $overlap_map->{$y_id}++;
      } else {
        $overlap_map->{$y_id} = 0;
      }

    }

    # if (exists $inv_index->{$w}) {
      push($inv_index->{$w},[$x_id,$i]);
    # } else {
      # $inv_index->{$w} = [[$x_id,$i]];
    # }
  }

  foreach my $y_id (keys $overlap_map) {
    next unless $overlap_map->{$y_id} > 0;
    my @y = @{$read_tokens->{$y_id}};
    my $p_x = $x_prefix;
    my $p_y = @y - ceil($j_thresh * @y) + 1;
    my $w_x = $x[$p_x-1];
    my $w_y = $y[$p_y-1];
    my $overlap = $overlap_map->{$y_id};
    my $alpha = ceil($j_thresh/(1+$j_thresh)*(@x + @y));

    if ($total_tokens->{$w_x} < $total_tokens->{$w_y}) {
      my $ubound = $overlap + @x - $p_x;
      if ($ubound >= $alpha) {
        # print "x: @x\ny: @y\n";
        my $lc = List::Compare->new([@x[$p_x..$#x]],[@y[$overlap..$#y]]);
        $overlap = $overlap + $lc->get_intersection;
      }
    } else {
      my $ubound = $overlap + @y - $p_y;
      if ($ubound >= $alpha) {
        # print "x: @x\ny: @y\n";
        my $lc = List::Compare->new([@x[$overlap..$#x]],[@y[$p_y..$#y]]);
        $overlap = $overlap + $lc->get_intersection;
      }
    }
    if ($overlap >= $alpha) {
      $dup_pairs->{$y_id} = [] unless exists $dup_pairs->{$y_id};
      push($dup_pairs->{$y_id},$x_id);
    }
  }

}

my $tmpfile = $readsfile.".tmp";

rename $readsfile, $tmpfile;

open TMP, "<", $tmpfile;

open READS, ">", $readsfile;

while (<TMP>) {
  chomp;
  my @read = split("\t");
  if (exists $dup_pairs->{$read[1]}) {
    push(@read, $dup_pairs->{$read[1]}->[ $#{$dup_pairs->{$read[1]}} ]);
  }
  print READS join("\t",@read)."\n";
}

close TMP;
close READS;



#
# End of program
#


sub parse_command_line {
  my $help;

  usage() if (scalar @ARGV == 0);

  my $result = GetOptions ( 
                            "thresh" => \$j_thresh ,
                            "start" => \$tstart ,
                            "end" => \$tend ,
                            "help" => \$help

                          );
  
  usage() if ($help);


  $readsfile = shift(@ARGV);
  $mutsfile = shift(@ARGV);

  #Check options

  croak "Error: must specifiy reads file" unless (defined $readsfile);
  croak "Error: must specifiy muts file" unless (defined $mutsfile);
  exit unless $result;
}


sub usage()
{
print<<EOF;
SHMSanger.pl, by Robin Meyers, 05mar2013

This program analyzes somatic hypermutation from Sanger sequencing data.

Usage: $0 readsfile mutsfile reffile [--opts]

Arguments (defaults in parentheses):

$arg{"readsfile"," "}
$arg{"mutsfile"," "}
$arg{"--thresh","Jaccard threshold for calling a pair of reads duplicates",$j_thresh}
$arg{"--start","ignore mutations before this base on the reference"}
$arg{"--end","ignore mutations after this base on the reference"}

$arg{"--help","This helpful help screen."}


EOF

exit 1;
}