#!/usr/bin/perl

use strict;
use warnings;

my @array = (10,15,18,22,104);

my ($a) = grep { $array[$_] ~~ 10 } 0 .. $#array;
my ($b) = grep { $array[$_] ~~ 15 } 0 .. $#array;
my ($c) = grep { $array[$_] ~~ 16 } 0 .. $#array;
print "10: $a\n" if ($a);
print "15: $b\n" if ($b);
print "16: $c\n" if ($c);


my $num = 21;
@array = 10..20;
if ( $num ~~ @array ){
  print "$num: true\n";
}