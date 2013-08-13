#!/usr/bin/perl

use strict;
use warnings;
use Carp;

require "/Users/robin/Working/LSYeap/lib/mutationVDJSwitchHelper.pl";

sub listItemsInDir ($)
{
  my $path = shift;
  croak "Error: path $path not found" unless -d $path;
  my ($record, @data, @filenames);
  my @list = `ls -l $path`;

  foreach $record (@list) {
    @data = split(/\s+/,$record);
    next unless defined($data[8]);
    push(@filenames, join("/",$path,join(" ",@data[8..$#data])));
  }
  
  carp "Warning: nothing found in path $path" unless @filenames > 0;

  return @filenames;
}


my $dir = shift;
croak "No directory" unless defined $dir && -d $dir;

my @items = listItemsInDir($dir);
print(join("\n",@items)."\n");

foreach my $item (@items) {
  next unless $item =~ /\s/;

  (my $olditem = $item) =~ s/\s+/\\ /g;
  print "$olditem\n";
  (my $newitem = $item) =~ s/\s+/_/g;

  System("mv $olditem $newitem"); 
}