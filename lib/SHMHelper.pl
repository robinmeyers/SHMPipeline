#!/usr/bin/perl

use strict;
use warnings;
use Carp;
use File::Basename;
use IO::File;
use Text::CSV;
use List::Util qw(min max sum);
use List::MoreUtils qw(pairwise);


sub argument {
	my $var = shift;
	my $description = shift;
	my $default = shift;

	return sprintf("  \%\-16s - %s\n",
		(defined $default ? sprintf("%s (%s)",$var,$default) : sprintf("%s",$var)),
		$description);
}

sub listFilesInDir ($)
{
	my $path = shift;
	croak "Error: path $path not found" unless -d $path;
	my ($record, @data, @filenames);
	my @list = `ls -l $path`;

	foreach $record (@list) {
		next if $record =~ /^d/;
		@data = split(/\s+/,$record);
		next unless defined($data[8]);
		push(@filenames, join("/",$path,$data[8]));
	}
	
	carp "Warning: nothing found in path $path" unless @filenames > 0;

	return @filenames;
}

sub parseFilename ($) {
	my $fullname = shift;
	my ($name, $path, $ext) = fileparse($fullname, qr/\.\w{2,5}$/);
	return ($path, $name, $ext);
}

#invert exit status of system call
sub System ($;$) {
	my $cmd = shift;
	my $quiet = shift;
	$quiet = 0 unless defined $quiet;
	print "$cmd\n" unless $quiet;
	my $status = system($cmd);
	return !$status;
}

sub Capture ($;$) {
  my $cmd = shift;
  my $quiet = shift;
  $quiet = 0 unless defined $quiet;
  print "$cmd\n" unless $quiet;
  my $output = capture($cmd);
  return $output; 
}

sub mean {
    return sum(@_)/@_;
}

sub read_fasta ($) {

	my $fh = shift;
	croak "Error: invalid file handle" unless $fh->opened;
	my $seq_name = $fh->getline();
	return () unless defined $seq_name;
	chomp($seq_name);
	croak "Error: unexpected format - missing '\>' symbol in sequence header $seq_name" unless $seq_name =~ s/^\>\s*(\S.*)$/$1/;
	my $seq_bases = "";
	while (my $line = $fh->getline()) {
		if ($line =~ /^>/) {
			seek($fh,-length($line),1);
			last;
		}
		chomp($line);
		croak "Error: unexpected base in sequence $seq_bases" unless $line =~ /^[AGCTagctNnXx\-\.]*$/;
		$seq_bases .= $line;
	}
	
	return ($seq_name,$seq_bases);
}

sub write_fasta ($$$) {
	croak "Error: too few arguments to 'write_fastq' method" if scalar @_ < 3;

	my $fh = shift;
	croak "Error: invalid file handle" unless $fh->opened;
	my $seq_name = shift;
	my $seq_bases = shift;
 

	$seq_name = ">" . $seq_name unless $seq_name =~ /^\>.*$/;
	croak "Error: sequence contains unexpected base" unless $seq_bases =~ /^[AGCTagctNnXx\-\.]+$/;

	print $fh $seq_name,"\n";
	print $fh $seq_bases,"\n";


}

sub read_qual ($) {

	my $fh = shift;
	croak "Error: invalid file handle" unless $fh->opened;
	my $seq_name = $fh->getline();
	return () unless defined $seq_name;
	chomp($seq_name);
	croak "Error: unexpected format - missing '\>' symbol in sequence header $seq_name" unless $seq_name =~ s/^\>\s*(\S.*)$/$1/;
	my @seq_quals = ();
	while (my $line = $fh->getline()) {
		if ($line =~ /^>/) {
			seek($fh,-length($line),1);
			last;
		}
		chomp($line);
		croak "Error: unexpected character in quality $line" unless $line =~ /^[\d\s]*$/;
		push(@seq_quals,split(/\s+/,$line));
	}
	
	return ($seq_name,\@seq_quals);
}

sub write_qual ($$$;$) {
	croak "Error: too few arguments to 'write_fastq' method" if scalar @_ < 3;

	my $fh = shift;
	croak "Error: invalid file handle" unless $fh->opened;
	my $seq_name = shift;
	my $quals = shift;
	my @seq_quals = @$quals;
	my $line_size = shift;
	$line_size = 20 unless defined $line_size;

	$seq_name = ">" . $seq_name unless $seq_name =~ /^\>.*$/;

	print $fh $seq_name,"\n";

	while (scalar @seq_quals > 0) {
		my $line = join(" ",splice(@seq_quals,0,$line_size));
		croak "Error: qualities contains unexpected character" unless $line =~ /^[\d\s]*$/;
		print $fh $line,"\n";
	}

}

sub reverseComplement ($) {
	my $seq = shift;
	(my $rc = reverse($seq)) =~ tr/ACGTacgtNn/TGCAtgcaNn/;
	return $rc;
}


sub manage_program_options ($$) {

  my $defaultoptstr = shift;
  my $useroptstr = shift;

  return $defaultoptstr unless $useroptstr =~ /\S/;
  return $useroptstr unless $defaultoptstr =~ /\S/;

  my @defaultopt = split(/\s+/,$defaultoptstr);
  my @useropt = split(/\s+/,$useroptstr);

  my %opt_hash;
  my $optkey;
  my $optval;

  my @return_opt = ();

  my $i = 0;

  while ($i < @defaultopt) {
    if ( defined $defaultopt[$i+1] && $defaultopt[$i+1] !~ /^-/ ) {
      $opt_hash{$defaultopt[$i]} = $defaultopt[$i+1];
      $i += 2;
    } else {
      $opt_hash{$defaultopt[$i]} = "";
      $i += 1;
    }
  }

  $i = 0;

  while ($i < scalar @useropt) {
    if ( defined $useropt[$i+1] && $useropt[$i+1] !~ /^-/ ) {
      if ($useropt[$i+1] eq "OFF") {
        delete $opt_hash{$useropt[$i]};
      } else {
        $opt_hash{$useropt[$i]} = $useropt[$i+1];
      }
      $i += 2;
    } else {
      $i += 1;
    }
  }
  
  foreach $optkey (sort keys %opt_hash) {
    $optval = $opt_hash{$optkey};
    push(@return_opt,$optkey,$optval);
  }

  return(join(" ",@return_opt));


}

sub manage_old_program_options ($$) {

	my $defaultoptstr = shift;
	my $useroptstr = shift;

	return $defaultoptstr unless $useroptstr =~ /\S/;
	return $useroptstr unless $defaultoptstr =~ /\S/;

	my @defaultopt = split(/\s+-/,$defaultoptstr);
	my @useropt = split(/\s+-/,$useroptstr);

	my %opt_hash;
	my $optkey;
	my $optval;

	my @return_opt = ();

	my $i = 1;
	while ($i < scalar @defaultopt) {
		if ($defaultopt[$i] =~ /^\d/) {
			$defaultopt[$i-1] .= " -".$defaultopt[$i];
			splice(@defaultopt,$i,1);
			next;
		}
		$i++;
	}
	$i = 1;
	while ($i < scalar @useropt) {
		if (defined $useropt[$i] && $useropt[$i] =~ /^\d/) {
			$useropt[$i-1] .= " -".$useropt[$i];
			splice(@useropt,$i,1);
			next;
		}
		$i++;
	}

	foreach (@defaultopt) {
		if ( /^-?(\S+)\s+(\S+)$/ ) {
			$opt_hash{$1} = $2;
		} elsif ( /^-?(\S+)$/ ) {
			$opt_hash{$1} = "";
		} else {
			croak "Error: incorrect default program option string $defaultoptstr";
		}
	}

	foreach (@useropt) {
		if ( /^-?(\S+)\s+(\S+)$/ ) {
			$opt_hash{$1} = $2;
			delete $opt_hash{$1} if $2 eq "off";
		} elsif ( /^-?(\S+)$/ ) {
			$opt_hash{$1} = "";
		} else {
			croak "Error: incorrect program option string $useroptstr";
		}
	}

	foreach $optkey (sort keys %opt_hash) {
		$optval = $opt_hash{$optkey};
		if ( $optval eq "" ) {
			push(@return_opt,"-$optkey");
		} else {
			push(@return_opt,"-$optkey $optval");
		}
	}

	return(join(" ",@return_opt));


}

sub sam_file_is_empty ($) {
  my $file = shift;
  my $fh = IO::File->new("<$file");

  while (my $line = $fh->getline) {
    next unless $line =~ /\S/;
    next if $line =~ /^@/;
    $fh->close;
    return(0);
  }
  $fh->close;
  return(1);

}




1;