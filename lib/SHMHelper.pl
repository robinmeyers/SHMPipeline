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

# sub smith_water_to_reference ($$$) {
# 	my $expt_id = shift;
# 	my $expt_hash = shift;
# 	my $swopt = shift;

# 	$expt_hash->{water} = $expt_hash->{exptdir} . "/${expt_id}.water";
# 	$expt_hash->{waternice} = $expt_hash->{exptdir} . "/${expt_id}_nice.water";

# 	(my $niceswopt = $swopt) =~ s/-aformat markx10//;
# 	System(join(" ","water","-asequence",$expt_hash->{reference},"-bsequence",$expt_hash->{fa},"-outfile",$expt_hash->{water},$swopt));
# 	System(join(" ","water","-asequence",$expt_hash->{reference},"-bsequence",$expt_hash->{fa},"-outfile",$expt_hash->{waternice},$niceswopt));

# }

# sub parse_smith_water_file ($$$$) {
# 	my $expt_id = shift;
# 	my $expt_hash = shift;
# 	my $min_score = shift;
# 	my $min_qual = shift;

# 	print "\n$expt_id: Parsing smith-waterman alignments...\n";

# 	my $swfh = IO::File->new("<".$expt_hash->{water}) or croak "Error: could not open water file for parsing";
# 	my $outfh = IO::File->new(">".$expt_hash->{mutfile}) or croak "Error: could not write to output";
# 	my $statfh = IO::File->new(">".$expt_hash->{clonefile}) or croak "Error: could not write to stats file";
# 	my $bxfh = IO::File->new(">".$expt_hash->{base_ex}) or croak "Error: could not write to base exchange file";

# 	my %q = %{$expt_hash->{q}};

# 	while ( my $line = $swfh->getline ) {
# 		chomp($line);
# 		last if $line =~ /^#-+$/;

# 		my $score;
# 		my $ident;
# 		my $overlap;
# 		my $refid;
# 		my $tstart;
# 		my $tend;
# 		my $tseq = "";
# 		my $clone;
# 		my $qstart;
# 		my $qend;
# 		my $qseq = "";

# 		if ($line =~ /^\>\>#\d+/) {


# 			until ($line =~ /^\>(\S+)\s\.\./) {
# 				last unless defined $line;
# 				$score = $1 if $line =~ /^; sw_score: (\d+)/ ;
# 				$ident = $1 if $line =~ /^; sw_ident: ([\d\.]+)/;
# 				$overlap = $1 if $line =~ /^; sw_overlap: (\d+)/;


# 				$line = $swfh->getline;
# 				chomp($line);
# 			}

# 			$refid = $1 if $line =~ /^\>(\S+)\s\.\./;

# 			$line = $swfh->getline;
# 			chomp($line);

# 			until ($line =~ /^\>(\S+)\s\.\./) {
# 				last unless defined $line;
# 				$tstart = $1 if $line =~ /^; al_start: (\d+)/;
# 				$tend = $1 if $line =~ /^; al_stop: (\d+)/;
# 				$tseq .= $1 if $line =~ /^([ACGTacgtNnXx-]+)/;

# 				$line = $swfh->getline;
# 				chomp($line);
# 			}

# 			$clone = $1 if $line =~ /^\>(\S+)\s\.\./;

# 			$line = $swfh->getline;
# 			chomp($line);

# 			until ($line =~ /^#/) {
# 				last unless defined $line;
# 				$qstart = $1 if $line =~ /^; al_start: (\d+)/;
# 				$qend = $1 if $line =~ /^; al_stop: (\d+)/;
# 				$qseq .= $1 if $line =~ /^([ACGTacgtNnXx-]+)/;

# 				$line = $swfh->getline;
# 				chomp($line);
# 			}

# 			#print join(" ",$clone,$score,$qstart,$qend,$tstart,$tend)."\n";
# 			if ($score >= $min_score) {
# 				unless ($clone =~ /\.[fr]$/) {

# 					$q{$clone}->{score} = $score;
# 					$q{$clone}->{ident} = $ident;
# 					$q{$clone}->{overlap} = $overlap;

# 					$q{$clone}->{refid} = $refid;
# 					$q{$clone}->{tstart} = $tstart;
# 					$q{$clone}->{tend} = $tend;
# 					$q{$clone}->{tseq} = $tseq;

# 					$q{$clone}->{qstart} = $qstart;
# 					$q{$clone}->{qend} = $qend;
# 					$q{$clone}->{qseq} = $qseq;

# 				} else {
# 					($clone, my $read) = $clone =~ /^(.*)\.([fr])$/;
# 					my $ori = $read eq 'f' ? 'for' : 'rev' ;
# 					$q{$clone}->{$ori}->{score} = $score;
# 					$q{$clone}->{$ori}->{ident} = $ident;
# 					$q{$clone}->{$ori}->{overlap} = $overlap;

# 					$q{$clone}->{$ori}->{refid} = $refid;
# 					$q{$clone}->{$ori}->{tstart} = $tstart;
# 					$q{$clone}->{$ori}->{tend} = $tend;
# 					$q{$clone}->{$ori}->{tseq} = $tseq;

# 					$q{$clone}->{$ori}->{qstart} = $qstart;
# 					$q{$clone}->{$ori}->{qend} = $qend;
# 					$q{$clone}->{$ori}->{qseq} = $qseq;
# 				}
# 			}
# 		}

# 	}
# 	$swfh->close;



# 	my $reffh = IO::File->new("<".$expt_hash->{reference}) or croak "Error: could not open reference file";
# 	my ($refid,$refseq) = read_fasta($reffh);
# 	$reffh->close;

# 	$statfh->print(join("\t",qw(ID Bp Subs Dels DelBp Ins InsBp RefA RefC RefG RefT RefN Coords Notes))."\n");

# 	foreach my $clone (sort keys %q) {
# 		my ($bps,$subs,$del,$delbp,$ins,$insbp) = parse_sw_alignment($outfh,$expt_id,$q{$clone},$min_qual);
# 		if (defined $q{$clone}->{bps}) {
# 			my @bases = qw(A C G T N);
# 			$statfh->print(join("\t",$clone,$q{$clone}->{bps},
# 																		$q{$clone}->{sub},
# 																		$q{$clone}->{del},
# 																		$q{$clone}->{delbp},
# 																		$q{$clone}->{ins},
# 																		$q{$clone}->{insbp},
# 																		$q{$clone}->{$bases[0]},
# 																		$q{$clone}->{$bases[1]},
# 																		$q{$clone}->{$bases[2]},
# 																		$q{$clone}->{$bases[3]},
# 																		$q{$clone}->{$bases[4]},
# 																		$q{$clone}->{coords})."\n");

# 			$bxfh->print(join("\t",$clone,@bases)."\n");
# 			foreach my $base (@bases) {
# 				next if $base eq "N";
# 				$bxfh->print(join("\t",$base,$q{$clone}->{bx}->{$base}->{$bases[0]},
# 																		 $q{$clone}->{bx}->{$base}->{$bases[1]},
# 																		 $q{$clone}->{bx}->{$base}->{$bases[2]},
# 																		 $q{$clone}->{bx}->{$base}->{$bases[3]})."\n");
# 			}
# 			$bxfh->print("\n");

# 		} else {
# 			$statfh->print(join("\t", $clone, 0,"","","","","","","","","","","","bad alignment")."\n");
# 		}
# 	}

# 	$statfh->close;
# 	$bxfh->close;

# }


# sub merge_alignments ($$) {

# 	my $q = shift;
# 	my $min_qual = shift;

# 	#print "\nmerging alignments for " . $q->{clone} . "\n";


# 	return if defined $q->{score};
# 	my $ori;
# 	if (!defined $q->{for}->{score}) {
# 		$ori = "rev";
# 	}
# 	if (!defined $q->{rev}->{score}) {
# 		$ori = "for";
# 	}
# 	if (defined $ori) {
# 		$q->{score} = $q->{$ori}->{score};
# 		$q->{tseq} = $q->{$ori}->{tseq};
# 		$q->{qseq} = $q->{$ori}->{qseq};
# 		$q->{contigqual} = $q->{$ori}->{quals};
# 		$q->{tstart} = $q->{$ori}->{tstart};
# 		$q->{tend} = $q->{$ori}->{tend};
# 		$q->{qstart} = $q->{$ori}->{qstart};
# 		$q->{qend} = $q->{$ori}->{qend};

# 		# print "only using $ori orientation\n";
# 		# print "tseq: ".$q->{tstart}."-".$q->{tend}."\n";
# 		# print "qseq: ".$q->{qstart}."-".$q->{qend}."\n";
# 	} else {
		
# 		# print "using both alignments\n";
# 		# print "forward: tseq ".$q->{for}->{tstart}."-".$q->{for}->{tend}." qseq ".$q->{for}->{qstart}."-".$q->{for}->{qend}."\n";
# 		# print "tseq: ".$q->{for}->{tseq}."\n";
# 		# print "qseq: ".$q->{for}->{qseq}."\n";
# 		# print "reverse: tseq ".$q->{rev}->{tstart}."-".$q->{rev}->{tend}." qseq ".$q->{rev}->{qstart}."-".$q->{rev}->{qend}."\n";
# 		# print "tseq: ".$q->{rev}->{tseq}."\n";
# 		# print "qseq: ".$q->{rev}->{qseq}."\n";

# 		$q->{score} = $q->{for}->{score} + $q->{rev}->{score};
# 		$q->{tstart} = min($q->{for}->{tstart},$q->{rev}->{tstart});
# 		$q->{tend} = max($q->{for}->{tend},$q->{rev}->{tend});
		
# 		my $i_f;
# 		my $i_r;
# 		my $tpos = $q->{tstart};
# 		my $qseq = "";
# 		my $tseq = "";
# 		my @quals = ();

# 		my $qpos_f;
# 		my $qpos_r;

# 		while ($q->{tend} >= $tpos) {

# 			if ($tpos == $q->{for}->{tstart}) {
# 				$qpos_f = $q->{for}->{qstart};
# 				$i_f = 0;
# 			}
# 			if ($tpos == $q->{rev}->{tstart}) {
# 				$qpos_r = $q->{rev}->{qstart};
# 				$i_r = 0;
# 			}


# 			if ($tpos > $q->{for}->{tend}) {
# 				$qpos_f = undef;
# 				$i_f = undef;
# 			}
# 			if ($tpos > $q->{rev}->{tend}) {
# 				$qpos_r = undef;
# 				$i_r = undef;
# 			}

# 			my $tbase_f;
# 			my $tbase_r;
# 			my $qbase_f;
# 			my $qbase_r;
# 			my $qual_f;
# 			my $qual_r;

# 			if (defined $qpos_f) {
# 				$tbase_f = substr($q->{for}->{tseq},$i_f,1);
# 				$qbase_f = substr($q->{for}->{qseq},$i_f,1);
# 				$qual_f = $q->{for}->{quals}->[$qpos_f-1];
# 			}
# 			if (defined $qpos_r) {
# 				$tbase_r = substr($q->{rev}->{tseq},$i_r,1);
# 				$qbase_r = substr($q->{rev}->{tseq},$i_r,1);
# 				$qual_r = $q->{rev}->{quals}->[$qpos_r-1];
# 			}

# 			if (defined $qbase_f && $qual_f >= $min_qual && defined $qbase_r && $qual_r >= $min_qual) {
# 				#compare both bases
# 				if ($tbase_f eq "-" xor $tbase_r eq "-") {
# 					if ($tbase_f eq "-") {
# 						$i_f++;
# 						$qpos_f++;
# 					}
# 					if ($tbase_r eq "-") {
# 						$i_r++;
# 						$qpos_r++;
# 					}
# 				} else {
# 					$tseq .= $tbase_f;
# 					$tpos++ unless $tbase_f eq "-";
# 					if ($qbase_f eq $qbase_r) {
# 						$qseq .= $qbase_f;
# 						push(@quals,$qual_f+$qual_r) unless $qbase_f eq "-";
# 					} else {
# 						$qseq .= "N";
# 						push(@quals,0);
# 					}
# 					$i_f++;
# 					$i_r++;
# 					$qpos_f++ unless $qbase_f eq "-";
# 					$qpos_r++ unless $qbase_r eq "-";
# 				}
# 			} elsif (defined $qbase_f && !(defined $qbase_r && $qual_r >= $min_qual)) {
# 				#use forward read
# 				$qseq .= $qbase_f;
# 				$tseq .= $tbase_f;
# 				push(@quals,$qual_f) unless $qbase_f eq "-";
# 				$tpos++ if $tbase_f ne "-";
# 				$qpos_f++ if $qbase_f ne "-";
# 				$qpos_r++ if (defined $qpos_r && $qbase_r ne "-");
# 				$i_f++;
# 				$i_r++ if defined $i_r;
# 			} elsif (defined $qbase_r && !(defined $qbase_f && $qual_f >= $min_qual)) {
# 				#use reverse read
# 				$qseq .= $qbase_r;
# 				$tseq .= $tbase_r;
# 				push(@quals,$qual_r) unless $qbase_r eq "-";
# 				$tpos++ if $tbase_r ne "-";
# 				$qpos_f++ if (defined $qpos_f && $qbase_f ne "-");
# 				$qpos_r++ if $qbase_r ne "-";
# 				$i_f++ if defined $qbase_f;
# 				$i_r++;
# 			} else {
# 				$tseq .= "N";
# 				$qseq .= "N";
# 				push(@quals,0);
# 				$tpos++;
# 				$qpos_f++ if defined $qpos_f;
# 				$qpos_r++ if defined $qpos_r;
# 				$i_f++ if defined $i_f;
# 				$i_r++ if defined $i_r;
# 			}



# 		}

# 		$q->{qstart} = 1;
# 		$q->{qend} = scalar @quals;

# 		$q->{tseq} = $tseq;
# 		$q->{qseq} = $qseq; 
# 		$q->{contigqual} = \@quals;

# 		# print "finished with tseq ".$q->{tstart}."-".$q->{tend}." qseq ".$q->{qstart}."-".$q->{qend}."\n";
# 		# print "tseq: ".$q->{tseq}."\n";
# 		# print "qseq: ".$q->{qseq}."\n";
# 		# print "quals: ".join(" ",@quals)."\n";
# 	}

		
# }


# sub parse_sw_alignment ($$$$) {

# 	my $fh = shift;
# 	my $expt= shift;
# 	my $q = shift;
# 	my $min_qual = shift;

# 	my $clone = $q->{clone};

# 	if (! defined $q->{score}) {
# 		if ( !defined $q->{for}->{score} && !defined $q->{rev}->{score} ) {
# 			return;
# 		} else {
# 			merge_alignments($q,$min_qual);
# 		}
# 	}

# 	if (defined $q->{score}) {

# 		my @tseq = split("",$q->{tseq});
# 		my @qseq = split("",$q->{qseq});
# 		my @quals = @{$q->{contigqual}};

# 		my $tstart = $q->{tstart};
# 		my $tend = $q->{tend};
# 		my $qstart = $q->{qstart};
# 		my $qend = $q->{qend};


# 		$q->{bps} = grep {$_ =~ /[ACGT]/} @tseq;
# 		$q->{sub} = 0;
# 		$q->{del} = 0;
# 		$q->{delbp} = 0;
# 		$q->{ins} = 0;
# 		$q->{insbp} = 0;

# 		my @pos_analyzed = ();
# 		my $pos = $tstart;
# 		my $start;
# 		my $end;
# 		foreach my $j (0 .. $#tseq) {
# 			next if $tseq[$j] eq "-";
# 			$start = $pos if (!defined $start && $tseq[$j] ne "N");
# 			$end = $pos if $tseq[$j] ne "N";
# 			if ($tseq[$j] eq "N") {
# 				push(@pos_analyzed,$start."-".$end) if (defined $start && defined $end);
# 				$start = undef;
# 				$end = undef;
# 			}
# 			$pos++;
# 		}
# 		push(@pos_analyzed,$start."-".$tend) if defined $start;
# 		$q->{coords} = join(",",@pos_analyzed);

# 		my @bases = qw(A C G T N);
# 		foreach my $base1 (@bases) {
# 			$q->{$base1} = 0;
# 			foreach my $base2 (@bases) {
# 				if ($base1 eq $base2) {
# 					$q->{bx}->{$base1}->{$base2} = "-";
# 				} else {
# 					$q->{bx}->{$base1}->{$base2} = 0;
# 				}
# 			}
# 		}

# 		my $i = 0;
# 		my $tpos = $tstart;
# 		my $qpos = $qstart;

# 		while ($i < scalar @qseq) {

# 			if ($qseq[$i] eq "-") {
# 				my $del_tstart = $tpos;
# 				my $del_qstart = $qpos-1;
# 				do { $i++; $tpos++; } while ($qseq[$i] eq "-");

# 				my @testquals = @quals[max($qstart-1,$del_qstart-3)..min($qend-1,$qpos+1)];
# 				#print(join(" ",$clone,"del",$del_tstart,$del_qstart-2,$qpos+2,@testquals)."\n");
# 				if (mean(@testquals) >= $min_qual) {
# 					$fh->print(join("\t",$expt,$clone,$del_tstart,"del","","",$tpos-$del_tstart,$tpos-1)."\n");
# 					$q->{del}++;
# 					$q->{delbp} += $tpos-$del_tstart;
# 				}
# 				next;
# 			}

# 			if ($tseq[$i] eq "-") {
# 				my $ins_qstart = $qpos;
# 				do { $i++; $qpos++; } while ($tseq[$i] eq "-");
# 				my @ins = @qseq[($i-$qpos+$ins_qstart)..($i-1)];
# 				my @testquals = @quals[($ins_qstart-1)..($qpos-2)];
# 				#print(join(" ",$clone,'ins',$tpos,$ins_qstart,$qpos-1,@testquals)."\n");

# 				my $size = $qpos-$ins_qstart;
# 				if (mean(@testquals) >= $min_qual) {
# 					$fh->print(join("\t",$expt,$clone,$tpos,"ins","","",$size,,"",join("",@ins))."\n");
# 					$q->{ins}++;
# 					$q->{insbp} += $size;
# 				}
# 				next;
# 			}


# 			if ($tseq[$i] ne $qseq[$i] && $tseq[$i] ne "N" && $qseq[$i] ne "N" && $qseq[$i] ne "X" ) {
# 				if ($quals[$qpos] >= $min_qual) {
# 						$fh->print(join("\t",$expt,$clone,$tpos,"sub",$tseq[$i],$qseq[$i])."\n");
# 						$q->{sub}++;
# 						$q->{bx}->{$tseq[$i]}->{$qseq[$i]}++;
# 					}
# 			}


# 			$q->{$tseq[$i]}++;
# 			$i++;
# 			$qpos++;
# 			$tpos++;
# 		}

# 	} 
# 	return;

# }

# sub phrap_sequence_pairs($$$) {
	
# 	my $expt_id = shift;
# 	my $expt_hash = shift;
# 	my $phrapopt = shift;

# 	print "\n$expt_id: Running phrap between read pairs using options $phrapopt\n";


# 	my %q;

# 	my $phredfh = IO::File->new("<".$expt_hash->{phred}) or croak "Error: could not open phred reads file";
# 	my $qualfh = IO::File->new("<".$expt_hash->{qual}) or croak "Error: could not open phred qual file";

# 	my $R1fh = Bio::SeqIO->new(-file => "<".$expt_hash->{R1}, -format => "fastq");
# 	my $R2fh = Bio::SeqIO->new(-file => "<".$expt_hash->{R2}, -format => "fastq");

# 	while (my $R1seq = $R1fh->next_seq()) {
# 		my $R2seq = $R2fh->next_seq();

# 		my $factory = Bio::Tools::Run::Phrap->new( -penalty => -2, -raw => 1 );

# 	}




# 	while (my ($head,$seq) = read_fasta($phredfh)) {
# 		my ($qual_head,$qualref) = read_qual($qualfh);

# 		my ($id,$len,$qualstart,$quallen,$type) = ($head =~ /^(\S+)\s+(\d+)\s+([\d\-]+)\s+(\d+)\s+(\S+)/) or croak "Error: unexpected sequence header $head";
# 		$qual_head =~ /^$id\s+/ or croak "Error: sequence and quality files don't match $head $qual_head";

# 		my ($path,$name,$ext) = parseFilename($id);

# 		my ($cloneid,$primer) = ($name =~ /^(\S+)_(\S+)$/);

# 		#my ($clone) = $cloneid =~ /([A-Za-z]-\d+)/;
# 		my $clone = $cloneid;

# 		my $ori;
# 		if ($primer eq $expt_hash->{forward}) {
# 			$ori = "for";
# 		} elsif ($primer eq $expt_hash->{reverse}) {
# 			$ori = "rev";
# 		} else {
# 			croak "Error: bad sequence id $id does not have either primer label";
# 		}
# 		$q{$clone}->{clone} = $clone;

# 		$q{$clone}->{$ori}->{seq} = $seq;
# 		$q{$clone}->{$ori}->{quals} = $qualref;

# 	}
# 	$phredfh->close;
# 	$qualfh->close;

# 	my $reffh = IO::File->new("<".$expt_hash->{reference});
# 	my ($head,$refseq) = read_fasta($reffh);
# 	my ($path,$refid,$ext) = parseFilename($expt_hash->{reference});
# 	$reffh->close;

# 	$expt_hash->{contig} = $expt_hash->{exptdir}."/${expt_id}_contigs.fa";
# 	$expt_hash->{contigqual} = $expt_hash->{exptdir}."/${expt_id}_contigs.fa.qual";

# 	my $contigfh = IO::File->new(">".$expt_hash->{contig}) or croak "Error: cannot write to contig file ";
# 	my $contigqualfh = IO::File->new(">".$expt_hash->{contigqual}) or croak "Error: cannot write to contig qual file";

# 	foreach my $clone (sort keys %q) {
# 		$q{$clone}->{fa} = $expt_hash->{phrapdir}."/$clone.fa";
# 		$q{$clone}->{qual} = $q{$clone}->{fa}.".qual";

# 		my $fafh = IO::File->new(">".$q{$clone}->{fa}) or croak "Error: cannot write to fasta file";
# 		my $qualfh = IO::File->new(">".$q{$clone}->{qual}) or croak "Error: cannot write to quals file";

# 		if (defined $q{$clone}->{for}->{seq} && length($q{$clone}->{for}->{seq}) > 0) {
# 			write_fasta($fafh,$clone.".f",$q{$clone}->{for}->{seq});
# 			write_qual($qualfh,$clone.".f",$q{$clone}->{for}->{quals});
# 		}
# 		if (defined $q{$clone}->{rev}->{seq} && length($q{$clone}->{rev}->{seq}) > 0) {
# 			write_fasta($fafh,$clone.".r",$q{$clone}->{rev}->{seq});
# 			write_qual($qualfh,$clone.".r",$q{$clone}->{rev}->{quals});
# 		}

# 		$fafh->close;
# 		$qualfh->close;

# 		if ((defined $q{$clone}->{for}->{seq} && length($q{$clone}->{for}->{seq}) > 0)
# 					|| (defined $q{$clone}->{rev}->{seq} && length($q{$clone}->{rev}->{seq}) > 0)) {

# 			my $phrap_cmd = join(" ","phrap",$q{$clone}->{fa},$phrapopt,">>",$expt_hash->{exptdir}."/phrap.out 2>&1");
# 			System($phrap_cmd,1);

# 			my $phrapfh = IO::File->new("<".$q{$clone}->{fa}.".contigs") or croak "Error: cannot open contig file";
# 			my $phrapqualfh = IO::File->new("<".$q{$clone}->{fa}.".contigs.qual") or croak "Error: cannot open qual file";

# 			my ($contigid,$contig) = read_fasta($phrapfh);
# 			my ($contigqualid,$contigqual) = read_qual($phrapqualfh);

# 			$phrapfh->close;
# 			$phrapqualfh->close;

# 			if (defined $contig) {
# 				write_fasta($contigfh,$clone,$contig);
# 				write_qual($contigqualfh,$clone,$contigqual);
# 				$q{$clone}->{contig} = $contig;
# 				$q{$clone}->{contigqual} = $contigqual;
# 			} else {

# 				if (defined $q{$clone}->{for}->{seq} && length($q{$clone}->{for}->{seq}) > 0) {
# 					write_fasta($contigfh,$clone.".f",$q{$clone}->{for}->{seq});
# 					write_qual($contigqualfh,$clone.".f",$q{$clone}->{for}->{quals});
# 				}
				
# 				if (defined $q{$clone}->{rev}->{seq} && length($q{$clone}->{rev}->{seq}) > 0) {
# 					$q{$clone}->{rev}->{seq} = reverseComplement($q{$clone}->{rev}->{seq});

# 					my @tmpquals = @{$q{$clone}->{rev}->{quals}};
# 					@tmpquals = reverse(@tmpquals);

# 					$q{$clone}->{rev}->{quals} = \@tmpquals;

# 					write_fasta($contigfh,$clone.".r",$q{$clone}->{rev}->{seq});
# 					write_qual($contigqualfh,$clone.".r",$q{$clone}->{rev}->{quals});
# 				}
# 			}
# 		}

# 	}

# 	$contigfh->close;
# 	$contigqualfh->close;

# 	$expt_hash->{q} = \%q;





# }




1;