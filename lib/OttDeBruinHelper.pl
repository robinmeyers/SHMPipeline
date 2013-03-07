#!/usr/bin/perl

use strict;
use warnings;
use Carp;
use File::Basename;
use IO::File;


sub argument {
	my $var = shift;
	my $description = shift;
	my $default = shift;

	return sprintf("  \%\-16s - %s\n",
		(defined $default ? sprintf("%s (%s)",$var,$default) : sprintf("%s",$var)),
		$description);
}

sub parseFilename ($) {
	my $fullname = shift;
	my ($name, $path, $ext) = fileparse($fullname, qr/\.\w{2,5}$/);
	return ($path, $name, $ext);
}

#invert exit status of system call
sub System ($) {
	my $cmd = shift;
	print "$cmd\n";
	my $status = system($cmd);
	return !$status;
}


sub manage_blast_options ($$) {

	my $defaultoptstr = shift;
	my $useroptstr = shift;

	return $defaultoptstr unless $useroptstr =~ /\S/;

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
			croak "Error: incorrect default blat option string $defaultoptstr";
		}
	}

	foreach (@useropt) {
		if ( /^-?(\S+)\s+(\S+)$/ ) {
			$opt_hash{$1} = $2;
			undef $opt_hash{$1} if $2 eq "off";
		} elsif ( /^-?(\S+)$/ ) {
			$opt_hash{$1} = "";
		} else {
			croak "Error: incorrect blat option string $useroptstr";
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


sub blast_to_reference_file ($$$) {
	my $expt_id = shift;
	my $expt_hash = shift;
	my $blastopt = shift;


	my $blastdir = $expt_hash->{exptdir} . "/blasts";
	unless (-d $blastdir) {
		mkdir $blastdir or croak "Error: could not create alignments directory for $expt_id";
	}

	$expt_hash->{blast} = "$blastdir/${expt_id}_blast.txt";
	
	System(join(" ","blastn","-subject",$expt_hash->{reference},"-query",$expt_hash->{raw},"-out",$expt_hash->{blast},$blastopt))
									or croak "Error: $expt_id failed during blast to reference sequence";
}


sub parse_blast_output ($$$$) {
	my $expt_id = shift;
	my $expt_hash = shift;
	my $stats_hash = shift;
	my $blast_header = shift;

	print "\n$expt_id: Parsing blast output...\n";

	my $infh = IO::File->new("<".$expt_hash->{blast});
	my $outfh = IO::File->new(">".$expt_hash->{output});
	$outfh->print(join("\t",qw(QueryID SubjectID DelStarts DelLengths InsStarts InsLengths QuerySeq SubjectSeq))."\n");

	my $csv = Text::CSV->new({sep_char => "\t"});
	$csv->column_names(@$blast_header);

	# For each blast alignment
	while( my $aln = $csv->getline_hr($infh) ) {
		# split the query and subject sequences into an arrays
		my @qseq = split("",$aln->{qseq});
		my @sseq = split("",$aln->{sseq});
		# move along unless the aligment is long enough
		next unless $aln->{length} > $expt_hash->{bindstart};

		my $sstart = $aln->{sstart} - 1;
		my $i = 0;
		my $gaps = 0;
		my ($bindstart,$cutstart,$cutend,$bindend);
		# Since gaps have been inserted into the reference sequence
		# (and maybe not all of it aligned) we have to re-calibrate the coordinates

		# bindstart
		while ( $i <= $expt_hash->{bindstart}-$sstart && $i+$gaps < scalar @sseq) {
			if ($sseq[$i+$gaps] eq "-") {
				$gaps++;
			} else {
				$i++;
			}
		}
		$bindstart = $i + $gaps - 1;
		# cutstart
		while ( $i <= $expt_hash->{cutstart}-$sstart && $i+$gaps < scalar @sseq ) {
			if ($sseq[$i+$gaps] eq "-") {
				$gaps++;
			} else {
				$i++;
			}
		}
		$cutstart = $i + $gaps - 1;
		#cutend
		while ( $i <= $expt_hash->{cutend}-$sstart && $i+$gaps < scalar @sseq ) {
			if ($sseq[$i+$gaps] eq "-") {
				$gaps++;
			} else {
				$i++;
			}
		}
		$cutend = $i + $gaps - 1;
		# Only count reads if the alignment extends through the cut site
		next unless scalar @sseq > $cutend;
		# bindend
		while ( $i <= $expt_hash->{bindend}-$sstart && $i+$gaps < scalar @sseq ) {
			if ($sseq[$i+$gaps] eq "-") {
				$gaps++;
			} else {
				$i++;
			}
		}
		$bindend = $i + $gaps - 1;

		# Initialize arrays to hold the starting location and length
		# of each indel
		my @del_starts = ();
		my @del_lens = ();
		my @ins_starts = ();
		my @ins_lens = ();

		# Start at the cutstart
		$i = $cutstart;
		while ($i <= $cutend) {
			# move along unless there's a gap in either the subject or query
			unless ( $qseq[$i] eq "-" || $sseq[$i] eq "-" ) {
				$i++;
				next;
			}
			# gap on the query sequence is a deletion
			if ( $qseq[$i] eq "-" ) {
				# if the first base in the cut region is a gap we have to backtrack
				# to find the start of that indel
				if ($i == $cutstart) {
					while ($qseq[--$i] eq "-") {}
					$i++;
				}
				# push the start point onto the array
				push(@del_starts,$i+1);
				my $j = 0;
				# move forward until the gap ends
				while ($qseq[$i+(++$j)] eq "-") {}
				# push the length onto the array
				push(@del_lens,$j);
				$i += $j;
				next;
			}
			# gap on the subject sequence is an insertion
			if ( $sseq[$i] eq "-" ) {
				if ($i == $cutstart) {
					while ($sseq[--$i] eq "-") {}
					$i++;
				}
				# push the start point onto the array
				push(@ins_starts,$i+1);
				my $j = 0;
				# move forward until the gap ends
				while ($sseq[$i+(++$j)] eq "-") {}
				# push the length onto the array
				push(@ins_lens,$j);
				$i += $j;
				next;
			}
		}
		# print the results out in the output file
		$outfh->print(join("\t",	$aln->{qseqid},
															$aln->{sseqid},
															join(",",@del_starts),
															join(",",@del_lens),
															join(",",@ins_starts),
															join(",",@ins_lens),
															substr($aln->{qseq},$bindstart,$bindend-$bindstart),
															substr($aln->{sseq},$bindstart,$bindend-$bindstart))."\n");

		# update our summary stats hash

		if (scalar @del_lens == 0 && scalar @ins_lens == 0) {
			# no insertions or deletions
		} elsif (scalar @del_lens == 1 && $del_lens[0] == 1 && scalar @ins_lens == 0) {
			# single bp deletion, no insertion
			$stats_hash->{singledel}++;
		} elsif (scalar @ins_lens == 1 && $ins_lens[0] == 1 && scalar @del_lens == 0) {
			# single bp insertion, no deletion
			$stats_hash->{singleins}++;
		} elsif (scalar @del_lens == 1 && scalar @ins_lens == 0) {
			# multi bp deletion, no insertion
			$stats_hash->{dels}++;
		} elsif (scalar @del_lens == 0 && scalar @ins_lens == 1) {
			# multi bp insertion, no deletion
			$stats_hash->{ins}++;
		} else {
			# otherwise multiple indels - complex
			$stats_hash->{complex}++;
		}

		$stats_hash->{total}++;
	}
	$infh->close;
	$outfh->close;
}

1;