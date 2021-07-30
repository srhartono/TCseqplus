package Mutationlib;

use strict; use warnings; use mitochy;
use vars qw(@EXPORT);
use parent "Exporter";

our @EXPORT = qw(
colorize
blank_if_undef
parse_fasta
parse_aln
get_mh
get_igtype
get_sampletype
fix_pos
fix_pos2
);

sub fix_pos {
   my ($poses, $seqprev, $seqcurr, $outBigLog) = @_;
	my @poses = split(",", $poses);
	my @respos;
	foreach my $pos (sort {$a <=> $b} @poses) {
		if ($pos < 0) {
			push(@respos, $pos);
			next;
		}
	   my @seqprev = split("", $seqprev);
	   my @seqcurr = split("", $seqcurr);
	   my $ind = -1;
	   my $seqp = "";
	   for (my $i = 0; $i < @seqprev; $i++) {
	      $ind ++ if $seqprev[$i] ne "-";
	      $seqp .= $seqprev[$i] if $seqprev[$i] ne "-";
	   #   print "$i. ind=$ind poswant=$pos seq=$seqp (seqprev[i] = $seqprev[$i])\n" if $i eq $pos or $i == @seqprev - 1;
	      last if $i eq $pos;
	   }
	   #print "new indwant=$ind\n";
	   $pos = $ind;
	   my $pos2 = -1;
	   $ind = -1;
	   my $seqc = "";
	   for (my $i = 0; $i < @seqcurr; $i++) {
	      $ind ++ if $seqcurr[$i] ne "-";
	      $seqc .= $seqcurr[$i] if $seqcurr[$i] ne "-";
	      $pos2 = $i if $ind eq $pos and $seqcurr[$i] ne "-";
	      last if $pos2 ne -1;
	   }
		push(@respos, $pos2);
	}

	LOG($outBigLog, "\nfix_pos:\n$poses\n" . join(",", @respos) . "\n-------------------\n\n") if defined $outBigLog;
	return(join(",", @respos));
}

sub fix_pos2 {
   my ($poses, $seqprev, $seqcurr, $outBigLog) = @_;
	my @origpos;
	my @poses = split(",", $poses);
	my @respos;
	my %pos;
	my %ind;
	foreach my $pos (sort {$a <=> $b} @poses) {
		$pos{$pos}{check} = 1;
		if ($pos < 0) {
			$pos{$pos}{ind0} = $pos;
			$pos{$pos}{ind1} = $pos;
		}
	}
#	foreach my $pos (sort {$a <=> $b} @poses) {
#		if ($pos < 0) {
#			push(@respos, $pos);
#			next;
#		}
	   my @seqprev = split("", $seqprev);
	   my @seqcurr = split("", $seqcurr);
	   my $ind = -1;
	   my $seqp = "";
		#print "poswant=$pos\n";
	   for (my $i = 0; $i < @seqprev; $i++) {
	      $ind ++ if $seqprev[$i] ne "-";
	      $seqp .= $seqprev[$i] if $seqprev[$i] ne "-";
	   #   print "$i. ind=$ind poswant=$pos seq=$seqp (seqprev[i] = $seqprev[$i])\n" if $i eq $pos or $i == @seqprev - 1;
			
			if (defined $pos{$i}) {
				$ind{$ind} = $i if not defined $pos{$i}{ind0};
				$pos{$i}{ind0} = $ind if not defined $pos{$i}{ind0};
			}


#	      last if $i eq $pos;

	   }
	   #print "new indwant=$ind\n";


#	   $pos = $ind;
#	   my $pos2 = -1;
	   $ind = -1;
	   my $seqc = "";
	   for (my $i = 0; $i < @seqcurr; $i++) {
	      $ind ++ if $seqcurr[$i] ne "-";
	      $seqc .= $seqcurr[$i] if $seqcurr[$i] ne "-";
##	   #   print "$i. ind=$ind indwant=$pos seq=$seqc (seqcurr[i] = $seqcurr[$i])\n" if $pos2 ne -1;
			if (defined $ind{$ind} and $seqcurr[$i] ne "-") {
				my $previ = $ind{$ind};
				$pos{$previ}{ind1} = $i;
			}
#	      $pos2 = $i if $ind eq $pos and $seqcurr[$i] ne "-";
#	      last if $pos2 ne -1;
	   }
	foreach my $previ (sort {$a <=> $b} keys %pos) {
		my $ind0 = $pos{$previ}{ind0};
		my $ind1 = $pos{$previ}{ind1};
		push(@origpos, $previ);
		push(@respos, $ind1);
	}
	LOG($outBigLog, "\n\n" . join(",", @origpos) . "\n" . join(",", @respos) . "\n\n------------------------------->\n") if defined $outBigLog;
	return(join(",", @respos));
}


sub get_sampletype {
	my ($sampleID) = @_; #W1
	my ($sample0, $treat0, $cut) = $sampleID =~ /^([A-Z])([0-9]+)(cut)?$/;
	$cut = " (B. CUT)" if defined $cut;
	$cut = " (A. NOT_CUT)" if not defined $cut;

	my $sample = $sample0 eq "W" ? "1. WT" :
					 $sample0 eq "S" ? "2. Setx KO" :
					 $sample0 eq "R" ? "3. RNH2b KO" :
					 $sample0 eq "D" ? "4. Setx RNH2b DKO" : "5. UNK SAMPLEID $sample0";

	my $treat = $treat0 eq "1" ? "1. IgG1 undig germ$cut" :
 					 $treat0 eq "2" ? "2. IgG3 undig rich$cut" :
					 $treat0 eq "3" ? "3. IgG1 diges rich$cut" : "4. UNK TREAT $treat0$cut";

	return($sample, $treat);
}
sub get_igtype {
	my ($junc1pos, $junc2pos, $chr, $beg, $end, $strand, $outLog) = @_;
#my @beg0ind =   (114650000,114590000 ,114564000 ,114540000  ,114520000  ,114508000,114490000,114000000         ,114675001);
#my @end0ind =   (114675000,114649999 ,114589999 ,114563999  ,114539999  ,114519999,114507999,114489999         ,114999999);
my @beg0ind =   (114656769,114594435 ,114563451 ,114541176  ,114523536  ,114507471,114493041,114000000         ,114675001);
my @end0ind =   (114670052,114609628 ,114582590 ,114559110  ,114539403  ,114520744,114508309,114489999         ,114999999);
my @junctypes = ("01. IgM","02. IgG3","03. IgG1","04. IgG2b","05. IgG2c","06. IgE","07. IgA","08. IgA Down100k","09. IgM Up250k","10. Chr12","11. Other Chrs");

	my ($junctype, $junc1dist, $junc2dist, $length1, $length2);
	$junc1dist = $junc1pos - $beg0ind[0];
	$length1 = $end0ind[0] - $beg0ind[0];

#      my ($igtype, $junc1dist, $junc2dist) = get_igtype($junc1, $junc2, $chr2, $beg2, $end2, $strand2, $outLog);

	if ($chr ne "chr12") {
		$junctype = $junctypes[10];
		$junc2dist = $end - $beg;
		$length2 = $end - $beg < 10000 ? 10000 : $end - $beg;
#		print "1. junc1pos=$junc1pos, junc2pos=$junc2pos, chr=$chr:$beg-$end ($strand), junctype=$junctype, junc1dist=$junc1dist, junc2dist=$junc2dist, length1=$length1, length2=$length2\n";
		return($junctype, $junc1dist, $junc2dist, $length1, $length2);
	}
	else {
		for (my $i = 0; $i < @beg0ind; $i++) {
			my $beg0 = $beg0ind[$i];
			my $end0 = $end0ind[$i];
	
			$junctype = $junctypes[$i];
			$junc2dist = $strand eq "-" ? $junc2pos - $beg0 : $end0 - $junc2pos;
			$length2 = $end0 - $beg0;
	
			my $junc = $strand eq "-" ? $beg : $end;
	
			my ($int) = intersect($beg0, $end0, $junc, $junc);
			#die "$int ($beg0-$end0, $junc)\n" if $int == 1;
			if ($int == 1) {
#				print "2. junc1pos=$junc1pos, junc2pos=$junc2pos, beg0-end0=$beg0-$end0, junc=$junc, chr=$chr:$beg-$end ($strand), junctype=$junctype, junc1dist=$junc1dist, junc2dist=$junc2dist, length1=$length1, length2=$length2\n";
				return($junctype, $junc1dist, $junc2dist, $length1, $length2);
			}
	
		}
		$junctype = $junctypes[9];
		$junc2dist = -1; #$end - $beg;
		$length2 = -1; #$end - $beg < 10000 ? 10000 : $end - $beg;
#		print "3. junc1pos=$junc1pos, junc2pos=$junc2pos, chr=$chr:$beg-$end ($strand), junctype=$junctype, junc1dist=$junc1dist, junc2dist=$junc2dist, length1=$length1, length2=$length2\n";
		return($junctype, $junc1dist, $junc2dist, $length1, $length2);
	}


#	if ($junc >= 114490000 and $junc <= 114675000) {	
#		DIELOG($outLog, "$chr:$beg-$end (using $LGN$junc$N) wasn't counted as IgN (gonna be \"$junctypes[9]\")\n");
#	}
}

sub colorize {
	my ($seqs, $juncpos, $reverse) = @_;
	die "Not defined seqs!\n" if not defined $seqs;
	my $reverseprint = defined $reverse ? "REVERSEPRINTCOLOR" : "colorrev undef";
	($seqs) = mitochy::revcomp($seqs) if defined $reverse;
	my @seq = split("", $seqs);
	my $total_del = $seqs =~ tr/\-/\-/;
	my $ind = defined $reverse ? @seq - $total_del + 1 : 0;
	my $print = "";
	for (my $i = 0; $i < @seq; $i++) {
		die "undef i=$i seq=$seqs\n" if not defined $seq[$i];
#		if ($seq[$i] ne "-") {
#		}
		if (defined $juncpos and $juncpos ne "NA" and $ind >= $juncpos - 2 and $ind <= $juncpos + 3) {
#		if (defined $juncpos and $juncpos ne "NA" and $i >= $juncpos - 2 and $i <= $juncpos + 3) {
			if (($ind == $juncpos or $ind == $juncpos + 1) and $seq[$i] ne "-") {
#			if (($i == $juncpos or $i == $juncpos + 1)) {
				$print .= "$LCY$seq[$ind]$N" if $ind == $juncpos;
				$print .= "$LGN$seq[$ind]$N" if $ind == $juncpos + 1;
			}
			else {
				$print .= "$N$seq[$ind]$N";
			}
#			$print .= "$N$seq[$i]$N" if $ind >= $juncpos - 2 and $ind < $juncpos;
#			$print .= "$N$seq[$i]$N" if $ind <= $juncpos + 3 and $ind > $juncpos + 1;
		}
		else {
			die "undef i=$ind seq=$seqs\n" if not defined $seq[$ind];
			$print .= $LGN . $seq[$ind] if $seq[$ind] =~ /^A$/;
			$print .= $YW . $seq[$ind] if $seq[$ind] =~ /^C$/;
			$print .= $LPR . $seq[$ind] if $seq[$ind] =~ /^G$/;
			$print .= $LRD . $seq[$ind] if $seq[$ind] =~ /^T$/;
			$print .= $N . $seq[$ind] if $seq[$ind] =~ /^a$/;
			$print .= $N . $seq[$ind] if $seq[$ind] =~ /^c$/;
			$print .= $N . $seq[$ind] if $seq[$ind] =~ /^g$/;
			$print .= $N . $seq[$ind] if $seq[$ind] =~ /^t$/;
			$print .= "$N$seq[$ind]" if $seq[$ind] !~ /^[acgtACGT]$/;
		}
		$ind = defined $reverse ? $ind - 1 : $ind + 1;
	}
	$print .= "$N";
	return($print);
}

sub blank_if_undef {
	my @strings = @_;
	for (my $i = 0; $i < @strings; $i++) {
		$strings[$i] = "" if not defined $strings[$i];
	}
	return(@strings);
}

sub parse_fasta {
	my ($reses, $resHASH, $restype, $Xbuffer, $outLog, $namewant) = @_;
	my @res = @{$reses};
	my $def;
	my %res;
	for (my $i = 0; $i < @res; $i ++) {
		my $line = $res[$i];
		chomp($line);
		if ($line =~ /^>/) {
			$def = $line;
			$def =~ s/^>//;
			$res{$def} = "";
		}
		else {
			$res{$def} .= $line;
		}
	}
	foreach my $def (sort keys %res) {
		$res{$def} =~ s/$Xbuffer//g;# if $res{$def} =~ /X$/;
		$res{$def} =~ s/$Xbuffer//g;# if $res{$def} =~ /^X/;
		$res{$def} =~ s/[UVWXYZ]//g;
		print $outLog date() . " Mutationlib.pm::parse_fasta: def $LCY$def$N has extra X:\n$res{$def}\n\n" if $res{$def} =~ /VWXYZ/;
		die "\n" if $res{$def} =~ /VWXYZ/;
		$resHASH->{$def}{$restype} = $res{$def};
	}
	return($resHASH);
}


sub get_mh {
	my ($outLog, $name, $seqQ1L, $seqQ1R, $seqC1L, $seqC1R, $seqC2L, $seqC2R, $junc, $muscleparam, $namewant) = @_;
	my $printmhdebug = "\n----------------------BELOW IS MH-------------------\n";
	$outLog = ${$outLog};
	#seqC1LL	seqC1LR	seqC1RL	seqC1RR
	#seqQ1LL	seqQ1LR	seqQ1RL	seqQ1RR
	#seqC2LL	seqC2LR	seqC2RL	seqC2RR

	my ($seqC1LL, $seqC1LR, $seqC1RL, $seqC1RR, 
		 $seqQ1LL, $seqQ1LR, $seqQ1RL, $seqQ1RR, 
		 $seqC2LL, $seqC2LR, $seqC2RL, $seqC2RR);

	($seqC1LL, $seqC1LR) = $seqC1L =~ /^(.*)(.{20})$/;
	($seqQ1LL, $seqQ1LR) = $seqQ1L =~ /^(.*)(.{20})$/;
	($seqC2LL, $seqC2LR) = $seqC2L =~ /^(.*)(.{20})$/;

	($seqC1RL, $seqC1RR) = $seqC1R =~ /^(.{20})(.*)$/;
	($seqQ1RL, $seqQ1RR) = $seqQ1R =~ /^(.{20})(.*)$/;
	($seqC2RL, $seqC2RR) = $seqC2R =~ /^(.{20})(.*)$/;

	($seqC1LL, $seqC1LR, $seqC1RL, $seqC1RR, 
		 $seqQ1LL, $seqQ1LR, $seqQ1RL, $seqQ1RR, 
		 $seqC2LL, $seqC2LR, $seqC2RL, $seqC2RR) = 
	blank_if_undef($seqC1LL, $seqC1LR, $seqC1RL, $seqC1RR, 
		 $seqQ1LL, $seqQ1LR, $seqQ1RL, $seqQ1RR, 
		 $seqC2LL, $seqC2LR, $seqC2RL, $seqC2RR);


	my @seqQ1LR = split("", $seqQ1LR);
	my @seqQ1RL = split("", $seqQ1RL);
	my @seqC1LR = split("", $seqC1LR);
	my @seqC1RL = split("", $seqC1RL);
	my @seqC2LR = split("", $seqC2LR);
	my @seqC2RL = split("", $seqC2RL);
	
	my ($mhtot, $mhbeg, $mhend, $mhL, $mhR, $lastmat, $printmat) = (0, -1, -1, -1, -1, 0, "");


	#seqC1LL	[seqC1LR]	seqC1RL	seqC1RR
	#seqQ1LL	[seqQ1LR]	seqQ1RL	seqQ1RR
	#seqC2LL	[seqC2LR]	seqC2RL	seqC2RR


	# find beg of microhomology chunk.
	# starts from junction of seqQ1 (seqQ1LR) to the left
	my $ind = 0;
	for (my $i = @seqQ1LR - 1; $i >= 0; $i --) {
		# if it matches exactly at junction break then flag as match (lastmat = 1)
		$lastmat = 1 if $i == @seqQ1LR - 1 and $seqQ1LR[$i] eq $seqC1LR[$i] and $seqQ1LR[$i] eq $seqC2LR[$i];

		# if match and prev is also match (lastmat == 1) then add to microhomology
		if ($lastmat == 1 and $seqQ1LR[$i] eq $seqC1LR[$i] and $seqQ1LR[$i] eq $seqC2LR[$i]) {
			$mhtot ++;
			$printmat = $seqQ1LR[$i] . $printmat;
#			$mhbeg = $i;
#			$mhend = $i if $mhend ne -1;
			$mhL = $ind; #how many bp left of junction mhbeg is
		}
		else { #otherwise mh stops there
			$printmat = "-" . $printmat;
			$lastmat = 0;
		}
		$ind ++;
	}
	$printmat .= ";"; #debug

	#seqC1LL	seqC1LR	[seqC1RL]	seqC1RR
	#seqQ1LL	seqQ1LR J	[seqQ1RL]	seqQ1RR
	#seqC2LL	seqC2LR	[seqC2RL]	seqC2RR

	# find end of microhomology chunk
	# starts from junction of seqQ1 (seqQ1RL at 0) to the right
	$lastmat = 0;
	for (my $i = 0; $i < @seqQ1RL; $i ++) {
		# if it matches right at junction break, then flag as match
		$lastmat = 1 if $i == 0 and $seqQ1RL[$i] eq $seqC1RL[$i] and $seqQ1RL[$i] eq $seqC2RL[$i];

		# if match and prev is also match (lastmat == 1) then add to microhomology
		if ($seqQ1RL[$i] eq $seqC1RL[$i] and $seqQ1RL[$i] eq $seqC2RL[$i] and $lastmat eq 1) {
			$mhtot ++;
#			$mhbeg = $i + $junc if $mhbeg eq -1;
#			$mhend = $i + $junc;
			$mhR = $i; #how many bp right of junction mhend is
#			$printmhdebug .= "i = $i, mhR = $mhR, seqQ1RL = $seqQ1RL[$i], same as seqC1RL $seqC1RL[$i] and seqC2RL $seqC2RL[$i], lastmat is $lastmat\n" if $name eq $namewant;
			$printmat = $printmat . "$seqQ1RL[$i]";
		}
		else { #otherwise stop here
			$printmat = $printmat . "-";
			$lastmat = 0;
		}
	}
#	$printmhdebug .= "ORIG: mhL = $mhL, mhR = $mhR\n";
	$mhbeg = $mhL eq -1 ? -1 : $junc - ($mhL);
	$mhend = $mhR eq -1 ? -1 : $junc + ($mhR + 1);
	#for printing purposes
	my $longest = length($seqC1LL);
	$longest = length($seqQ1LL) if $longest < length($seqQ1LL);
	$longest = length($seqC2LL) if $longest < length($seqC2LL);
	my ($seqC1LLspace, $seqQ1LLspace, $seqC2LLspace) = ("", "", "");
	if (length($seqC1LL) < $longest) {
		$seqC1LLspace = join("", (" ") x ($longest - length($seqC1LL)));
	}
	if (length($seqQ1LL) < $longest) {
		$seqQ1LLspace = join("", (" ") x ($longest - length($seqQ1LL)));
	}
	if (length($seqC2LL) < $longest) {
		$seqC2LLspace = join("", (" ") x ($longest - length($seqC2LL)));
	}


	$printmhdebug .= "\njunc=$junc\n";
	$printmhdebug .= join("", (" ") x length("1_con_$name")) . "\t" . join("", (" ") x $longest) . ";" . $printmat . ";" . "\n";
	$printmhdebug .= "1_con_$name\t" . $seqC1LLspace . $seqC1LL . "" . colorize($seqC1LR) . ";" . colorize($seqC1RL) . "" . $seqC1RR . "\n";
	$printmhdebug .= "1_seq_$name\t" . $seqQ1LLspace . $seqQ1LL . "" . colorize($seqQ1LR) . ";" . colorize($seqQ1RL) . "" . $seqQ1RR . "\n";
	$printmhdebug .= "2_con_$name\t" . $seqC2LLspace . $seqC2LL . "" . colorize($seqC2LR) . ";" . colorize($seqC2RL) . "" . $seqC2RR . "\n";
	######### printing purposes end

	$printmhdebug .= "mhL = $mhL, mhR = $mhR\n";

	my %mh;
	$mh{seqC1LL} = $mhL eq -1 ? $seqC1LL . $seqC1LR : $seqC1LL . join("", @seqC1LR[0..(@seqC1LR-1-($mhL+1))]);
	$mh{seqC1LR} = $mhL eq -1 ? ""                  : join("", @seqC1LR[(@seqC1LR-1-$mhL)..(@seqC1LR-1)]);
	$mh{seqC1RL} = $mhR eq -1 ? ""                  : join("", @seqC1RL[0..$mhR]);
	$mh{seqC1RR} = $mhR eq -1 ? $seqC1RL . $seqC1RR : join("", @seqC1RL[($mhR+1)..(@seqC1RL-1)]) . $seqC1RR;
	$mh{seqQ1LL} = $mhL eq -1 ? $seqQ1LL . $seqQ1LR : $seqQ1LL . join("", @seqQ1LR[0..(@seqQ1LR-1-($mhL+1))]);
	$mh{seqQ1LR} = $mhL eq -1 ? ""                  : join("", @seqQ1LR[(@seqQ1LR-1-$mhL)..(@seqQ1LR-1)]);
	$mh{seqQ1RL} = $mhR eq -1 ? ""                  : join("", @seqQ1RL[0..$mhR]);
	$mh{seqQ1RR} = $mhR eq -1 ? $seqQ1RL . $seqQ1RR : join("", @seqQ1RL[($mhR+1)..(@seqQ1RL-1)]) . $seqQ1RR;
	$mh{seqC2LL} = $mhL eq -1 ? $seqC2LL . $seqC2LR : $seqC2LL . join("", @seqC2LR[0..(@seqC2LR-1-($mhL+1))]);
	$mh{seqC2LR} = $mhL eq -1 ? ""                  : join("", @seqC2LR[(@seqC2LR-1-$mhL)..(@seqC2LR-1)]);
	$mh{seqC2RL} = $mhR eq -1 ? ""                  : join("", @seqC2RL[0..$mhR]);
	$mh{seqC2RR} = $mhR eq -1 ? $seqC2RL . $seqC2RR : join("", @seqC2RL[($mhR+1)..(@seqC2RL-1)]) . $seqC2RR;

	$printmhdebug .= "1_con_$name\t$seqC1LLspace$mh{seqC1LL}" . colorize($mh{seqC1LR}) . ";" . colorize($mh{seqC1RL}) . "$mh{seqC1RR}$N\n";
	$printmhdebug .= "1_seq_$name\t$seqQ1LLspace$mh{seqQ1LL}" . colorize($mh{seqQ1LR}) . ";" . colorize($mh{seqQ1RL}) . "$mh{seqQ1RR}$N\n";
	$printmhdebug .= "2_con_$name\t$seqC2LLspace$mh{seqC2LL}" . colorize($mh{seqC2LR}) . ";" . colorize($mh{seqC2RL}) . "$mh{seqC2RR}$N\n";
	my $Xbuffer = join("", ("XVWXY") x 3);
	my $cmdL = ">1_con\n$mh{seqC1LL}$Xbuffer\n>1_seq\n$mh{seqQ1LL}$Xbuffer\n>2_con\n$mh{seqC2LL}$Xbuffer\n";
	my $cmdR = ">1_con\n$Xbuffer$mh{seqC1RR}\n>1_seq\n$Xbuffer$mh{seqQ1RR}\n>2_con\n$Xbuffer$mh{seqC2RR}\n";

	my @resL = `echo '$cmdL' | muscle $muscleparam 2>/dev/null`;
	my @resR = `echo '$cmdR' | muscle $muscleparam 2>/dev/null`;
	my $resHASH;
	$resHASH = parse_fasta(\@resL, $resHASH, "L", $Xbuffer, $outLog, $namewant);
	$resHASH = parse_fasta(\@resR, $resHASH, "R", $Xbuffer, $outLog, $namewant);
	$mh{'1_con'} = $resHASH->{'1_con'}{L} . $mh{seqC1LR} . $mh{seqC1RL} . $resHASH->{'1_con'}{R};
	$mh{'1_seq'} = $resHASH->{'1_seq'}{L} . $mh{seqC1LR} . $mh{seqC1RL} . $resHASH->{'1_seq'}{R};
	$mh{'2_con'} = $resHASH->{'2_con'}{L} . $mh{seqC1LR} . $mh{seqC1RL} . $resHASH->{'2_con'}{R};
	$mh{'3_con'} = $resHASH->{'1_con'}{L} . $mh{seqC1LR} . $mh{seqC1RL} . $resHASH->{'2_con'}{R};

#	print "\ncmdL:\n" . $cmdL . "\n\n" if $name eq "M02034:489:000000000-CYYL8:1:1101:7360:4607";
#	print "\ncmdR:\n" . $cmdR . "\n\n" if $name eq "M02034:489:000000000-CYYL8:1:1101:7360:4607";
#	print "\nresL:\n" . join("", @resL) . "\n\n" if $name eq "M02034:489:000000000-CYYL8:1:1101:7360:4607";
#	print "\nresR:\n" . join("", @resR) . "\n\n" if $name eq "M02034:489:000000000-CYYL8:1:1101:7360:4607";

	$printmhdebug .= "\n\nmh=$mhL-$mhR\n\n";
	$printmhdebug .= "1_con_$name\t" . colorize($mh{'1_con'}) . "\n";
	$printmhdebug .= "1_seq_$name\t" . colorize($mh{'1_seq'}) . "\n";
	$printmhdebug .= "2_con_$name\t" . colorize($mh{'2_con'}) . "\n";


	
	$printmhdebug .= "\n1. JUNC=$junc, MHL=$mhL (mhbeg=junc-mhL=$mhbeg), MHR=$mhR (end=junc+mhR+1=$mhend)\n";
	print {$outLog} date() . "\n" . "\n$printmhdebug\n\n";
	#LOG($outLog, date() , "\n" . "MH MH MH\n$printmhdebug\n");# if $name eq $namewant;
	return($mhtot, $mhbeg, $mhend, \%mh);
	#seqC1LL			seqC1LR			seqC1RL	seqC1RR
	#seqC2LL	seqC2LR	seqC2RL			seqC2RR
}
#print "MH = $mhtot ($mhbeg-$mhend)\n";


###########################
# PARSE ALN #
###########################

sub parse_aln {
	my ($lines, $leader, $type, $muthash, $readorient, $lenMaxL, $lenMaxR, $name, $gap, $outBigLog, $outLog, $namewant, $reverse) = @_;
	my $reverseprint = defined $reverse ? $reverse : "rev undef";
	#print "reverse is $reverseprint, type is $type, leader is $leader, lines\n";
	my $die = 0;
	my ($beg0, $end0) = (-1,-1);
	my ($primer_pos_fix, $adapter_pos_fix, $mhbeg, $mhend) = (-1,-1,-1,-1);
	my $junc = 0;
	my $juncpos = 0;
	my ($mhbegorig, $mhendorig) = ($mhbeg, $mhend);
	if ($type =~ /^(primer|adapter)/) {
		($type, $beg0, $end0, $junc) = split(",", $type);
		$juncpos = $junc if $junc ne 0;
	}
	elsif ($type =~ /^none/) {
		($type, $beg0, $end0, $junc, $primer_pos_fix, $adapter_pos_fix, $mhbeg, $mhend) = split(",", $type);
#		print "mh=$mhbeg, $mhend\n";
		$juncpos = $junc if $junc ne 0;
		$mhbeg = -1 if not defined $mhbeg;
		$mhend = -1 if not defined $mhend;
		$mhbegorig = $mhbeg;
		$mhendorig = $mhend;
		$mhbeg = -1 if $mhbeg eq 0;
		$mhend = -1 if $mhend eq 0;

		if ($mhbeg ne -1 and $mhend eq -1) {
	#		$mhbeg = $juncpos - $mhbeg;
			$mhend = $juncpos;
		}
		elsif ($mhbeg eq -1 and $mhend ne -1) {
			$mhbeg = $juncpos + 1;
	#		$mhend = $juncpos + 1 + $mhend;
		}
		elsif ($mhbeg ne -1 and $mhend ne -1) {
	#		($mhbeg) = ($juncpos - $mhbeg);
	#		($mhend) = ($juncpos + 1 + $mhend);
		}
	}
	else {
		($beg0, $end0, $junc) = (-1, -1, 0);
	}


	my @line = @{$lines};
	my $begposorig = $beg0;
	my $endposorig = $end0;
	#fix end0 in case seq2 is shorter than seq1
	my %temp;
	my $tempdef = "INIT";
	for (my $i = 0; $i < @line; $i++) {
#		print "$i. $line[$i]\n";
		my $line = $line[$i]; chomp($line);
		if ($line =~ /^>/) {
			$tempdef = $line; $tempdef =~ s/^>//;
			$temp{$tempdef} = "";
		}
		else {
			$temp{$tempdef} .= $line;
		}
	}
	
	my $short = -1;
	my $longbeg = -1;
	my %short;
	my %long;
	my $info = "\n>${LGN}INFO$N\n";
	if ($type =~ /none/) {
		foreach my $tempdef1 (sort keys %temp) {
			next if $tempdef1 !~ /1_neg/;
			my $tempseq1 = $temp{$tempdef1};
			$tempseq1 =~ s/([A-Za-z0-9])\-*$/$1/;
			my ($tempbegseq1) = $temp{$tempdef1} =~ /^(\-*)[A-Za-z0-9]/;
			$tempbegseq1 = "" if not defined $tempbegseq1;
			foreach my $tempdef2 (sort keys %temp) {
				next if $tempdef2 !~ /3_con/;
				my ($tempbegseq2) = $temp{$tempdef2} =~ /^(\-*)[A-Za-z0-9]/;
				$tempbegseq2 = "" if not defined $tempbegseq2;
				$longbeg = length($tempbegseq1) > length($tempbegseq2) ? length($tempbegseq1) : length($tempbegseq2);
#				print "$tempdef, length=$longbeg\n";

				if ($longbeg > $beg0 and not defined $long{$tempdef1}{$tempdef2} and (($tempdef1 =~ /1_neg/ and $tempdef2 =~ /3_con/) or ($tempdef1 =~ /3_con/ and $tempdef2 =~ /1_neg/))) {
					$info .= "BEG0 ($LGN$beg0$N): longbeg=$LGN$longbeg$N, has to modify beg0 into";
					$long{$tempdef1}{$tempdef2} = $longbeg;
					$long{$tempdef2}{$tempdef1} = $longbeg;
					$beg0 = $longbeg + 10 < $junc ? $longbeg + 10 : $longbeg;
					$info .= " $LGN$beg0$N; detail: $beg0 = longest beg non-dash ($longbeg, + 10 if less than junc $junc)\n";
#					$info .= "$tempbegseq2\n";
				}


				my $tempseq2 = $temp{$tempdef2};
				$tempseq2 =~ s/([A-Za-z0-9])\-*$/$1/;
				$short = length($tempseq1) < length($tempseq2) ? length($tempseq1) : length($tempseq2);
				LOG($outBigLog, "HERE!\ntempdef=$tempdef\nlength=$short\n");
				next if defined $short{$tempdef1}{$tempdef2};
#				$end0 = $short if ($end0 > $short and $tempdef1 =~ /1_neg/ and $tempdef2 =~ /3_con/) or ($tempdef2 =~ /1_neg/ and $tempdef1 =~ /3_con/);
				if ($end0 > $short and not defined $short{$tempdef1}{$tempdef2} and (($tempdef1 =~ /1_neg/ and $tempdef2 =~ /3_con/) or ($tempdef1 =~ /3_con/ and $tempdef2 =~ /1_neg/))) {
					$short{$tempdef1}{$tempdef2} = $short;
					$short{$tempdef2}{$tempdef1} = $short;
#				if ($end0 >= $short and $tempdef2 =~ /3_con/ and $short eq length($tempseq1) and $tempdef1 =~ /1_neg/) or ($tempdef1 =~ /3_con/ and $short eq length($tempseq2) and $tempdef2 =~ /1_neg/ and $end0 >= $short)) {
					$info .= "END0 ($LGN$endposorig$N): short=$LGN$short$N, has to modify end0 into";
					my $infomod = "short $short-25 = " . ($short-25);
					$end0 = $short - 25 if $short eq length($tempseq1) and $tempdef1 =~ /1_neg/;
#					$end0 = $short - 25 if $short eq length($tempseq2) and $tempdef2 =~ /1_neg/;
					if ($end0 < $junc) {
						$end0 = $short - 20 if $short eq length($tempseq1) and $tempdef1 =~ /1_neg/;
#						$end0 = $short - 20 if $short eq length($tempseq2) and $tempdef2 =~ /1_neg/;
						$infomod = "short $short-20 = " . ($short-20);
					}
					if ($end0 < $junc) {
						$end0 = $short - 15 if $short eq length($tempseq1) and $tempdef1 =~ /1_neg/;
#						$end0 = $short - 15 if $short eq length($tempseq2) and $tempdef2 =~ /1_neg/;
						$infomod = "short $short-15 = " . ($short-15);
					}
					if ($end0 < $junc) {
						$end0 = $short - 10 if $short eq length($tempseq1) and $tempdef1 =~ /1_neg/;
#						$end0 = $short - 10 if $short eq length($tempseq2) and $tempdef2 =~ /1_neg/;
						$infomod = "short $short-10 = " . ($short-10);
					}
					if ($end0 < $junc) {
						$end0 = $short - 5 if $short eq length($tempseq1) and $tempdef1 =~ /1_neg/;
#						$end0 = $short - 5 if $short eq length($tempseq2) and $tempdef2 =~ /1_neg/;
						$infomod = "short $short-5 = " . ($short - 5);
					}
					if ($end0 < $junc) {
						$end0 = $short;
						$infomod = "shrot $short";
					}
					if ($end0 <= $junc) {
						$end0 = $junc + 5;
						if ($end0 > $short) {
							$end0 = $short;
							$infomod = "short $short as even junc $junc + 5 is bigger than short";
						}	
						else {
							$infomod = "junc $junc + 5 (short $short is too close)";
						}
					}
					$info = $info . " $infomod\n";
				}
			}
		}
#		print "short= $short\n";
	}

	if ($type eq "none") {
		$info .= "junc = $junc" if defined $juncpos;
		$info .= ", beg_pos = $beg0 (orig = $begposorig)" if $beg0 ne 0;
		$info .= ", end_pos = $end0 (orig = $endposorig)" if $end0 ne 0;
		$info .= ", mhbeg = $mhbeg (orig = $mhbegorig)";# if $mhbeg ne -1;
		$info .= ", mhend = $mhend (orig = $mhendorig)";# if $mhend ne -1;
	}
	my ($def, $seq) = ("INIT", "");
	my %seq;
	my %type;
	my $printnamewant = 0;
	my %print;
	my ($def1name, $def3name) = ("", "");
	my $printhead = "";
	my $printdef1 = "";
	my $printdef1mat = "";
	my $printdef2 = "";
	my $printdef2mat = "";
	my $printdef3 = "";
	my $whitespace = "";
	my %mut;
	my @mut;
	
	for (my $i = 0; $i < @line; $i++) {
#		print "$i. $line[$i]\n";
		my $line = $line[$i]; chomp($line);
		if ($line =~ /^>/) {
			if ($def ne "INIT") {
				$printnamewant = 1 if defined $namewant and $def =~ /$namewant/i;
				$seq{$def} = $seq;
#				if (defined $namewant) {
					if ($def =~ /1_neg/) {
						my @whitespace = (" ") x (length($def) - length($name));
						@whitespace = (" ") x (length($def)) if $type ne "none";
						$whitespace = join("", @whitespace);
						my @number; my @number10; my @flag;
						for (my $l = 0; $l < (1+length($seq)); $l++) {
							push(@number10, $l-1) if $l % 10 == 1;
							push(@number10, (" ") x (10 - length($l))) if $l % 10 == 1;
							push(@number, (($l % 10)));
						}
						@mut = (" ") x @number if @mut == 0;
						# colorize number with beg0 and end0 and junc
						if ($type eq "none") {#defined $junc and ($type eq "none")) {
							@flag = (" ") x (@number);
							my @seq = split("", $seq);
							my $seqind = 0; my $flag = 0;
							for (my $k = 0; $k < @number; $k++) {
								my $seqchunk = $seq[$k];
								$seqind = $k;
								if ($seqind eq $primer_pos_fix - 1) {
									$number[$k] = "$LRD$number[$k]$N";
									$flag[$k] = "P";
									$flag = 1 if $flag eq 0;
								}
								elsif ($seqind eq $adapter_pos_fix + 1) {
									$number[$k] = "$LRD$number[$k]$N";
									$flag[$k] = "E";
									$flag = 0;
								}
								elsif ($seqind eq $beg0) {
									$number[$k] = "$LRD$number[$k]$N";
									$flag  = 2;
									$flag[$k] = $seq[$k];
								}
								elsif (($seqind eq $junc)) {
									$number[$k] = "$LCY$number[$k]$N" if $seq[$k] ne "-";
									$flag[$k] = $seq[$k];
								}
								elsif ($seqind eq $junc+1) {
									$number[$k] = "$LGN$number[$k]$N" if $seq[$k] ne "-";
									$flag[$k] = $seq[$k];
								}
								elsif ($seqind eq $end0) {
									$flag[$k] = $seq[$k];
									$number[$k] = "$LPR$number[$k]$N";
									$flag = 1 if $flag eq 2;
								}
								elsif ($flag eq 1) {
									$flag[$k] = "-";
								}
								elsif ($flag eq 2) {
									$flag[$k] = $seq[$k];
									$number[$k] = "$YW$number[$k]$N";
								}
#								$seqind ++ if $seqchunk ne "-";
								$die = "Undefined flag at k=$k def=$def typr=$type\n" if not defined $flag[$k];
							}
							$whitespace = join("", @whitespace);
							$printhead .= join("", @whitespace) . "\t" . colorize(join("", @flag), $junc) . "\n";
						}
			
						$printhead .= join("", @whitespace) . "\t" . join("", @number10) . "\n";
						$printhead .= join("", @whitespace) . "\t" . join("", @number) . "\n";
					}
					if ($def =~ /1_neg/) {
						$def1name = $def;
						$printdef1 .= "$def\t" . colorize($seq, $junc, $reverse) . "\n" if $type eq "none";# if defined $namewant;# and $printnamewant eq 1;
						$printdef1 .= "$def\t" . colorize($seq) . "\n" if $type ne "none";# if defined $namewant;# and $printnamewant eq 1;
					}
					elsif ($def =~ /3_con/) {
						$def3name = $def;
						$printdef3 .= "$def\t" . colorize($seq, "NA", $reverse) . "\n";# if defined $namewant;# and $printnamewant eq 1;
					}
					else {
						$print{$def} = "$def\t" . colorize($seq, "NA", $reverse) . "\n";# if defined $namewant;# and $printnamewant eq 1;
					}
				}
#			}
			$seq = "";
			$def = $line; $def =~ s/>//;
		}
		else {
			$seq .= $line;
		}
	}
	$seq{$def} = $seq if $def ne "INIT";
#	die "juncpos = $junc, junc=$junc\n" if $junc ne 0;
#	if (defined $namewant and $def ne "INIT") {
	if ($def ne "INIT") {#efined $namewant and $def ne "INIT") {
		if ($def =~ /1_neg/) {
			$printdef1 .= "$def\t" . colorize($seq, $junc, $reverse) . "\n" if $type eq "none";
			$printdef1 .= "$def\t" . colorize($seq) . "\n" if $type ne "none";
		}
		elsif ($def =~ /3_con/) {
			$printdef3 .= "$def\t" . colorize($seq, "NA", $reverse) . "\n";
		}
		else {
			$print{$def} .= "$def\t" . colorize($seq, "NA", $reverse) . "\n";
		}
	}
	if ($type ne "none") {
		LOG($outBigLog, "\n\n$printhead$printdef1");# if $name eq $namewant;
		foreach my $def (sort keys %print) {
			LOG($outBigLog, $print{$def} . "\n\n");# if $name eq $namewant;
		}
	}
	if ($type eq "none") {
#		print "\n";
#		print $printhead;
		$print{head} = $printhead;
		if ($type eq "none") {
			$print{$def3name} = $printdef3;
			$print{$def1name} = $printdef1;
#			print $printdef3;
#			print $printdef1;
#			print $printdef2;
		}
		else {
			$print{$def3name} = $printdef3;
			$print{$def1name} = $printdef1;
#			print $printdef1;
#			print $printdef3;
#			print $printdef2;
		}
	}
	undef $def;
	undef $seq;

	# length of cons

	$reverseprint = defined $reverse ? $reverse : "rev undef";
#	if ($type =~ /jun/) {
#		LOG($outBigLog, $print{$def1name});
#	}

=comment
	if ($type eq "none") {
		LOG($outBigLog, $print{$def3name});
		foreach my $def (sort keys %print) {
			LOG($outBigLog, $print{$def}) if $def =~ /2_neg/;
		}
		LOG($outBigLog, $print{$def1name});
		foreach my $def (sort keys %print) {
			LOG($outBigLog, $print{$def}) if $def =~ /3_neg/;
		}
		foreach my $def (sort keys %print) {
			LOG($outBigLog, $print{$def}) if $def =~ /4_neg/;
		}
		foreach my $def (sort keys %print) {
			LOG($outBigLog, $print{$def}) if $def =~ /^5_/;
		}
		foreach my $def (sort keys %print) {
			LOG($outBigLog, $print{$def}) if $def =~ /^6_/;
		}
	}
=cut
#	die "aaa\n" if $type eq "none";

	$type{tot}{len} = 0 if not defined $type{tot}{len};
	$type{tot}{tot} = 0 if not defined $type{tot}{tot};
	if ($type eq "none") {
		$type{tot}{len} += 1 if $gap ne 0;
	}

	my $defind = 0;
	foreach my $def1 (sort keys %seq) {
		next if $def1 !~ /$leader/;
		foreach my $def2 (sort keys %seq) {
			next if $def1 eq $def2;
			#def is is always 1_neg
			if ($type eq "none") {
				next if $def2 !~ /(3_con|2_neg|4_neg|3_neg)/;
			}
			@{$mut{$def2}} = @mut;
			my @seq1 = split("", $seq{$def1});
			my @seq2 = split("", $seq{$def2});
			#print "def1 = $def1, def2=$def2\n" if defined $namewant;
			my $seqz1 = $seq{$def1}; 
			my $seqz2 = $seq{$def2};
			if (defined $reverse) {
				$seqz1 = mitochy::revcomp($seqz1);
				$seqz2 = mitochy::revcomp($seqz2);
			}
			$seqz1 =~ s/([ACGTN])[\-]+$/$1/i;
			$seqz2 =~ s/([ACGTN])[\-]+$/$1/i;
#			my ($countlen1) = $seq{$def1} =~ tr/ACGTNacgtn/ACGTNacgtn/;
#			my ($countlen2) = $seq{$def2} =~ tr/ACGTNacgtn/ACGTNacgtn/;
			my @seqz1 = split("", $seqz1);
			my @seqz2 = split("", $seqz2);
			my $shortestname;#= @seq1 < @seq2 ? $def1 : $def2;
			my $shortest;# = @seq1 < @seq2 ? @seq1 : @seq2;
			my $newend0 = $end0;
			if ($type ne "none") {
				$shortestname = @seqz1 > @seqz2 ? $def1 : $def2;
				$shortest     = @seqz1 > @seqz2 ? @seqz1 : @seqz2;
#				print "LONGEST = $shortestname: $shortest\n";
			}
			else {
				$shortestname = @seqz1 < @seqz2 ? $def1 : $def2;
				$shortest = @seqz1 < @seqz2 ? @seqz1 : @seqz2;
				$newend0 = $shortest if $end0 > $shortest;
#				print "SHORTEST = $shortestname: $shortest, end0=$end0, newend0=$newend0\n";
			}
			my ($mat, $mis, $ins, $del, $mh, $currmat, $maxmat, $maxpos, $len1, $len2) = (0,0,0,0,0,0,0,0,0,0);
			my $prev = "NA";
			my ($currmat2, $maxmat2, $maxpos2, $currpos2, $prev2, $currins2, $currdel2) = (0,0,0,"NA","NA",0,0);
         my $currpos = "NA";
			my $currins = "";
			my $currdel = "";
         my $ind1 = 0;
         my $ind2 = 0;
			my $lasti = 0;
			my %type2;
			#my $print3 = $whitespace . "\t";
			my $prevlenindels = 0;
			my $mhseq = "";
			my $seqfinal1 = "";
			my $seqfinal2 = "";
			my $hasmatch = 0;
			my $begdash = 0;
			my $dashes = 0;
			for (my $i = 0; $i < $shortest; $i++) {
				$lasti = $i;
				$mut{$def2}[$i] = " " if not defined $mut{$def2}[$i];
				my $seqchunk1 = $seq1[$i]; 
				my $seqchunk2 = $seq2[$i]; 
				$seqfinal1 .= $seqchunk1;
				$seqfinal2 .= $seqchunk2;
				$seqchunk1 = "-" if not defined $seqchunk1;
				$seqchunk2 = "-" if not defined $seqchunk2;
#				my $seqchunk1print = defined $seqchunk1 ? $seqchunk1 : "UNDEF";
#				my $seqchunk2print = defined $seqchunk2 ? $seqchunk2 : "UNDEF";
#				print "i=$i, ind1=$ind1, ind2=$ind2, len1=$len1, len2=$len2, seqchunk2=$seqchunk2print, undef seqchunk1 at name=$name\n" if not defined $seq1[$i];
#				print "i=$i, ind1=$ind1, ind2=$ind2, len1=$len1, len2=$len2, seqchunk1=$seqchunk1print, undef seqchunk2 at name=$name\n" if not defined $seq2[$i];
#				print "seq1=" . colorize(join("",@seq1)) . "\n" if not defined $seq1[$i];
#				print "seq2=" . colorize(join("",@seq2)) . "\n" if not defined $seq2[$i];
#				die if not defined $seq1[$i] or not defined $seq2[$i];
				$len2 ++ if $seqchunk2 ne "-" and $seqchunk2 ne "x";
				$len1 ++ if $seqchunk1 ne "-" and $seqchunk1 ne "x";
            $ind1 ++ if $seqchunk1 ne "-" and $seqchunk1 ne "x";
            $ind2 ++ if $seqchunk2 ne "-" and $seqchunk2 ne "x";
				
#				if ($name eq $namewant and $type eq "none") {
#					if ($beg0 != -1 and $end0 != -1 and $ind1 >= $beg0 and $ind1 <= $end0) {
#						$print3 .= colorize($seqchunk2);
#					# "$type: beg0=$beg0, end0 = $end0, ind1=$ind1\n";
#					}
#					else {
#						$print3 .= " ";
#					}
#				}

				if ($seqchunk1 eq $seqchunk2 and $def2 !~ /3_con/) {
#					$mut{$def2}[$i] = "|";
#					$mut{$def2}[$i] = "-" if $seqchunk1 eq "-";
				}
#				if ($mhbeg ne -1 and $ind1 >= $mhbeg and $ind1 <= $mhend and $seqchunk1 ne "-" and $type eq "none") {
				if ($mhbeg ne -1 and $i >= $mhbeg and $i <= $mhend and $type eq "none") {
#					LOG($outBigLog, "mhbeg=$mhbeg, ind=$ind1, chunk=$seqchunk1, mhend=$mhend\n") if $name eq $namewant;
					$mut{$def2}[$i] = "H" if $seqchunk1 ne "-" and $seqchunk1 eq $seqchunk2;
					$mhseq .= $seqchunk1 if $seqchunk1 ne "-" and $seqchunk1 eq $seqchunk2;
#					print "ind=$ind1, i=$i, mut $def2 [$i] = H\n";
					if ($seqchunk1 ne "-" and $seqchunk2 ne "-" and $def2 =~ /3_con/ and $def1 =~ /1_neg/ and $type eq "none") {
						$type{tot}{len} -= 1;
					}
					if ($i == $mhend  and $type eq "none" and $def2 =~ /3_con/ and $mhseq ne "") {
						$type{mh}{$mhseq} ++;
						push(@{$type2{mhpos}{$mhseq}}, $mhbeg);
						$muthash->{$mhbeg}{$readorient} = "mh_$mhseq";
						if ($seqchunk1 ne "-" and $seqchunk2 ne "-" and $def2 =~ /3_con/ and $def1 =~ /1_neg/ and $type eq "none") {
							$type{tot}{len} += 1;
						}
					}
				}

#				if (($ind1 >= $beg0 or $beg0 eq -1) and ($ind1 <= $end0 or $end0 eq -1)) {
				if (($i >= $beg0 or $beg0 eq -1) and ($i < $newend0 or $newend0 eq -1)) {
					if ($seqchunk1 ne "-" and $seqchunk2 ne "-" and $def2 =~ /3_con/ and $def1 =~ /1_neg/ and $type eq "none") {
#						die "END0 isn't newend0 ($end0 != $newend0\n" if $end0 ne $newend0;
						$type{tot}{len} ++;
					}
					#next if $type eq "none" and $def2 =~ /3_neg/ and $ind1 <= $junc;
					#next if $type eq "none" and $def2 =~ /2_neg/ and $ind1 > $junc;
					if ($def2 !~ /3_con/ and $hasmatch ne 1) {
						$currpos2 = $i;
						$begdash = $i;
#						$currmat2 ++;
					}
					if ($def2 =~ /3_con/ and $type eq "none") {
						$type{nuc}{$seqchunk2} ++;
						$type{tot}{tot} ++;
					}
					if ($seqchunk1 eq $seqchunk2 and $seqchunk1 eq "-") {
						$mut{$def2}[$i] = " ";
#						$currpos2 = $i if $hasmatch eq 1;
						$currmat2 ++ if $hasmatch eq 1;
						$dashes ++ if $hasmatch eq 1;
###						print "$i. $seqchunk1 $seqchunk2 $mut{$def2}[$i] currpos=$LCY$currpos$N,$LCY$currpos2$N prevlenindels=$LCY$prevlenindels$N prev=$LCY$prev$N currmat=$LCY$currmat$N,$LCY$currmat2$N maxmat=$LCY$maxmat$N,$LCY$maxmat2$N, begdash=$begdash, dashes=$dashes\n" if $type =~ /none/;
					}
					elsif ($seqchunk1 eq $seqchunk2 and $seqchunk1 ne "-") {
						$hasmatch = 1;
###						print "$i. $seqchunk1 $seqchunk2 $mut{$def2}[$i] currpos=$LCY$currpos$N,$LCY$currpos2$N prevlenindels=$LCY$prevlenindels$N prev=$LCY$prev$N currmat=$LCY$currmat$N,$LCY$currmat2$N maxmat=$LCY$maxmat$N,$LCY$maxmat2$N, begdash=$begdash, dashes=$dashes\n" if $type eq "none";
						if ($def2 =~ /3_con/ and $type eq "none") {
							if ($mut{$def2}[$i] eq "H") {
								$mh ++ if $def2 =~ /3_con/;
								my $nuc2nuc1 = "$seqchunk2\_$seqchunk1";
							}
							else {
								$mut{$def2}[$i] = "|";
								$mat ++;
								my $nuc2nuc1 = "$seqchunk2\_$seqchunk1";
								$type{mat}{$nuc2nuc1} ++;# if $def2 =~ /3_con/;
								#my $matpos = $ind1 - $junc;
								my $matpos = $i - $junc;
								push(@{$type2{matpos}{$nuc2nuc1}}, $matpos);
								$muthash->{$i}{$readorient} = "mat_$nuc2nuc1";
							}
						}
						else {
							$mut{$def2}[$i] = $mut{$def2}[$i] ne "H" ? "|" : "${LGN}H$N";
							$mat ++;
							my $nuc2nuc1 = "$seqchunk2\_$seqchunk1";
							#my $matpos = $ind1 - $junc;
							my $matpos = $i - $junc;
							push(@{$type2{matpos}{$nuc2nuc1}}, $matpos);
						}
						if ($currpos ne "NA" and $prevlenindels > 0 and $prevlenindels <= 5) {
	#						print "0 mis. prevlenindels = $prevlenindels, i = $i, maxpos = $maxpos, maxmat = $maxmat, currpos = $currpos, currmat = $currmat\n" if defined $namewant;
							$currmat += $prevlenindels;
						}
	  					$currmat ++;
	  					$currmat2 ++;
						$currpos2 = $i if $currpos eq "NA";
						$currpos = $ind1 if $currpos eq "NA";
						$prevlenindels = 0;
						$prev = "mat";
#						print "        -> currpos=$LGN$currpos$N,$LGN$currpos2$N prevlenindels=$LGN$prevlenindels$N prev=$LGN$prev$N currmat=$LGN$currmat$N,$LGN$currmat2$N maxmat=$LGN$maxmat$N,$LGN$maxmat2$N\n\n";
					}
					elsif ($seqchunk1 ne $seqchunk2 and $seqchunk1 ne "-" and $seqchunk2 ne "-") {
						$mis ++;
						if (defined $mut{$def2}[$i]) {
							$mut{$def2}[$i] = "m";
						}
						my $nuc2nuc1 = "$seqchunk2\_$seqchunk1";
						if ($def2 =~ /3_con/) {
							$type{mis}{$nuc2nuc1} ++;
							$muthash->{$i}{$readorient} = "mis_$nuc2nuc1";
						}
#						my $mispos = $ind1 - $junc;
						my $mispos = $i - $junc;
						push(@{$type2{mispos}{$nuc2nuc1}}, $mispos);
						if ($currpos ne "NA" and $prevlenindels > 0 and $prevlenindels <= 5) {
#							print "0 mis. prevlenindels = $prevlenindels, i = $i, maxpos = $maxpos, maxmat = $maxmat, currpos = $currpos, currmat = $currmat\n" if defined $namewant;
							$currmat += $prevlenindels;
						}
	  					$currmat ++;
	  					$currmat2 ++;
						$currpos2 = $i if $currpos eq "NA";
						$currpos = $ind1 if $currpos eq "NA";
						$prevlenindels = 0;
						$prev = "mis";
					}
					elsif ($currpos ne "NA") {
						if ($prevlenindels > 5) {
#							LOG($outBigLog, "\t1. prevlenindels = $prevlenindels, i = $i, maxpos = $maxpos, maxmat = $maxmat, currpos = $currpos, currmat = $currmat\n") if $name eq $namewant and $type ne "none";
#							LOG($outBigLog, "\t1. prevlenindels = $prevlenindels, i = $i, maxpos = $maxpos2, maxmat = $maxmat2, currpos = $currpos2, currmat = $currmat2\n") if $name eq $namewant and $type ne "none";
	   	            $maxpos2 = $currpos2 if $maxmat < $currmat;
							$maxmat2 = $currmat2 if $maxmat < $currmat;
	   	            $maxpos = $currpos if $maxmat < $currmat;
							$maxmat = $currmat if $maxmat < $currmat;
				         $currmat2 = 1;
				         $currpos2 = "NA";
	   	            $currmat = 1;
	   	            $currpos = "NA";
							$prevlenindels = 0;
						}
					}
					if ($seqchunk1 eq "-" and $seqchunk2 ne "-") {
						$del ++;
#					print "ind1 = $ind1, beg0=$beg0, end0=$end0, del=$seqchunk2, currdel = $currdel\n";
						$currdel .= $seqchunk2;
						$prevlenindels ++ if $prev eq "mat" or $prev eq "mis" or $prevlenindels > 0;
						$prev = "del";
						if (defined $mut{$def2}[$i]) {
							if ($mut{$def2}[$i] eq "I") {$mut{$def2}[$i] = "B";} else {$mut{$def2}[$i] = "D"};
						}
					}
					elsif ($currdel ne "") {
#						if (length($currdel) < 20) {
#							my $delpos = $ind1 - $junc;
							my $delpos = $i - $junc;
							push(@{$type2{delpos}{$currdel}}, $delpos);
							if ($def2 =~ /3_con/) {
								$type{del}{$currdel} ++;
								$muthash->{$i}{$readorient} = "del_$currdel";
							}
#						}
#						print "INS: $currdel $type{del}{$currdel}\n";
						$currdel = "";
					}
					if ($seqchunk1 ne "-" and $seqchunk2 eq "-") {
						$ins ++;
#					print "ind1 = $ind1, beg0=$beg0, end0=$end0, ins=$seqchunk2, currins = $currins\n";
						$currins .= $seqchunk1;
						$prevlenindels ++ if $prev eq "mat" or $prev eq "mis" or $prevlenindels > 0;
						$prev = "ins";
						if (defined $mut{$def2}[$i]) {
							if ($mut{$def2}[$i] eq "D") {$mut{$def2}[$i] = "B";} else {$mut{$def2}[$i] = "I"};
						}
					}
					elsif ($currins ne "") {
#							my $inspos = $ind1 - $junc;
							my $inspos = $i - $junc;
							push(@{$type2{inspos}{$currins}}, $inspos);
							if ($def2 =~ /3_con/) {
								$type{ins}{$currins} ++;
								$muthash->{$i}{$readorient} = "ins_$currins";
							}
#						}
#						print "DEL: $currins $type{ins}{$currins}\n";
						$currins = "";
					}
				}
#				if ($mhbeg ne -1 and $i >= 0 and $i <= $mhbeg and $def2 =~ /3_neg/ and $type eq "none" ) {
				if ($mhbeg ne -1 and $mhbeg ne 0 and $i == 0 and $def2 =~ /3_neg/ and $type eq "none") {
					$mut{$def2}[$i] = "$YW" . $mut{$def2}[$i] if $i == 0;
					$mut{$def2}[$i] = $mut{$def2}[$i] . "$N" if $i == $mhbeg;
				}
#				if ($mhend ne -1 and $i >= $mhend + 1 and $i <= @{$mut{$def2}} and $def2 =~ /2_neg/ and $type eq "none") {
				if ($mhend ne -1 and $mhend ne 0 and $i == $mhend + 1 and $def2 =~ /2_neg/ and $type eq "none") {
					$mut{$def2}[$i] = "$YW" . $mut{$def2}[$i] if $i == $mhend + 1;
					$mut{$def2}[$i] = $mut{$def2}[$i] . $N if $i == @{$mut{$def2}};
				}
			}
			if ($type eq "none") {
				$mut{$def2}[$end0] = "F";
				my $printdef4 = "";
#				if ($def2 =~ /3_con/) {
#					for (my
#				}
				$printdef4 .= $whitespace . "\t" . join("", @{$mut{$def2}}) . "\n";
				$print{$def2} = $print{$def2} . $printdef4 if $def2 =~ /2_neg/;
				$print{$def2} = $printdef4 . $print{$def2} if $def2 =~ /3_neg/;
				$print{$def2} = $print{$def2} . $printdef4 if $def2 =~ /3_con/;
			}
			if ($currpos ne "NA") {
#				LOG($outBigLog, "\t2. maxpos = $maxpos2, maxmat = $maxmat2, currpos = $currpos2, currmat = $currmat2\n") if $name eq $namewant and $type ne "none";
				$maxpos2 = $currpos2 if $maxmat < $currmat;
				$maxmat2 = $currmat2 if $maxmat < $currmat;
				$maxpos = $currpos if $maxmat < $currmat;
				$maxmat = $currmat if $maxmat < $currmat;
	         $currmat2 = 1;
	         $currpos2 = "NA";
	         $currmat = 1;
	         $currpos = "NA";
			}
			if ($currdel ne "") {
#				my $delpos = $ind1 - $junc;
				my $delpos = $shortest - $junc;
				push(@{$type2{delpos}{$currdel}}, $delpos);
				if ($def2 =~ /3_con/) {
					$type{del}{$currdel} ++;
				}
				$currdel = "";
			}
			if ($currins ne "") {
#				my $inspos = $ind1 - $junc;
				my $inspos = $shortest - $junc;
				push(@{$type2{inspos}{$currins}}, $inspos);
				if ($def2 =~ /3_con/) {
					$type{ins}{$currins} ++;
				}
				$currins = "";
			}
			if ($type eq "none" and $def2 =~ /3_con/) {
				foreach my $type1 (sort keys %type2) {
					foreach my $type2 (sort keys %{$type2{$type1}}) {
						my $mean = int(100*mitochy::mean(@{$type2{$type1}{$type2}}))/100;
						my $sd = int(mitochy::sd(@{$type2{$type1}{$type2}}))/100;
						if ($def2 =~ /3_con/) {
							$type{$type1}{$type2} = "$mean,$sd";
							#LOG($outBigLog, "type1=$type1, typ2=$type2, mean=$mean, sd=$sd, value=@{$type2{$type1}{$type2}}\n");
						}
					}
				}
			}
			if ($type ne "none") {
				LOG($outBigLog, "type=$type, def1=$LCY$def1$N,def2=$LGN$def2$N, $info, len1 = $len1, len2= $len2, lasti = $lasti\n");
			}
			my $percmaxmat = $len2 == 0 ? 0 : int(1000*$maxmat/$len2)/10;
			my $wantpos = -1;
			my $lenz1 = $junc;
#			my $lenz2 = $len1 - $junc;
			my $lenz2 = $shortest - $junc;
#			$maxmat = $maxmat2;
#			$maxpos = $maxpos2;
			if ($type ne "none") {
#				my $print3 = "\t3. maxmat = $maxmat, maxpos = $maxpos, perc = $percmaxmat, final pos =" . ($maxpos + $maxmat) . "\n" if $name eq $namewant and $type ne "none";
				my $print3 .= "\t3. maxmat = $maxmat2, maxpos = $maxpos2, perc = $percmaxmat, final pos =" . ($maxpos2 + $maxmat2) . "\n" if $name eq $namewant and $type ne "none";
				LOG($outBigLog, $print3) if $name eq $namewant and $type ne "none";
			}
			if ($type eq "primer") {
				$wantpos = $percmaxmat > 90 ? $maxpos2 + $maxmat2 : -1;
				my $primer_pos = $wantpos;
				my $orig_primer_pos = $primer_pos;
				$primer_pos = 25 if $primer_pos eq -1;
				my $beg_pos = $lenz1 - $lenMaxL >= $primer_pos ? $lenz1 - $lenMaxL : $lenz1 >= $primer_pos ? $primer_pos : ($junc - 10 > 0 and $primer_pos > $junc) ? $junc - 10 : ($primer_pos > $junc) ? $junc : -99;
				LOG($outBigLog, "\t${LGN}FINAL BEG_POS=$beg_pos and FINAL PRIMER_POS=$primer_pos$N. (orig_primer_pos = $YW$orig_primer_pos$N, perc mat=$percmaxmat\%, thus we use $YW$primer_pos$N. Junc = $LGN$junc$N, so beg_pos = $LCY$beg_pos$N, begdash=$begdash, dashes=$dashes)\n\n");# if $name eq $namewant;
				return($wantpos, $primer_pos, $beg_pos, $begdash, $dashes);
			}
			if ($type eq "adapter") {
#			   my ($seqadp0) = $seq{$def1} =~ /^(.*)$seq{$def}2/;
#				if (defined $seqadp0) {
#					$wantpos = length($seqadp0);
#				}
#				else {
#					$wantpos = $percmaxmat > 90 ? $maxpos2 : -1;
#				}
				$wantpos = $percmaxmat > 90 ? $maxpos2 : -1;
				my $adapter_pos = $wantpos;
				my $orig_adapter_pos = $adapter_pos;
				LOG($outBigLog, "adapter wantpos = $wantpos\n");
				LOG($outBigLog, "\t${LGN}FINAL END_POS=-99 and FINAL ADAPTER_POS=$adapter_pos$N. (orig_adapter_pos = $YW$orig_adapter_pos$N, perc mat=$percmaxmat\%, thus we use $YW$adapter_pos$N. len1=$lenz1, len2=$lenz2, lenMaxL=$lenMaxL, lenMaxR=$lenMaxR, junc = $LGN$junc$N, so end pos = $LPR-99$N, begdash=$begdash, dashes=$dashes)$N\n\n") if $name eq $namewant and $adapter_pos ne -1 and $adapter_pos < $lenz1;
				return(-99, -99, -99, $begdash, $dashes) if $adapter_pos < $lenz1 and $adapter_pos ne -1;
				$adapter_pos = $lenz1 + $lenz2 - 25 if $adapter_pos eq -1;
				$adapter_pos = $lenz1 + $lenz2 - $adapter_pos; #25 or 20
				my $end_pos = $lenz2 >= $lenMaxR + $adapter_pos ? $lenz1 + $lenMaxR : $lenz2 >= $adapter_pos ? $lenz1 + ($lenz2 - $adapter_pos) : -99;
				LOG($outBigLog, "\t${LGN}FINAL END_POS=$end_pos and FINAL ADAPTER_POS=$adapter_pos$N. (orig_adapter_pos = $YW$orig_adapter_pos$N, perc mat=$percmaxmat\%, thus we use $YW$adapter_pos$N. len1=$lenz1, len2=$lenz2, lenMaxL=$lenMaxL, lenMaxR=$lenMaxR, junc = $LGN$junc$N, so end pos = $LPR$end_pos$N, begdash=$begdash, dashes=$dashes)$N\n\n");# if $name eq $namewant;
				return($wantpos, ($len1 - $adapter_pos), $end_pos, $begdash, $dashes);
			}
			if ($type eq "junc") {
#				print "$def1\t$seqfinal1\n$def2\t$seqfinal2\ni=$lasti/$shortest, ind1=$ind1, ind2=$ind2";
				$wantpos = $percmaxmat > 80 ? $maxpos2 + $maxmat2 : -1;
				return($wantpos, $begdash, $dashes);
			}
		}
	}
	my $NA = $name eq $namewant ? undef : "NA";
	my $print_aln = "";
	if ($type =~ /none/) {
		$print_aln .= "\n\n";
		$print_aln .= "------------------------------------------------\n";
		$print_aln .= "1. CONSENSUS (CON) vs. READ (SEQ):\n";
		$print_aln .= "$YW$name$N (${LGN}R$readorient$N)\n";
		$print_aln .= "------------------------------------------------\n";
		$print_aln .= "type{tot}len} = $type{tot}{len}\n";
		$print_aln .= "type{tot}tot} = $type{tot}{tot}\n\n";
		$print_aln .= "$print{head}";
		foreach my $defz (sort keys %print) {
			if ($defz =~ /3_con/) {
				my $printz = $print{$defz};
				$printz =~ s/3_con/  CON/;
				$printz =~ s/_$name//g;
				$print_aln .= "$printz";
			}
		}
		foreach my $defz (sort keys %print) {
			if ($defz =~ /1_neg/) {
				my $printz = $print{$defz};
				$printz =~ s/1_neg/  SEQ/;
				$printz =~ s/_$name//g;
				$print_aln .= "$printz";
			}
		}
		$print_aln .= "\n\n";
		$print_aln .= "------------------------------------------------\n";
		$print_aln .= "2. IgM BAIT (IgM) vs. READ (SEQ) vs FAR JOINED (FAR):\n";
		$print_aln .= "$YW$name$N (${LGN}R$readorient$N)\n";
		$print_aln .= "------------------------------------------------\n\n";
		$print_aln .= "$print{head}";
		foreach my $defz (sort keys %print) {
			if ($defz =~ /2_neg/) {
				my $printz = $print{$defz};
				$printz =~ s/2_neg/  IgM/;
				$printz =~ s/_$name//g;
				$print_aln .= "$printz";
			}
		}
		foreach my $defz (sort keys %print) {
			if ($defz =~ /1_neg/) {
				my $printz = $print{$defz};
				$printz =~ s/1_neg/  SEQ/;
				$printz =~ s/_$name//g;
				$print_aln .= "$printz";
			}
		}
		foreach my $defz (sort keys %print) {
			if ($defz =~ /3_neg/) {
				my $printz = $print{$defz};
				$printz =~ s/3_neg/  FAR/;
				$printz =~ s/_$name//g;
				$print_aln .= "$printz";
			}
		}
		$print_aln .= "\n";
		$print_aln .= "------------------------------------------------$N\n" if $readorient eq "1";
		$print_aln .= "$YW<------------------------------------------------>$N\n\n" if $readorient eq "	2";
	}
	$print_aln .= "$YW$name$N\n" . $info . "\n\n";
	LOG($outLog, $print_aln, "NA");
	LOG($outBigLog, $print_aln);
	die if not defined $type{tot}{tot};
	LOG($outBigLog, $die) if $die ne 0;
	return($muthash, \%type, $beg0, $end0);
}

1;

__END__
Qname   JuncID  Rname   Junction        Strand  Rstart  Rend    B_Rname B_Rstart        B_Rend  B_Strand        B_Qstart        B_Qend  Qstart  Qend    Qlen    B_Cigar Cigar   Seq     J_Seq      Barcode unaligned       baitonly        uncut   misprimed       freqcut largegap        mapqual breaksite       sequential      repeatseq       duplicate
M02034:489:000000000-CYYL8:1:1101:16291:1246    1       chr12   114663366       -1      114663157       114663366       chr12   114664820       114664909       -1      7       96      96310      343     90M     34M5I7M1X7M1X35M1X3M1X1M1X2M1X21M1X2M1X14M1X6M1X27M1X4M1X1M1X2M3X14M1X7M1X5M    ACAGTGCACACAAAGACTCTGGACCTCTCCGAAACCAGGCACCGCAAATGGTAAGCCAGAGGCAGCCACAGCTGTGGCTGCTGCTCTTAAAGCTTGCTGAGCTGGGGTGAGCTGGGGTGAGCTGAGCTGAGCTGGGGTGAGATGAGCTGTGCTGGGGTGAGCTGAGCTGGGGTGAGCTGAGCTGATCTGAGCTGGGCTGAGCTGAGCTGAGCTGAGGTGGGCTGGGGTGAGCTGGGCTGAGGTGAGCTGGGGTGAGCTGGGGTGAGCTGGGCTGGGGTGAGCTGAGCTGAGCTGGGGTGAGCTGAGGTGACCACGCGTGCTCTACAACTTTCGTAAGATCGTA     

MUSCLE (3.8) multiple sequence alignment


2_neg_M02034:489:000000000-CYYL8      -----------------GTTGGA------GGATATGGGGGA--------GGCGAGCATGA
