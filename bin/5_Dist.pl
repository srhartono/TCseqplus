#!/usr/bin/perl
##
use strict; use warnings; use mitochy; use Getopt::Std; use FAlite; use Cwd qw(abs_path); use File::Basename qw(dirname);
use vars qw($opt_v $opt_i $opt_m $opt_n);
getopts("vm:i:n:");

BEGIN {
   my $libPath = dirname(dirname abs_path $0) . '/lib';
   push(@INC, $libPath);
   #print "\n- Pushed $libPath into perl lib path INC\n";
}

use Mutationlib;

my ($inputFile) = $opt_i;
die "\nusage: $YW$0$N -i $CY<input1>$N\n\n" unless defined $opt_i and -e $opt_i;

my $namewant = $opt_n;

my $outLogFile = "$inputFile.5_Distlog";
open (my $outLog, ">", $outLogFile) or die;
LOG($outLog, "\n\n-i = $LCY$opt_i$N\n\nEXAMPLE:\n\n");
parse_inputFile_short($inputFile, "$inputFile.juncdist", $namewant, $outLog);
#my ($junc1pos, $junc2pos, $chr, $beg, $end, $strand, $outLog) = @_;
close $outLog;

sub junc2check {
   my ($beg1, $end1, $beg2, $end2, $strand1, $strand2, $junc1, $junc2, $linecount, $line, $outLog) = @_;
   my $junc2check = $strand2 eq "-" ? $end2 : $beg2;
	DIELOG($outLog, "Died ${LCY}$inputFile${N} at linecount=$linecount: junc2check $junc2check isn't the same as junc2 $junc2
   junc1 ($junc1) = strand1 ($strand1) eq + ? end1 ($end1) : beg1 ($beg1)
   junc2check ($junc2check) = strand2 ($strand2) eq + ? beg2 ($beg2) : end2 ($end2)
   junc2 = $junc2\n\nline:\n$line\n\n") if $junc2check ne $junc2;
}


sub parse_inputFile_short {
	my ($inputFile, $outDistFile, $namewant, $outLog) = @_;
	my ($log, $die) = ("", 0);
	my $linecount = 0;
	my $total_printed = 0;
	my ($data, @def);
	my ($folder1, $fileName1) = mitochy::getFilename($inputFile, "folderfull");

	open (my $in1, "<", $inputFile) or return($data, $log . "\nCannot read from $LCY$inputFile$N: $!\n", 1);
	open (my $outDist, ">", "$outDistFile") or return($data, $log . "\nCannot write into $LCY$outDistFile$N: $!\n", 1);
	print $outDist "name,sampleID,sample,treat,isotype,strand,junc,juncperc,juncdist,junclen,type\n";

	my ($sampleID) = $fileName1 =~ /^\w*([WDRS][123](cut)?)_.+$/;
	DIELOG($outLog, "Can't parse sampleID from fileName1=$LCY$fileName1$N\n") if not defined $sampleID;
	my ($sample, $treat) = get_sampletype($sampleID);
	my ($total_line_inputFile) = `wc -l $inputFile` =~ /^(\d+)/;

	$log .= " -> Parsing $LGN$total_line_inputFile$N lines from inputFile=$LCY$inputFile$N\n";
	while (my $line = <$in1>) {
		chomp($line);
		next if $line =~ /^#/;
		$log .= "  - parsed ${LGN}$linecount${N} / ${LCY}$total_line_inputFile${N} on "  . date() . "\n" if $linecount % 1e3 == 0;
		my @arr = split("\t", $line);
		if ($line =~ /^Qname\t/) {
			@def = @arr;
#			$linecount ++;
			next;
		}
		$linecount ++;
		LOG($outLog, "#$LCY$line$N\n") if $linecount eq 1;
		my ($name, $juncid, $chr2, $junc2, $strand2, $beg2, $end2, $chr1, $beg1, $end1, $strand1, $begT1, $endT1, $begT2, $endT2, $lenTN, $cigarQ1, $cigarQ2, $seqTN, $defseqjunc, $defbarcode, $unaligned, $baitonly, $uncut, $misprimed, $freqcut, $largegap, $mapqual, $breaksite, $sequential, $repeatseq, $duplicate) = @arr;

		# change strand from -1/1 into -/+
		$strand1 = $strand1 eq "-1" ? "-" : $strand1 eq "1" ? "+" : $strand1 eq "-" ? "-" : $strand1 eq "+" ? "+" : return($data, $log . "\nCan't define strand1 (strand1=$strand1) not 1 or -1!\n\n$line\n\n", 1);
		$strand2 = $strand2 eq "-1" ? "-" : $strand2 eq "1" ? "+" : $strand2 eq "-" ? "-" : $strand2 eq "+" ? "+" : return($data, $log . "\nCan't define strand2 (strand2=$strand2) not 1 or -1!\n\n$line\n\n", 1);

		$lenTN = length($seqTN) if $lenTN !~ /^[0-9]+$/;
		if (defined $opt_n) {
			next if $namewant ne $name;
		}
		# save order so output can be sorted with same order as inputFile
		my ($order) = scalar(keys(%{$data})) + 1;
		$namewant = $name if not defined $namewant and not defined $opt_n;# and $strand2 eq "+"; #to print example use 1st order name, or use -n if defined

		my $junc1 = $strand1 eq "+" ? $end1 : $beg1;
		#check if junc2 is the same as end2 (-) or beg2 (+)
		junc2check($beg1, $end1, $beg2, $end2, $strand1, $strand2, $junc1, $junc2, $linecount, $line, $outLog);

#		$data->{$name}{junc1pos} = $junc1;
#		$data->{$name}{junc2pos} = $junc2;
		
		LOG($outLog, "#-> beg1=$beg1-$end1, beg2=$beg2-$end2, junc1=$junc1, junc2=$junc2, chr2=$chr2, strand2=$strand2\n") if $name eq $namewant;
		my ($igtype, $junc1dist, $junc2dist, $junc1len, $junc2len) = get_igtype($junc1, $junc2, $chr2, $beg2, $end2, $strand2, $outLog);
#																	   my ($junc1pos, $junc2pos, $chr, $beg, $end, $strand, $outLog) = @_;
		my $junc1perc = $junc1len == 0 ? $junc1dist : int($junc1dist / $junc1len * 1000+0.5)/10;
		my $junc2perc = $junc2len == 0 ? $junc2dist : int($junc2dist / $junc2len * 1000+0.5)/10;
		$junc1perc = $junc1perc >= 100 ? 100 : $junc1perc <= 0 ? 0 : $junc1perc;
		$junc2perc = $junc2perc >= 100 ? 100 : $junc2perc <= 0 ? 0 : $junc2perc;
		push(@{$data->{$igtype}{junc1perc}}, $junc1perc);
		push(@{$data->{$igtype}{junc2perc}}, $junc2perc);
		push(@{$data->{$igtype}{junc1dist}}, $junc1dist);
		push(@{$data->{$igtype}{junc2dist}}, $junc2dist);
		push(@{$data->{$igtype}{junc1len}}, $junc1len);
		push(@{$data->{$igtype}{junc2len}}, $junc2len);
		print $outDist "$name,$sampleID,$sample,$treat,$igtype,$strand1,junc1,$junc1perc,$junc1dist,$junc1len,BAIT\n";
		print $outDist "$name,$sampleID,$sample,$treat,$igtype,$strand2,junc2,$junc2perc,$junc2dist,$junc2len,PREY\n";
#		print $outDist "$name,$sampleID,$sample,$treat,$igtype,$strand1,$strand2,junc1,$junc2,$junc1perc,$junc2perc,$junc1dist,$junc2dist,$junc1len,$junc2len\n";
		$total_printed ++;
		LOG($outLog, "#-> $name,$junc1,$junc2,$igtype,$junc1dist,$junc2dist,$junc1len,$junc2len\n") if $name eq $namewant;
		last if defined $opt_n and $name eq $namewant;
	}
	LOG($outLog, "\n#${LGN}Done$N! Printed $LGN$total_printed$N / $YW$linecount$N total reads\n#Output=\n#$LCY$outDistFile$N\n");
	close $outDist;

	open (my $outDistALL, ">", "$outDistFile.ALL") or return($data, $log . "\nCannot write into $LCY$outDistFile.ALL$N: $!\n", 1);
	print $outDistALL "sampleID,sample,treat,isotype,strand,total_read,junc1percmean,junc2percmean,junc1distmean,junc2distmean,junc1lenmean,junc2lenmean,junc1percse,junc2percse,junc1distse,junc2distse,junc1lense,junc2lense\n";
	LOG($outLog, "#sampleID\tsample\ttreat\tisotype\ttotal_read\tjunc1percmean\tjunc2percmean\tjunc1distmean\tjunc2distmean\tjunc1lenmean\tjunc2lenmean\tjunc1percse\tjunc2percse\tjunc1distse\tjunc2distse\tjunc1lense\tjunc2lense\n");
	foreach my $igtype (sort keys %{$data}) {
		foreach my $strand2 (sort keys %{$data->{$igtype}}) {
			my $total_read = @{$data->{$igtype}{junc1perc}};
			my $junc1percmean = myformat(mitochy::mean(@{$data->{$strand2}{$igtype}{junc1perc}}));
			my $junc1percse = myformat(mitochy::se(@{$data->{$strand2}{$igtype}{junc1perc}}));
			my $junc2percmean = myformat(mitochy::mean(@{$data->{$strand2}{$igtype}{junc2perc}}));
			my $junc2percse = myformat(mitochy::se(@{$data->{$strand2}{$igtype}{junc2perc}}));
			my $junc1distmean = myformat(mitochy::mean(@{$data->{$strand2}{$igtype}{junc1dist}}));
			my $junc1distse = myformat(mitochy::se(@{$data->{$strand2}{$igtype}{junc1dist}}));
			my $junc2distmean = myformat(mitochy::mean(@{$data->{$strand2}{$igtype}{junc2dist}}));
			my $junc2distse = myformat(mitochy::se(@{$data->{$strand2}{$igtype}{junc2dist}}));
			my $junc1lenmean = myformat(mitochy::mean(@{$data->{$strand2}{$igtype}{junc1len}}));
			my $junc1lense = myformat(mitochy::se(@{$data->{$strand2}{$igtype}{junc1len}}));
			my $junc2lenmean = myformat(mitochy::mean(@{$data->{$strand2}{$igtype}{junc2len}}));
			my $junc2lense = myformat(mitochy::se(@{$data->{$strand2}{$igtype}{junc2len}})); 
			$junc1percse = 0 if $junc1percse < 1e-6;
			$junc2percse = 0 if $junc2percse < 1e-6;
			$junc1distse = 0 if $junc1distse < 1e-6;
			$junc2distse = 0 if $junc2distse < 1e-6;
			$junc1lense = 0 if $junc1lense < 1e-6;
			$junc2lense = 0 if $junc2lense < 1e-6;
			print $outDistALL "$sampleID,$sample,$treat,$igtype,$strand2,$total_read,$junc1percmean,$junc2percmean,$junc1distmean,$junc2distmean,$junc1lenmean,$junc2lenmean,$junc1percse,$junc2percse,$junc1distse,$junc2distse,$junc1lense,$junc2lense\n";
			LOG($outLog, "$sampleID\t$sample\t$treat\t$igtype\t$strand2,$total_read\t$junc1percmean\t$junc2percmean\t$junc1distmean\t$junc2distmean\t$junc1lenmean\t$junc2lenmean\t$junc1percse\t$junc2percse\t$junc1distse\t$junc2distse\t$junc1lense\t$junc2lense\n");
		}
	}
	close $outDistALL;
}

sub myformat {
   my ($val) = @_;
   $val = $val == 0 ? 0 : $val > 10 ? int(10*$val+0.5)/10 : $val > 1 ? int(10*$val+0.5)/10 : $val > 0.001 ? int(10000*$val+0.5)/10000 : $val;
   if ($val =~ /^\-?0\.0+[1-9]\d*$/) {
      my ($val2) = $val =~ /^(\-?0\.0+[1-9]\d?)\d*$/;
      print "Can't reformat val ($val) into val2\n" if not defined $val2;
      $val = $val2 if defined $val2;
   }
   if ($val =~ /^\-?\d+\.?\d*e\-?\d+\.?\d*$/) {
      my ($val2, $val3) = $val =~ /^(\-?\d+\.?\d?\d?)\d*(e.+)$/;
      print "Can't reformat val ($val) into val2\n" if not defined $val2;
      $val = $val2 . $val3 if defined $val2 and defined $val3;
   }
   if ($val =~ /^\-?[1-9]\d*\.[1-9]\d*$/) {
      my ($val2) = $val =~ /^(.+\.\d\d?)\d*$/;
      print "Can't reformat val ($val) into val2\n" if not defined $val2;
      $val = $val2 if defined $val2;
   }
   return $val;
}


__END__

		# defined junc1 pos (IgM) which is beg1 if strand is +, or end1 if -
		# PREY - strand [ beg2<==PREY==]end2(junc2) ------ (junc1)beg1<==BAIT==]end1 ] BAIT - strand
		# PREY + strand [ end2<==PREY==]beg2(junc2) ------ (junc1)beg1<==BAIT==]end1 ] BAIT - strand
		# PREY - strand [ beg2<==PREY==]end2(junc2) ------ (junc1)end1<==BAIT==]beg1 ] BAIT + strand
		# PREY + strand [ end2<==PREY==]beg2(junc2) ------ (junc1)end1<==BAIT==]beg1 ] BAIT + strand

		my $junc1 = $strand1 eq "+" ? $end1 : $beg1;
		#check if junc2 is the same as end2 (-) or beg2 (+)
		junc2check($beg1, $end1, $beg2, $end2, $strand1, $strand2, $junc1, $junc2, $linecount, $line, $outLog);

		$data->{$name}{junc1pos} = $junc1;
		$data->{$name}{junc2pos} = $junc2;
		$data->{$name}{line} .= "\n$dashhead\n>${LGN}INPUT LINE:${N}\n";
		my $myprintcigarseq = "";
		for (my $i = 0; $i < @def; $i++) {
			if ($def[$i] =~ /^(B_Cigar|Cigar|Seq)$/i) {
				my ($space) = join("", (" ") x (15 - length($def[$i])));
				my $val = $arr[$i];
				$val = colorize($val) if $def[$i] =~ /Seq/i;
				$myprintcigarseq .= "$def[$i]$space = $val\n";
			}
			else {
				$data->{$name}{line} .= "$def[$i]=${LGN}$arr[$i]${N}";
				$data->{$name}{line} .= ","  if ($i != 0 and $i % 5 != 0 and $i != @def - 1);
				$data->{$name}{line} .= "\n" if ($i != 0 and $i % 5 == 0 or ($i == @def - 1));
			}
		}
		$data->{$name}{line} .= $myprintcigarseq . "\n$dashhead\n";
		# TCseq uses both 1 based for coordinate. fastaFromBed uses 0/1 (1 beg 0 end). So turn ALL TCseq end coordinate into 1 based (beg-1)
		# starts from 1, 
		$data->{$name}{order}   = $order;
		$data->{$name}{chr1}    = $chr1;
		$data->{$name}{chr2}    = $chr2;
		$data->{$name}{beg1}    = $beg1 - 1;
		$data->{$name}{beg2}    = $beg2 - 1;
		$data->{$name}{end1}    = $end1;
		$data->{$name}{end2}    = $end2; #turn into 1 based as fastaFromBed is 1 based
		$data->{$name}{begJ1}   = $junc1 - 1;     #begJ1 = junc1 pos, already 0 based.
		$data->{$name}{endJ1}   = $junc1; #endjunc = begJ1 + 1 to make it 1 based
		$data->{$name}{begJ2}   = $junc2 - 1;     #again already 0 based
		$data->{$name}{endJ2}   = $junc2; #endjunc = begJ2 + 1 to make it 1 based
		$data->{$name}{strand1} = $strand1;
		$data->{$name}{strand2} = $strand2;
		$data->{$name}{len1}    = $data->{$name}{end1} - $data->{$name}{beg1};
		$data->{$name}{len2}    = $data->{$name}{end2} - $data->{$name}{beg2};
		$data->{$name}{lenC1L}  = $data->{$name}{end1} - $data->{$name}{beg1};
		$data->{$name}{lenC2R}  = $data->{$name}{end2} - $data->{$name}{beg2};
#		$log .= "$data->{$name}{beg1}-$data->{$name}{end1}, len=$data->{$name}{lenC1L}\n$data->{$name}{beg2}-$data->{$name}{end2}, len=$data->{$name}{lenC2R}\n") if $name eq $namewant;
		$data->{$name}{cigarQ1}  = $cigarQ1;	
		$data->{$name}{cigarQ2}  = $cigarQ2;	
		$data->{$name}{begT1}   = $begT1 - 1;	   #7, 0 based = 6, 1 based = 7
		$data->{$name}{endT1}   = $endT1;	 #156, 0 based = 155, 1 based = 156
		$data->{$name}{lenT1}   = $data->{$name}{endT1} - $data->{$name}{begT1}; #1 based 156-7+1 = 150, 0/1based 6to 156 = 156-6 = 150, same as cigar=150M
		$data->{$name}{begT2}   = $begT2 - 1;	       #154; 0 based 153, 1 based 154
		$data->{$name}{endT2}   = $endT2 ;     #183; 0 based 182, 1 based 183
		$data->{$name}{lenT2}   = $data->{$name}{endT2} - $data->{$name}{begT2}; #1 based = 183 - 154 + 1 = 30, 0/1based 153-183 = 30M, same as cigar=30M
		$data->{$name}{seqTN}    = $seqTN;	
		$data->{$name}{begTN}    = 1 - 1; #orig=1, mod=0; 0based=0, 1 based=1
		$data->{$name}{endTN}    = $data->{$name}{begTN} + $lenTN - 1; #orig=217; 0 based 216, 1 based 217.
		$data->{$name}{lenTN}    = $lenTN;     #217; 0 based 216-0+1 = 217; 0/1 based 217-0 = 217
		$data->{$name}{mhbeg}   = -1;
		$data->{$name}{mhend}   = -1;
		$data->{$name}{mhpos}   = -1;
		$data->{$name}{mhseq}   = "";
		$data->{$name}{gap}     = 0;
		#if Qstart < B_Qend, then remove (B_Qend-Qstart) bp from the Qstart (beg2 if strand is +, or end2 if strand is -)

		my $beg1orig = $data->{$name}{beg1};
		my $end1orig = $data->{$name}{end1};
		my $len1orig = $data->{$name}{len1};
		my $beg2orig = $data->{$name}{beg2};
		my $end2orig = $data->{$name}{end2};
		my $len2orig = $data->{$name}{len2};

		my $beg2orend2extra = $strand2 eq "+" ? $beg2orig : $end2orig;
		my $origextra = "";
		my $printmhseqlog = "\n\n--------------- MICROHOMOLOGY! --------------\n";

		if ($data->{$name}{begT2} - $data->{$name}{endT1} < 0) {
#			my $mhlenz = $data->{$name}{begT2} - $data->{$name}{endT1}; # <0 if mh, > 0 if gap
			my $mhlenz = abs($data->{$name}{begT2} - $data->{$name}{endT1});
			$data->{$name}{beg2} += $mhlenz if $strand2 eq "+"; #mhlenz < 0 if mh
			$data->{$name}{end2} -= $mhlenz if $strand2 eq "-";
#			$data->{$name}{beg2} -= $mhlenz if $strand2 eq "+";
#			$data->{$name}{end2} += $mhlenz if $strand2 eq "-";
			$data->{$name}{lenC2R}  = $data->{$name}{end2} - $data->{$name}{beg2};
			my $beg2extra = $data->{$name}{beg2};
			my $end2extra = $data->{$name}{end2};
			my $len2extra = $data->{$name}{end2} - $data->{$name}{beg2};
#			   $mhlenz = "+ $mhlenz" if $strand2 eq "+";
#			   $mhlenz = "- $mhlenz" if $strand2 eq "-";
			$origextra .= "<original beg2> ($LCY$beg2orig$N) + <length of MH> ($YW$mhlenz$N bp) = $LCY$beg2extra$N (len $YW$len2extra$N bp) " if $strand2 eq "+";
			$origextra .= "<original end2> ($LGN$end2orig$N) - <length of MH> ($YW$mhlenz$N bp) = $LGN$end2extra$N (len $YW$len2extra$N bp) " if $strand2 eq "-";
			$beg2orend2extra = $strand2 eq "+" ? $beg2extra : $end2extra;
			$data->{$name}{mhbeg} = $data->{$name}{begT2};
			$data->{$name}{mhend} = $data->{$name}{endT1};
			$data->{$name}{mhseq} = substr($data->{$name}{seqTN}, $data->{$name}{begT2}, ($data->{$name}{endT1} - $data->{$name}{begT2}));
			my @seqTN = split("", $data->{$name}{seqTN});
			my $mhseq = "";
			$printmhseqlog .= "\n\n-------------------------------------\n";
			$printmhseqlog = ">${LGN}MICROHOMOLOGY${N}:\n\nnum    = ";
			for (my $i = 0; $i < @seqTN; $i++) {
				my $ind = (($i+1) % 10);
				if ($i == $data->{$name}{begT2}) {
					$ind = "${LCY}$ind";
					$seqTN[$i] = "${LCY}$seqTN[$i]";
					$mhseq .= $seqTN[$i]
				}
				elsif ($i == $data->{$name}{endT1}-1) { #-1 as it's 1 based
					$ind = "$ind${N}";
					$seqTN[$i] = "$seqTN[$i]${N}";
					$mhseq .= $seqTN[$i];
				}
				elsif ($i >= $data->{$name}{begT2} and $i <= $data->{$name}{endT1}-1) {
					$mhseq .= $seqTN[$i];
				}
				else {
					$seqTN[$i] = ".";
				}
				$printmhseqlog .= $ind;
			}
			$data->{$name}{endT1} = $data->{$name}{begT2};
			$data->{$name}{mhpos} = $data->{$name}{mhend} - $data->{$name}{mhbeg};
			$printmhseqlog .= "\nseqT   = $data->{$name}{seqTN}\nseqTmh = " . join("", @seqTN) . "\n\nMHbeg after mod from base1/1 to 0/1 = $data->{$name}{mhbeg}-$data->{$name}{mhend}, mhseq saved=$data->{$name}{mhseq}, mhseq=$mhseq, mhpos = $data->{$name}{mhpos}\n";
			$data->{$name}{mhlen} = $data->{$name}{mhend} - $data->{$name}{mhbeg};
		}
		else {
			$printmhseqlog .= "\nseqT   = seqT2\n" . colorize($data->{$name}{seqTN}) . "\n\nRead=${LCY}$name${N} has no MH!\n";
			$printmhseqlog .= "\n\n-------------------------------------\n\n";
			$data->{$name}{gap}   = $data->{$name}{begT2} - $data->{$name}{endT1};
#			$data->{$name}{endT1} = $data->{$name}{begT2};

		}
		
		$data->{$name}{line} .= $printmhseqlog;

		# BAIT - strand: ADD PREY length bp of sequence beg1 - (end2-beg2)

		# BAIT (-): B (beg1-(end2-beg2)) <==end2-beg2==]              =J= (junc1/beg1) <====BAIT=====] (end1) E
		# or
		# BAIT (+): E (end1+(end2-beg2)) <==end2-beg2==]              =J= (junc1/end1) <====BAIT=====] (beg1) B

		# PREY (-): B            (beg2) <====PREY=====] (end2/junc2) =J=              <==end1-beg1==] (end2+(end1-beg1)) E
		# or 
		# PREY (+): E            (end2) <====PREY=====] (beg2/junc2) =J=              <==end1-beg1==] (beg2-(end1-beg1)) B

		$data->{$name}{beg1_mod} = $data->{$name}{strand1} eq "-" ? $data->{$name}{beg1} - $data->{$name}{lenC2R} :  $data->{$name}{beg1};
		$data->{$name}{end1_mod} = $data->{$name}{strand1} eq "-" ? $data->{$name}{end1}                        : $data->{$name}{end1} + $data->{$name}{lenC2R};
		$data->{$name}{beg2_mod} = $data->{$name}{strand2} eq "-" ? $data->{$name}{beg2}                        : $data->{$name}{beg2} - $data->{$name}{lenC1L};
		$data->{$name}{end2_mod} = $data->{$name}{strand2} eq "-" ? $data->{$name}{end2} + $data->{$name}{lenC1L}: $data->{$name}{end2};
		$data->{$name}{len1_mod} = $data->{$name}{end1_mod} - $data->{$name}{beg1_mod};
		$data->{$name}{len2_mod} = $data->{$name}{end2_mod} - $data->{$name}{beg2_mod};

		print $outBed "$data->{$name}{chr1}\t$data->{$name}{beg1_mod}\t$data->{$name}{end1_mod}\t$name\t0\t$data->{$name}{strand1}\n";
		print $outBed "$data->{$name}{chr2}\t$data->{$name}{beg2_mod}\t$data->{$name}{end2_mod}\t$name\t0\t$data->{$name}{strand2}\n";

		$data->{$name}{begorig1} = $data->{$name}{beg1};
		$data->{$name}{endorig1} = $data->{$name}{end1};
		$data->{$name}{begorig2} = $data->{$name}{beg2};
		$data->{$name}{endorig2} = $data->{$name}{end2};
		$data->{$name}{beg1} = $data->{$name}{beg1_mod};
		$data->{$name}{end1} = $data->{$name}{end1_mod};
		$data->{$name}{beg2} = $data->{$name}{beg2_mod};
		$data->{$name}{end2} = $data->{$name}{end2_mod};
		$data->{$name}{len1} = $data->{$name}{end1} - $data->{$name}{beg1};
		$data->{$name}{len2} = $data->{$name}{end2} - $data->{$name}{beg2};
		my $lenC1L = $data->{$name}{lenC1L};
		my $lenC2R = $data->{$name}{lenC2R};
		($data->{$name}{lenC2L}, $data->{$name}{lenC1R}) = ($lenC1L, $lenC2R);
		($data->{$name}{lenC1L}, $data->{$name}{lenC2R}) = ($data->{$name}{len1_mod} - $data->{$name}{lenC1R}, $data->{$name}{len2_mod} - $data->{$name}{lenC2L});

		if ($name eq $namewant) {
			$log .= "\n -> ${LGN}BED COORDINATE EXAMPLE:${N} (0/1 based)\n";
			$log .= "  - LINE        = $line\n";
			$log .= "  - NAME        = $LGN$name$N\n";
			$log .= "  - REF1 FINAL  = $data->{$name}{chr1}:$LCY$data->{$name}{beg1}$N-$data->{$name}{end1} (strand $LCY$data->{$name}{strand1}$N, len $YW$data->{$name}{len1}$N bp)\n";
			$log .= "  - REF2 FINAL  = $data->{$name}{chr2}:$data->{$name}{beg2}-$LGN$data->{$name}{end2}$N (strand $LCY$data->{$name}{strand2}$N, len $YW$data->{$name}{len2}$N bp)\n";
			$log .= "  - REF1 L/R    = L=$YW$data->{$name}{lenC1L}$N, R=$YW$data->{$name}{lenC1R}$N\n";
			$log .= "  - REF2 L/R    = L=$YW$data->{$name}{lenC2L}$N, R=$YW$data->{$name}{lenC2R}$N\n";
			$log .= "  - SEQ TLX LEN = " . length($data->{$name}{seqTN}) . "\n";
			$log .= "  - MH          = $LGN$data->{$name}{mhbeg}$N - $LGN$data->{$name}{mhend}$N, bp from junction=$LCY$data->{$name}{mhpos}$N\n";
			$log .= "  - GAP         = $LGN$data->{$name}{gap}$N\n";
			$log .= "\n";
			$log .= "  - REF1 MOD INFO:\n";
			$log .= "     - [ORIGINAL] = $data->{$name}{chr1}:$LCY$beg1orig$N-$LGN$end1orig$N (strand $LCY$data->{$name}{strand1}$N, len $YW$len1orig$N bp)\n";
			$log .= "     - ${LCY}[MOD]$N      = <original beg1> ($LCY$beg1orig$N) - <orig/MH length of REF2> ($YW$data->{$name}{lenC2R}$N) = $LCY$data->{$name}{beg1}$N\n";
			$log .= "     - [FINAL]    = $data->{$name}{chr1}:$LCY$data->{$name}{beg1}$N-$data->{$name}{end1} (strand $LCY$data->{$name}{strand1}$N, len $YW$data->{$name}{len1}$N bp)\n";
			$log .= "\n";
			$log .= "  - REF2 MOD INFO:\n";
			$log .= "     - [ORIGINAL] = $data->{$name}{chr2}:$LCY$beg2orig$N-$LGN$end2orig$N (strand $LCY$data->{$name}{strand2}$N, len $YW$len2orig$N bp)\n";
			$log .= "     - ${LCY}[MH]$N       = $origextra\n" if $origextra ne "";
			$log .= "     - ${LCY}[MOD]$N      = <orig/MH beg2>  ($LCY$beg2orend2extra$N) + <original length of REF1> ($YW$len1orig$N) = $LCY$data->{$name}{beg2}$N\n" if $strand2 eq "+";
			$log .= "     - ${LCY}[MOD]$N      = <orig/MH end2>  ($LGN$beg2orend2extra$N) - <original length of REF1> ($YW$len1orig$N) = $LGN$data->{$name}{end2}$N\n" if $strand2 eq "-";
# finalcoor1 = beg1=$LCY$data->{$name}{beg1}$N, end1=$data->{$name}{end1} (len1=$YW$data->{$name}{len1}$N bp)\n";
			$log .= "     - [FINAL]    = $data->{$name}{chr2}:$data->{$name}{beg2}-$LGN$data->{$name}{end2}$N (strand $LCY$data->{$name}{strand2}$N, len $YW$data->{$name}{len2}$N bp)\n";

		}

		if (defined $opt_n and defined $data->{$namewant}) {
			last;
		}
		elsif ($debug_num ne 0 and $data->{$name}{order} >= $debug_num) {
			last;
		}
	}
	close $in1;
	close $outBed;

	$log .= "\n ------------\n ${LGN}Done parsing inputFile$N\n\n Parsed $LGN" . (keys %{$data}) . "$N reads (linecount = $LGN$linecount$N)\n\n";
	return($data, $namewant, $log, $die);
}


__END__
__END__
   my ($junc1pos, $junc2pos, $chr, $beg, $end, $strand, $outLog) = @_;

return($junctype, $junc1dist, $junc2dist);

