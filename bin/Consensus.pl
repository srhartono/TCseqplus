#!/usr/bin/perl

use strict; use warnings; use mitochy; use Getopt::Std; use FAlite; use Cwd qw(abs_path); use File::Basename qw(dirname);
use vars qw($opt_v $opt_i);
getopts("vi:");

BEGIN {
   my $libPath = dirname(dirname abs_path $0) . '/lib';
   push(@INC, $libPath);
   #print "\n- Pushed $libPath into perl lib path INC\n";
}

use Mutationlib;

my ($input1) = ($opt_i);
die "\nusage: $YW$0$N -i $CY<input1>$N\n\n" unless defined $opt_i and -e $opt_i;

my ($folder1, $fileName1) = mitochy::getFilename($input1, "folderfull");
my %data;
my $inputFA = $input1;
open (my $outLog, ">", "$inputFA.LOG") or die;
open (my $outBigLog, ">", "$inputFA.BIGLOG") or die;
my ($gene) = $inputFA =~ /(IgM|IgG1|IgG3)/;
die "Cannot parse gene (IgM or IgG1 or IgG3) from inputFA $LCY$inputFA$N\n" if not defined $gene;
open (my $inFASTA, "<", $inputFA) or die "Cannot read from $inputFA: $!\n";
my $fasta = new FAlite($inFASTA);
while (my $entry = $fasta->nextEntry()) {
	my $seq = $entry->seq;
	my $def = $entry->def;
	$def =~ s/^>//;
	my ($chr, $beg, $end, $genome) = $def =~ /^(.+)\:(\d+)\-(\d+)_(B6|S129|MG1G3)$/;
	die "Cannot parse genome (B6 or S129) from def=$LGN$def$N\nfastaFile=$LCY$inputFA$N\n\n" if not defined $genome;
	$def = "$gene\_$genome";
#	$seq = mitochy::revcomp($seq);# if $genome eq "B6";
	push(@{$data{def}}, $def);
	push(@{$data{seq}}, $seq);
	push(@{$data{chr}}, $chr);
	push(@{$data{beg}}, $beg);
	push(@{$data{end}}, $end);
	push(@{$data{genome}}, $genome);
	push(@{$data{gene}}, $gene);
}
close $inFASTA;

my $def1 = $data{def}[0];
my $seq1 = $data{seq}[0];
my $def2 = $data{def}[1];
my $seq2 = $data{seq}[1];

#my $seq1revcomp = mitochy::revcomp($seq1);
#my $seq2revcomp = mitochy::revcomp($seq2);

LOG($outBigLog,  colorize($seq1) . "\n","NA");
LOG($outBigLog,  colorize($seq2) . "\n","NA");

my ($len0, $len1) = (0,0);
my ($prev0, $prev1) = (0,0);
my ($curr0, $curr1) = (0,0);
my ($seq10, $seq11, $seq20, $seq21) = (0,0,0,0);
if ($seq1 =~ /^\-/) {
	my ($edgetemp) = $seq1 =~ /^([\-]+)[A-Z]/;
	my $lentemp = (not defined $edgetemp) ? 0 : length($edgetemp);
	$prev0 = $lentemp;
	$len0 = $lentemp;
	$seq10 = 1;
}
if ($seq1 =~ /\-$/) {
	my ($edgetemp) = $seq1 =~ /[A-Z]([\-]+)$/;
	my $lentemp = (not defined $edgetemp) ? 0 : length($edgetemp);
	$prev1 = $lentemp;
	$len1 = $lentemp;
	$seq11 = 1;
}
if ($seq2 =~ /^\-/) {
	my ($edgetemp) = $seq2 =~ /^([\-]+)[A-Z]/;
	my $lentemp = (not defined $edgetemp) ? 0 : length($edgetemp);
	$curr0 = $lentemp;
	$len0 = $lentemp;
	$seq20 = 1;
}
if ($seq2 =~ /\-$/) {
	my ($edgetemp) = $seq2 =~ /[A-Z]([\-]+)$/;
	my $lentemp = (not defined $edgetemp) ? 0 : length($edgetemp);
	$curr1 = $lentemp;
	$len1 = $lentemp;
	$seq21 = 1;
}

my ($chr1, $beg1, $end1, $genome1) = ($data{chr}[0], $data{beg}[0], $data{end}[0], $data{genome}[0]);
my ($chr2, $beg2, $end2, $genome2) = ($data{chr}[1], $data{beg}[1], $data{end}[1], $data{genome}[1]);
my $strand1 = "+";
$strand1 = "-" if $genome1 eq "B6";
my $strand2 = "+";
$strand2 = "-" if $genome2 eq "B6";



LOG($outBigLog, "$chr1\t$beg1\t$end1\t$gene\_$genome1\t0\t$strand1\n");
LOG($outBigLog, "$chr2\t$beg2\t$end2\t$gene\_$genome2\t0\t$strand2\n");
if ($len0 > 0) {
	
	if ($seq1 =~ /^\-{$len0}[A-Z]/) {
	}
	else {
		$beg1 += $len0 if $strand1 eq "+";
		$end1 -= $len0 if $strand1 eq "-";
	}

	if ($seq2 =~ /^\-{$len0}[A-Z]/) {
	}
	else {
		$beg2 += $len0 if $strand2 eq "+";
		$end2 -= $len0 if $strand2 eq "-";
	}

	$seq1 =~ s/^.{$len0}//;
	$seq2 =~ s/^.{$len0}//;
}
if ($len1 > 0) {

	if ($seq1 =~ /[A-Z]\-{$len1}$/) {
	}
	else {
		$end1 -= $len1 if $strand2 eq "+";
		$beg1 += $len1 if $strand2 eq "-";
	}

	if ($seq2 =~ /[A-Z]\-{$len1}$/) {
	}
	else {
		$end2 -= $len1 if $strand2 eq "+";
		$beg2 += $len1 if $strand2 eq "-";
	}

	$seq1 =~ s/.{$len1}$//;
	$seq2 =~ s/.{$len1}$//;
}


LOG($outBigLog, colorize($seq1) . "\n","NA");
LOG($outBigLog, colorize($seq2) . "\n","NA");
$strand1 = "+";
$strand1 = "-" if $genome1 eq "B6";
$strand1 = "-" if $genome1 eq "MG1G3";
$strand2 = "+";
$strand2 = "-" if $genome2 eq "B6";
$strand2 = "-" if $genome2 eq "MG1G3";
LOG($outBigLog, "$chr1\t$beg1\t$end1\t$gene\_$genome1\t0\t$strand1\n");
LOG($outBigLog, "$chr2\t$beg2\t$end2\t$gene\_$genome2\t0\t$strand2\n");
open (my $outBedFixB6, ">", "$inputFA.B6.fix.bed") or die;
open (my $outBedFixS129, ">", "$inputFA.S129.fix.bed") or die;
print $outBedFixB6 "$chr1\t$beg1\t$end1\t$gene\_$genome1\t" . ($end1 - $beg1) . "\t$strand1\n";
print $outBedFixS129 "$chr2\t$beg2\t$end2\t$gene\_$genome2\t" . ($end2 - $beg2) . "\t$strand2\n";
close $outBedFixB6;
close $outBedFixS129;
system("fastaFromBed -fi ~/Bowtie2_indexes/B6IGH/B6IGH.fa -bed $inputFA.B6.fix.bed -fo $inputFA.B6.fix.fa -s");
system("fastaFromBed -fi ~/Bowtie2_indexes/S129IGH/S129IGH.fa -bed $inputFA.S129.fix.bed -fo $inputFA.S129.fix.fa -s");
#system("cat $inputFA.B6.fix.fa $inputFA.S129.fix.fa > $inputFA.B6S129.fix.fa.temp");
#system("muscle -in $inputFA.B6S129.fix.fa.temp -clw > $inputFA.B6S129.fix.fa.aln 2>&1");

LOG($outBigLog, "\n$LCY\nless $inputFA.B6S129.fix.fa.aln\n$N\n");

#$seq1 = mitochy::revcomp($seq1) if $strand1 eq "-";
#$seq2 = mitochy::revcomp($seq2) if $strand2 eq "-";

#my ($chr1, $beg1, $end1, $genome1) = ($data{chr}[0], $data{beg}[0], $data{end}[0], $data{genome}[0]);
#my ($chr2, $beg2, $end2, $genome2) = ($data{chr}[1], $data{beg}[1], $data{end}[1], $data{genome}[1]);
#my ($edge0) = $seq1 =~ /^([\-]+)[A-Z]/;
#my $len0 = (not defined $edge0) ? 0 : length($edge0);
my @seq1 = split("", $seq1);
my @seq2 = split("", $seq2);
#$seq1 = substr($seq1,0,100);
#$seq2 = substr($seq2,0,100);
$seq1 .= "N";
$seq2 .= "N";
my @res;
$res[0] = ">3_con_NAME";
$res[1] = $seq1;#mitochy::revcomp($seq1);
$res[2] = ">1_neg_NAME";
$res[3] = $seq2;#mitochy::revcomp($seq2);
$res[4] = ">2_neg_NAME";
$res[5] = $seq1;#mitochy::revcomp($seq1);
$res[6] = ">3_neg_NAME";
$res[7] = $seq2;#mitochy::revcomp($seq2);
my $beg_pos = 0;
my $end_pos = length($seq1)-1; #100;
my $junc_pos = int($end_pos/2); #50;
my $primer_pos_fix = $beg_pos;
my $adapter_pos_fix = $end_pos;
my $mhbeg = -1;
my $mhend = -1;
my $gap = 0;
my $lenMaxL = $junc_pos;
my $lenMaxR = ($end_pos - $junc_pos);
my $namewant = "NAME";
my $name = "NAME";
LOG($outBigLog, "beg_pos=$LCY$beg_pos$N, end_pos=$LCY$end_pos$N, junc_pos=$LGN$junc_pos$N, lenMaxL=$LGN$junc_pos$N, lenMaxR=$LGN$lenMaxR$N\n");
#   ($muthash, $typehash2, $beg_pos2, $end_pos2) = parse_aln(\@resTN2, "1_", "none,$beg_pos,$end_pos,$juncpos,$primer_pos_fix,$adapter_pos_fix,$DATA->{$name}{mhbeg},$DATA->{$name}{mhend}", $muthash, 2, $lenMaxL, $lenMaxR, ", 0, $outBigLog, $outLog, $namewant);
my ($muthash, $typehash1, $beg_pos2, $end_pos2);
($muthash, $typehash1, $beg_pos2, $end_pos2) = parse_aln(\@res, "1_", "none,$beg_pos,$end_pos,$junc_pos,$primer_pos_fix,$adapter_pos_fix,$mhbeg,$mhend", $muthash, 1, $lenMaxL, $lenMaxR, $name, $gap, $outBigLog, $outLog, $namewant);
my $seq1_nodash = $seq1; $seq1_nodash =~ s/\-//g;
my $seq2_nodash = $seq2; $seq2_nodash =~ s/\-//g;
#for (my $i = 0; $i < @

=comment
open (my $outMut, ">", "$inputFA.mut") or die;
foreach my $pos (sort {$a <=> $b} keys %{$muthash}) {
	foreach my $R (sort keys %{$muthash->{$pos}}) {
		my $mutcurr = $muthash->{$pos}{$R};
#		next if $muthash->{$pos}{$R} =~ /^mat/ and $pos < (keys %{$muthash}) - 300 and $pos > 300 and $pos !~ /^\d*3\d\d$/ and $pos !~ /^10[0-9][0-9][0-9]$/;
#		next if $mutcurr =~ /^mat/ and ($pos < 8000 or $pos >= 12000);
		next if $mutcurr =~ /^mat/ and (($pos > 1000 and $pos < 4000) or ($pos >= 6000 and $pos < length($seq1) - 1000));
		next if $R ne 1;
		my $posbeg1 = $pos;
		my $posend1 = $pos + 1;
		my $posbeg2 = $pos;
		my $posend2 = $pos + 1;
		if ($mutcurr =~ /^(mis|mat)_\w+_\w+$/) {
			my ($mut, $nuc1, $nuc2) = split("_", $mutcurr);
		}
		elsif ($mutcurr =~ /^(ins|del)_\w+$/) {
			my ($mut, $nuc1) = split("_", $mutcurr);
			my $nucl = length($nuc1) - 1;
			if ($mutcurr =~ /ins/) {
				$posbeg2 = $pos;
				$posend2 = $pos + $nucl + 1;
			}
			elsif ($mutcurr == /del/) {
				$posbeg1 = 
			}
			$posbeg2 = $pos;
			$posend2 = $pos + $nucl + 1 if $mutcurr =~ /ins_/;
#			$posfix = $pos - $nucl;
		}
		my ($posbeg1_nodash) = fix_pos($posbeg1, $seq1, $seq1_nodash);
		my ($posend1_nodash) = fix_pos($posend1, $seq1, $seq1_nodash);
		my ($posbeg2_nodash) = fix_pos($posbeg2, $seq2, $seq2_nodash);
		my ($posend2_nodash) = fix_pos($posend2, $seq2, $seq2_nodash);
		my ($pos1_nodash) = fix_pos($pos, $seq1, $seq1_nodash);
		my ($pos2_nodash) = fix_pos($pos, $seq2, $seq2_nodash);
#		my $pos1fixpos = length($seq1_nodash) - $pos1_nodash;
#		my $pos2fixpos = length($seq2_nodash) - $pos2_nodash;
		#my $genomeposfix1_nodash = $strand1 eq "+" ? $beg1 + $posfix1_nodash : $end1 - $posfix1_nodash;
		#my $genomeposfix2_nodash = $strand2 eq "+" ? $beg2 + $posfix2_nodash : $end2 - $posfix2_nodash;
		#my $genomepos1_nodash = $strand1 eq "+" ? $beg1 + $pos1_nodash : $end1 - $pos1_nodash;
		#my $genomepos2_nodash = $strand2 eq "+" ? $beg2 + $pos2_nodash : $end2 - $pos2_nodash;
		#my $outprint = "$pos\t$posfix\t$"
		my $outprint = "$pos\t$mutcurr\t$posbeg1\t$posbeg1_nodash\t$chr1\n";#$genomeposfix1_nodash\t$chr2\t$genomeposfix2_nodash\tR$R\t$mutcurr\t$pos1_nodash\t$pos2_nodash\t$posfix1_nodash\t$posfix2_nodash\t$genomepos1_nodash\t$genomepos2_nodash\n";
#		LOG($outBigLog, "$pos $posfix $R $mutcurr\n");
#		LOG($outLog, "$pos $posfix $R $mutcurr\n");
		LOG($outBigLog, $outprint);
		LOG($outLog, $outprint,"NA");
		print $outMut $outprint;
#		print $outprint;
#		print "$pos $posfix $R $mutcurr\n";
	}
}
close $outMut;
=cut
my @flag = (" ") x @seq1;
my @count = (0) x @seq1;
my @nuc = ("") x @seq1;
my @beg;
my @end;
my ($chunk1, $chunk2, $flagchunk, $numchunk, $countchunk) = ("","","","","");
#my $ind1 = $strand1 eq "+" ? $beg1 : $end1;
#my $ind2 = $strand2 eq "+" ? $beg2 : $end2;
#my $add1 = $strand1 eq "+" ? 1 : -1;
#my $add2 = $strand2 eq "+" ? 1 : -1;
my $ind1 = $beg1;
my $ind2 = $beg2;
my $add1 = 1;
my $add2 = 1;
open (my $out1, ">", "$fileName1.out") or die "Cannot write to $fileName1.out: $!\n";
open (my $outmut1, ">", "$fileName1.mut") or die "Cannot write to $fileName1.out: $!\n";
#@seq1 = split("", mitochy::revcomp($seq1));
#@seq2 = split("", mitochy::revcomp($seq2));
#my $mutchunk = "";
my $lasti = 0;
my $lastmut = -1;
my @nucz;
for (my $i = 0; $i < @seq1; $i++) {
	my $mutcurr = $muthash->{$i}{1};
	my $ind = $i + 1;
	my $nuc1 = $seq1[$i];
	my $nuc2 = $seq2[$i];
	$ind1 += $add1 if $nuc1 ne "-";
	$ind2 += $add2 if $nuc2 ne "-";
	$mutcurr = "UNKNOWN" if not defined $mutcurr;
	print $outmut1 "$i\t$ind1\t$ind2\t$mutcurr\n";
	if ($nuc1 eq $nuc2 and $nuc1 ne "-") {
		$lasti = $i;
		$flag[$i] = " ";
		$lastmut = "M";
		$nucz[$i] = "mat_$nuc1";
		$count[$i] ++;
	}
	elsif ($nuc1 eq "-" and $nuc2 ne "-") {
		if ($lastmut eq "I") {
			$nucz[$lasti] .= $nuc2;
			$count[$lasti] ++;
		}
		else {
			$lasti = $i;
			$lastmut = "I";
			$flag[$i] = "I";
			$nucz[$i] = "ins_$nuc2";
			$count[$i] ++;
		}
	}
	elsif ($nuc1 ne "-" and $nuc2 eq "-") {
		if ($lastmut eq "D") {
			$nucz[$lasti] .= $nuc1;
			$count[$lasti] ++;
		}
		else {
			$lasti = $i;
			$lastmut = "D";
			$flag[$i] = "D";
			$nucz[$i] = "del_$nuc1";
			$count[$i] ++;
		}
	}
	elsif ($nuc1 ne $nuc2 and $nuc1 ne "-" and $nuc2 ne "-") {
		$flag[$i] = "m";
		$lasti = $i;
		$lastmut = "m";
		$nucz[$i] = "mis_$nuc1\_$nuc2";
		$count[$i] ++;
	}
	else {
		die "can't determine flag at i=$LCY$i$N nuc1=$LGN$nuc1$N, nuc2=$LPR$nuc2$N\n";
	}
}
my ($indbeg1, $indbeg2) = (0,0);
my ($indend1, $indend2) = (0,0);
my $printmut;
for (my $i = 0; $i < @seq1; $i++) {
	my $ind = $i;#+1
	my $nuc1 = $seq1[$i];
	my $nuc2 = $seq2[$i];
	my $mutcurr = $nucz[$i];
	if (defined $mutcurr and $mutcurr ne "") {
#		$indbeg1 ++ if $nuc1 ne "-";#= $add1 if $nuc1 ne "-";
#		$indbeg2 ++ if $nuc2 ne "-";#+= $add2 if $nuc2 ne "-";
		$mutcurr = $nucz[$i]; # if $flag[$i] =~ /^(D|I)$/;#"UNKNOWN" if not defined $mutcurr;
		if ($mutcurr =~ /ins/) {
			$indend1 += 0;
			$indend2 += $count[$i];
		}
		elsif ($mutcurr =~ /del/) {
			$indend2 += 0;
			$indend1 += $count[$i];
		}
		else {
			$indend1 ++;
			$indend2 ++;
		}
		$printmut .= "$i\t$indbeg1\t$indend1\t$indbeg2\t$indend2\t$mutcurr\tcount=$count[$i]\n";
		if ($mutcurr =~ /ins/) {
			$indbeg1 += 0;
			$indbeg2 = $indend2;
		}
		elsif ($mutcurr =~ /del/) {
			$indbeg2 += 0;
			$indbeg1 = $indend1;
		}
		else {
			$indbeg1 ++;
			$indbeg2 ++;
		}
	}
=comment
	if ($nuc1 eq $nuc2 and $nuc1 ne "-") {
#		$lasti = $i;
#		$flag[$i] = " ";
#		$lastmut = "M";
	}
	elsif ($nuc1 eq "-" and $nuc2 ne "-") {
#		if ($lastmut eq "I") {
#			$count[$lasti] ++;
#		}
#		else {
#			$lasti = $i;
#			$lastmut = "I";
#			$flag[$i] = "I";
#			$count[$i] ++;
#		}
#	}
	elsif ($nuc1 ne "-" and $nuc2 eq "-") {
		if ($lastmut eq "D") {
			$count[$lasti] ++;
		}
		else {
			$lasti = $i;
			$lastmut = "D";
			$flag[$i] = "D";
			$count[$i] ++;
		}
	}
	elsif ($nuc1 ne $nuc2 and $nuc1 ne "-" and $nuc2 ne "-") {
		$flag[$i] = "m";
		$count[$i] ++;
		$lasti = $i;
		$lastmut = "m";
	}
	else {
		die "can't determine flag at i=$LCY$i$N nuc1=$LGN$nuc1$N, nuc2=$LPR$nuc2$N\n";
	}
=cut
	$chunk1 .= $nuc1;
	$chunk2 .= $nuc2;
	$countchunk .= $count[$i];
	$countchunk =~ s/0/ /g;
	$flagchunk .= $flag[$i];
	$numchunk .= $ind % 10;
	if ($i ne 0 and ($ind+1) % 100 == 0) {
		print $out1 ">" . ($ind-99) . "\n$numchunk\n" . colorize($chunk1) . "\n" . colorize($chunk2) . "\n$flagchunk\n$countchunk\n\n";
		($chunk1, $chunk2, $flagchunk, $numchunk, $countchunk) = ("","","","","");
		print $out1 $printmut if defined $printmut;
		undef $printmut;
	}
#	print $out1 "$i,$nuc1,$nuc2,
#	if ($edge0 eq 0 and $nuc1 eq "-" or $nuc2 eq "-") {
#	}
}
if ($chunk1 ne "") {
	print $out1 ">" . (int(@seq1/100) * 100) . "\n$numchunk\n" . colorize($chunk1) . "\n" . colorize($chunk2) . "\n$flagchunk\n$countchunk\n\n";
	($chunk1, $chunk2, $flagchunk, $numchunk, $countchunk) = ("","","","","");
	print $out1 $printmut if defined $printmut;
	undef $printmut;
}
close $out1;

=comment
my $flag = join("", @flag);
LOG($outBigLog, "\n\nB6   " . colorize($seq1) . "\nflag $flag\nS129 " . colorize($seq2) . "\n\n");
LOG($outBigLog, "prevlen = $prev0, $prev1\ncurrlen = $curr0, $curr1\nusedlen = $len0, $len1\n");

system
=cut
