#!/usr/bin/perl

use strict; use warnings; use mitochy; use Getopt::Std; use FAlite; use Cwd qw(abs_path); use File::Basename qw(dirname);
use vars qw($opt_v $opt_i $opt_e $opt_o);
getopts("vi:e:o:");

BEGIN {
   my $libPath = dirname(dirname abs_path $0) . '/lib';
   push(@INC, $libPath);
        print "\n- Pushed $libPath into perl lib path INC\n";
}

use Mutationlib;

die "
Usage: $YW$0$N ${LGN}[Optional: -e ext -o outFile]$N $LCY-i <folder containing .all>$N

${LGN}Optional$N:

-e extension (e.g. 50.50.final.tsv.all) [.all]
-o output file [inputfolder/results.csv]

" unless defined $opt_i and -d $opt_i;

# inputFolder
my $inputFolder = $opt_i;
$inputFolder = getFullpath($inputFolder);
$inputFolder =~ s/\/+$//;
$inputFolder =~ s/\/\.$//;

# ext
my $ext     = defined $opt_e ? "$opt_e" : "*.final.tsv.all";
if ($ext =~ /^.+.tlx\..+.final.tsv.all$/) {
	($ext) = $ext =~ /^.+.tlx\.(.+)$/;
}

# inputFiles
my @inputFiles = <$inputFolder/*$ext>;
my $inputFilesPrint = ""; #for log
for (my $i = 0; $i < @inputFiles; $i++) {
	my $ind = $i+1 . ".";
	$ind .= $ind < 10 ? " " : "";
	$inputFiles[$i] =~ s/^\.\///;
	my ($foldertemp, $inputFilestemp) = getFilename($inputFiles[$i], "folderfull");
	my ($sampleID) = $inputFilestemp =~ /^([A-Z][0-9])_/;
	($sampleID) = $inputFilestemp =~ /([A-Z][0-9](cut)?)_result_filter/ if not defined $sampleID;
	$sampleID = "UNKNOWN" if not defined $sampleID;
	my ($lenMaxL, $lenMaxR) = $inputFilestemp =~ /^.+.tlx.(\d+)\.(\d+)\.final.tsv.all$/;
	($lenMaxL) = ("UNKNOWN") if not defined ($lenMaxL);
	($lenMaxR) = ("UNKNOWN") if not defined ($lenMaxR);
	$inputFilesPrint .= "\t$LGN$ind$N <inputFolder>/$LCY$inputFilestemp$N (sampleID $LGN$sampleID$N, lenMaxL $LGN$lenMaxL$N, lenMaxR $LGN$lenMaxR$N)\n";
	die "\n\n$LRD ERROR!!! $N: Cannot parse sampleID ($LGN$sampleID$N) (e.g. W1 or D3) from inputFile $LCY$inputFilestemp$N\n\n" if not defined $sampleID or $sampleID eq "UNKNOWN";
	#die "SmapleID = $LGN$sampleID$N\n";
}

# outFile
my $outFile = defined $opt_o ? $opt_o : "results.csv";
	$outFile =~ s/\/\.\//\//;
	$outFile =~ s/\/+/\//g;
my $outFilePrint = "<inputFolder>/$outFile" if not defined $opt_o;
$outFile = "$inputFolder/$outFile";
my $total_inputFiles = @inputFiles;

# log parameters
print "
$YW-------------------------
0. Input Parameters
-------------------------$N

-i inputFolder : $LCY$inputFolder$N
-e extension   : $LCY$ext$N
-o outputFile  : $LCY$outFilePrint$N

Total input files: $LGN$total_inputFiles$N

Input files inside ${LGN}inputFolder$N ($LCY$inputFolder$N):
$inputFilesPrint
";


print "
$YW-------------------------
1. Processing .all files
-------------------------$N

";

my %data;
my %type;
my %isotype;
my @def;
for (my $i = 0; $i < @inputFiles; $i++) {

	my ($input1) = $inputFiles[$i];
	my ($folder1, $fileName1) = mitochy::getFilename($input1, "folderfull");

	my ($sampleID) = $fileName1 =~ /([WRSD][123](cut)?)_/;
	die "Cannot get sample from $LCY$fileName1$N!\n" if not defined $sampleID;

	my ($lenMaxL, $lenMaxR) = $fileName1 =~ /^.+.tlx.(\d+)\.(\d+)\.final.tsv.all$/;
	($lenMaxL) = ("UNKNOWN") if not defined ($lenMaxL);
	($lenMaxR) = ("UNKNOWN") if not defined ($lenMaxR);

	my $linecount = 0;
	print "$LGN$i.$N Doing sample=$LGN$sampleID$N lenMaxL=$LGN$lenMaxL$N lenMaxR=$LGN$lenMaxR$N file=$LCY$input1$N\n";
	open (my $in1, "<", $input1) or die "Cannot read from $input1: $!\n";
	my %typecur;
	my %isotypecur;
	while (my $line = <$in1>) {
		chomp($line);
		$linecount ++;
		my @arr = split("\t", $line);

		if ($line =~ /^#/) {
			$arr[0] =~ s/^#//;
			if (@def == 0) {
				@def = @arr;
				for (my $j = 0; $j < @def; $j++) {
					$def[$j] = "type" if $def[$j] eq "type1";
					$def[$j] = "isotype" if $def[$j] eq "igtype";
				}
			}
		}
		else {
			my $totaldef = @def;
			my $totalarr = @arr;
			die "def array ($LGN$totaldef$N) is less or more than row array ($LGN$totalarr$N) at $LCY$input1$N line $LGN$linecount$N! Def and line:\n\n" . join("\t", @def) . "\n$line\n\n" if $totaldef ne $totalarr;
			my $type = $arr[2];
			my $isotype = $arr[1];
			$type{$type} = 1;
			$isotype{$isotype} = 1;
			$typecur{$type} = 1;
			$isotypecur{$isotype} = 1;
			print date() . "    sample=$LGN$sampleID$N,isotype=$LGN$isotype$N,type=$LGN$type$N,\n" . date() . "    " if $linecount == 2;
			for (my $j = 0; $j < @arr; $j++) {
				next if $j eq 2 or $j eq 0 or $j eq 1; #0 is sample, 2 is type, 1 is isotype
				my $def = $def[$j];
				my $val = $arr[$j];
				die "file $input1 line $linecount: Undefined def!\n\n" . join("\t", @def) . "\n$line\n\n" if not defined $def;
				die "file $input1 line $linecount: Undefined val!\n\n" . join("\t", @arr) . "\n$line\n\n" if not defined $val;
				$data{$sampleID}{$isotype}{$type}{$def} = $val;
				print "$def=$LGN$val$N" if $linecount == 2;
				if ($linecount == 2 and $j == 5 and $j != @arr-1) {
					print ",\n" . date() . "    ";
				}
				elsif ($linecount == 2 and $j != @arr-1) {
					print ",";
				}
			}
			print "\n" if $linecount == 2;
		}
	}
	
	my $total_isotypecur = (keys %isotype);
	my $total_typecur    = (keys %type);
	print date() . "  ${LGN}DONE!$N -> Processed $LGN$linecount$N lines, $LGN$total_isotypecur isotypes, $LGN$total_typecur$N types\n\n";
	close $in1;
}

print "
$YW-------------------------
2. Printing into outFile
-------------------------$N

";
#
my $ind = 0;
open (my $out1, ">", "$outFile") or die "Cannot write to $outFile: $!\n";
foreach my $sampleID (sort keys %data) {
	my ($sample, $treat) = get_sampletype($sampleID);
	foreach my $isotype (sort keys %isotype) {
		foreach my $type (sort keys %type) {
			my @defp = qw(sampleID sample treat isotype type);
			my @valp = ($sampleID, $sample, $treat, $isotype, $type);
			my $total_read = $data{$sampleID}{$isotype}{$type}{total_read};
			$total_read = 0 if not defined $total_read;
			my $average_length = $data{$sampleID}{$isotype}{$type}{average_length};
			$average_length = 0 if not defined $average_length;
			foreach my $def (@def[0..@def-1]) {
				next if $def eq "sample" or $def eq "type" or $def eq "isotype";
				$data{$sampleID}{$isotype}{$type}{$def} = 0 if not defined $data{$sampleID}{$isotype}{$type}{$def};
				my $val = $data{$sampleID}{$isotype}{$type}{$def};
				   $val = 0 if not defined $val;
				if ($def eq "count_atleast1") {
					my $perc_atleast1  = $total_read eq 0 ? 0 : 100 * $val / $total_read;
					  ($perc_atleast1) = myformat($perc_atleast1);
					push(@defp, $def); #count_atleast1
					push(@valp, $val); #count_atleast1
					push(@defp, "perc_atleast1");
					push(@valp, $perc_atleast1);
					
				}
				elsif ($def eq "count_total") {
					my $mean_total  = $total_read eq 0 ? 0 : $val / $total_read;
					my $mean_perbp_allread  = $total_read eq 0 ? 0 : $average_length eq 0 ? 0 : ($val / $total_read) / $average_length;
					push(@defp, $def); #count_total
					push(@valp, $val); #count_total
					push(@defp, "mean_total");
					push(@valp, myformat($mean_total));
					push(@defp, "mean_perbp_allread");
					push(@valp, myformat($mean_perbp_allread));
				}
				elsif ($def eq "percperbp") {
					$val = $val / 100; #dont use the perc
					my ($mean_perbp_eachread) = $total_read eq 0 ? 0 : $val / $total_read;
					push(@defp, "total_perbp_eachread"); #totalperbp_eachread
					push(@valp, myformat($val)); #totalperbp_eachread
					push(@defp, "mean_perbp_eachread");
					push(@valp, myformat($mean_perbp_eachread));
				}
				else {
					if ($val =~ /^\-?\d+\.?\d*$/) {
						($val) = myformat($val);
					}
					push(@defp, $def);
					push(@valp, $val);
				}
			}
			if ($ind < 10) {
				print $out1 join(",", @defp) . "\n" if $ind == 0;
				print date() . "Example output:\n\n" if $ind == 0;
				print join("\t", @defp) . "\n" if $ind == 0;
				print join("\t", @valp) . "\n";
			}
			print $out1 join(",", @valp) . "\n";
			$ind ++;
		}
	}
}
close $out1;

print "\n" . date() . "
$YW-------------------------
3. Done!
-------------------------$N
outFile:\n$LCY$outFile$N

cat $outFile | perl -pi -e 's/,/\\t/g' |myless

";

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
