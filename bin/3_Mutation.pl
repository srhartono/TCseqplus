#!/usr/bin/perl
##
use strict; use warnings; use mitochy; use Getopt::Std; use FAlite; use Cwd qw(abs_path); use File::Basename qw(dirname);
use vars qw($opt_v $opt_i $opt_m $opt_n $opt_a $opt_b $opt_o $opt_0 $opt_D $opt_p $opt_I);
getopts("vm:n:a:b:o:i:0D:p:I:");

BEGIN {
   my $libPath = dirname(dirname abs_path $0) . '/lib';
   push(@INC, $libPath);
	#print "\n- Pushed $libPath into perl lib path INC\n";
}

use Mutationlib;

die "
Usage: ${YW}$0${N} [Optional:${LGN} -m <metafile> -n <namewant>${N}] -I $LPR<sampleID>$N -i $CY<W3_result_filter.tlx>${N} 

${YW}IMPORTANT:$N
-I <sampleID> has to be [W/D/R/S][1/2/3].
-i HAS to be [W/D/R/S][1/2/3]_result_filter.tlx if -I <sampleID> is not given

-I will override sampleID from input -i (-i W1_result_filter.tlx -I W3 -> sampleID is W3 and not W1)

${LGN}Optional:${N}

-m <Metafile.TXT> ${YW}[not defined]${N}
	-m 1: ${LCY}/group/stella/Work/Data/Fastq/200323_SLIMS0323_BC_mm_TCseq_JackieBarlowTCseq/0_Fastq/0_META.TXT${N}

-n namewant: e.g.${LCY} M02034:489:000000000-CYYL8:1:1101:16291:1246${N} for W3 ${YW}[not defined]${N}

-o: Output file <default: fileName1.lenMaxL.lenMaxR.final.tsv>
-a: length (bp) of sequence to get left of junction (at IgM) ${YW}[50]${N}
-b: length (bp) of sequence to get right of junction (at far) ${YW}[50]${N}
-0: don't re-run if there's error. Useful to identify it [toggled=TRUE, default=FALSE/can rerun]

-D: DEBUG, process -D number of read [0 = all]

\n\n" unless defined $opt_i and -e $opt_i and defined $opt_m;

my @add = (0);
for (my $ki = 0; $ki < 100; $ki ++) {
	push(@add, -1 * $ki);
	push(@add, $ki);
}

my $inputFile = getFullpath($opt_i);
my ($folder1, $fileName1) = mitochy::getFilename($inputFile, "folderfull");
my ($sampleIDfull) = $fileName1 =~ /^(.*[WSDR][123](cut)?)_result_filter.tlx/;
my ($sampleID) = defined $opt_I ? $opt_I : $inputFile =~ /^.*([WSDR][123])(cut)?_result_filter.tlx/;
my ($sIDprt1, $sIDprt2)  = $inputFile =~ /^.*([A-Z][0-9])(_result_filter.tlx.*)$/;
my $sIDprt3 = defined $opt_I ? $opt_I : "<Not Defined>";
#$sIDprt2 = "<Not Defined>" if not defined $sIDprt2;
my $sIDprt  = defined $sIDprt1 ? $folder1 . "/$LGN" . $sIDprt1 . $N . $sIDprt2 : $inputFile;
my $optIprt = defined $opt_I ? $opt_I : "<Not Defined>";
die "\n${LRD}ERROR!$N:\nUndefined <sampleID>! Example <sampleID>: W1 or S3. Problem:\n   -I <sampleID> undefined\n   -i <inputFile> also not in this format: $LCY<folders>$N/${LGN}\([W/D/R/S][1/2/3]\)${N}_result_filter.tlx). InputFile:\n      $LCY$folder1$N/$LGN$fileName1$N\n \n  Solution: Either supply <sampleID> using -I (e.g. -I W1) or have inputFile in that format!\n\n" if not defined $sampleID;
die "\n${LRD}ERROR!$N:\nsampleID from -I ($LGN$optIprt$N) or -i ($sIDprt) has to be [W/D/R/S][1/2/3]! (current = $LCY$sampleID$N)\n\n" if $sampleID !~ /^[WDRS][123]$/;

# Define other parameters
my $metaFile    = $opt_m; 
	$metaFile    = "/group/stella/Work/Data/Fastq/200323_SLIMS0323_BC_mm_TCseq_JackieBarlowTCseq/0_Fastq/0_META.TXT" if $metaFile eq 1;
	$metaFile    = "/group/stella/Work/Data/Fastq/200903_SLIMS4168_BC_mm_TCseq_JackieBarlowTCseq/0_Fastq/0_META.TXT" if $metaFile eq 2;
	$metaFile    = "/group/stella/Work/Data/Fastq/220303_SLIMS4832_BC_mm_TCseq_JackieBarlowTCseq/0_Fastq/0_META.TXT.FULL" if $metaFile eq 3;
my ($metaFolder, $metaFilename) = getFilename($metaFile, "folderfull");

my $debug_num   = defined $opt_D ? $opt_D : 0;

my $lenMaxL     = defined $opt_a ? $opt_a : 50;

my $lenMaxR     = defined $opt_b ? $opt_b : 50;

my $namewant    = defined $opt_n ? $opt_n : undef;
my $namewantprt = defined $opt_n ? $opt_n : "<Not Defined>";

my $muscleparam = defined $opt_p ? $opt_p : "-maxiters 64";

# define output & log files
# -> logFile is for important only
# -> biglogfile is for all logs
# -> outbedfile is for bed file coordinate
my $outFile     = defined $opt_o ? $opt_o : "$folder1/$fileName1.$lenMaxL.$lenMaxR.final.tsv";
my $logFile     = $outFile . ".LOG";
my $biglogFile  = $outFile . ".BIGLOG";
my $outzfixFile = $outFile . ".outzfix.bed.temp";
my $outztempFile = $outFile . ".outztemp.bed.temp";
my $outzfaFile = $outFile . ".outztemp.fa";
open (my $out1,      ">", "$outFile") or die "Failed to write to ${LCY}$outFile${N}: $!\n\n";
open (my $outBigLog, ">", "$biglogFile") or die "Failed to write to ${LCY}$biglogFile${N}: $!\n\n";
open (my $outLog,    ">", "$logFile") or DIELOG($outBigLog, "Failed to write to ${LCY}$logFile${N}: $!\n\n");
open (my $outloc,    ">", "$outFile.loc") or DIELOG($outBigLog, "Failed to write to $LCY$outFile.loc$N: $!\n\n");
# Create tempFolder $folder1/.temp
my $tempFolder = "$folder1/.temp";
if (not -d "$tempFolder") {
	system("mkdir -p $tempFolder");
	LOG($outBigLog, " -> Created .temp. folder: $LCY$tempFolder$N temp\n\n");
}
else {
	LOG($outBigLog, " -> Found .temp folder: $LCY$tempFolder$N\n\n");
}


my $runScript = "$0 -i $inputFile";
$runScript .= " -I $opt_I" if defined $opt_I;
$runScript .= " -m $metaFile" if defined $metaFile;
$runScript .= " -a $opt_a" if defined $opt_a;
$runScript .= " -b $opt_b" if defined $opt_b;
$runScript .= " -o $opt_o" if defined $opt_o;
$runScript .= " -n $namewant" if defined $opt_n;
$runScript .= " -0" if defined $opt_0;
my $canrerun = defined $opt_0 ? "FALSE" : "TRUE";

my $linecount;
my $dashhead = ${LGN} . join("", ("-") x 50) . ${N};
my $DATA;
my $total_read = 0;


# -------------------------------------------
# 0. Print Parameters
# -------------------------------------------
my $param_print = "
-i InputFile    = ${LCY}$inputFile${N}
-I expl.samID   = ${LCY}$sIDprt3${N}
   SampleID     = ${LCY}$sampleID${N}
   SampleIDFull = ${LCY}$sampleIDfull${N}
-m metaFile     = ${LCY}$metaFile${N}
-a lenMaxL      = ${LCY}$lenMaxL${N}
-b lenMaxR      = ${LCY}$lenMaxR${N}
-n namewant     = ${LCY}$namewantprt${N}
-o outFile      = ${LCY}$outFile${N}
   logFile      = ${LCY}$logFile${N}
-0 can rerun    = ${LCY}$canrerun${N}

runScript=$runScript\n\n";

# header_print
my $header_print_0_param = header_print("0. PARAMETERS:", "", $param_print);
LOG($outBigLog, $header_print_0_param);

# Parse consensus fix file

my %cons = %{get_cons_fix_hash()};

sub get_cons_fix_hash {
	my %cons;
	my $IgMFile = "/home/mitochi/Work/Project/TCseqplus/JackieTCseq/FA/S129_IgM_fa.aln.fa.out";
	my $IgG1File = $IgMFile; $IgG1File =~ s/_IgM_/_IgG1_/;
	my $IgG3File = $IgMFile; $IgG3File =~ s/_IgM_/_IgG3_/;
	die "Cannot find IgMFile $LCY$IgMFile$N!\n" if not -e $IgMFile;
	die "Cannot find IgG1File $LCY$IgG1File$N!\n" if not -e $IgG1File;
	die "Cannot find IgG3File $LCY$IgG3File$N!\n" if not -e $IgG3File;

	$cons{IgM} = get_cons_fix_main($IgMFile);
	$cons{IgG1} = get_cons_fix_main($IgG1File);
	$cons{IgG3} = get_cons_fix_main($IgG3File);

	return(\%cons);
}

sub get_cons_fix_main {
	my ($file) = @_;
	my %cons;
	my $linecount = 0;
	open (my $in, "<", $file) or die "Failed to read from $file: $!\n";
	while (my $line = <$in>) {
		chomp($line);
		$linecount ++;
		my @arr = split("\t", $line);
		next if @arr ne 7;
		my ($pos, $beg1, $end1, $beg2, $end2, $type, $count) = @arr;
		$count =~ s/^count\=//;
		my ($type1, $type2) = $type =~ /^([A-Za-z0-9]+)_(.+)$/;
		$type2 = "$type2\_$type2" if $type1 eq "mat";
		if ($type1 eq "ins") {
			LOG($outBigLog, "file=$LCY$file$N, insertion at pos=$pos, beg1=$beg1, beg2=$beg2, count=$count\n");
			my @type2 = split("", $type2);
			for (my $i = 0; $i < @type2; $i++) {
				my $ind = $beg1 + $i;
				die "Already defined cons 2 beg1=$beg1 i=$i ind=$ind type1=$cons{2}{$beg1}{$ind}{type1}, type2 = $type2[$i]\n" if defined $cons{2}{$beg1}{$ind}{type1};
				$cons{2}{$beg1}{$ind}{pos} = $pos;
				$cons{2}{$beg1}{$ind}{beg2} = $ind;
				$cons{2}{$beg1}{$ind}{type1} = $type1;
				$cons{2}{$beg1}{$ind}{type2} = $type2[$i];
				$cons{2}{$beg1}{$ind}{count} = 1;
			}
		}
		if ($type1 ne "del") {
			$cons{1}{$beg1}{pos} = $pos;
			$cons{1}{$beg1}{beg2} = $beg2;
			$cons{1}{$beg1}{type1} = $type1;
			$cons{1}{$beg1}{type2} = $type2;
			$cons{1}{$beg1}{count} = $count;
		}
		else {
			LOG($outBigLog, "file=$LCY$file$N, deletion at pos=$pos, beg1=$beg1, beg2=$beg2, count=$count\n");
			my @type2 = split("", $type2);
			for (my $i = 0; $i < @type2; $i++) {
				my $ind = $beg1 + $i;
				die "Already defined cons 1 i=$ind type1=$cons{1}{$ind}{type1}, type2 = $type2[$i]\n" if defined $cons{1}{$ind}{type1};
				$cons{1}{$ind}{pos} = $pos;
				$cons{1}{$ind}{beg2} = $ind;
				$cons{1}{$ind}{type1} = $type1;
				$cons{1}{$ind}{type2} = $type2[$i];
				$cons{1}{$ind}{count} = 1;
			}
		}
		
	}
	close $in;
	return(\%cons);
}
# -------------------------------------------
# 1. Processing $metaFile <metaFile.txt>
# -------------------------------------------

# header print
my $metaFile_header_print = header_print("1. PROCESSING META file:", "   ${LCY}$metaFile${N}\n   TO GET BARCCODE+PRIMER & ADAPTER SEQUENCS");
LOG($outBigLog, $metaFile_header_print);

# parse
my ($seqPRI, $seqADP, $metaFile_log, $metaFile_die) = parse_metaFile($metaFile, $sampleID);
$metaFile_die ne 0 ? DIELOG($outBigLog, $metaFile_log) : LOG($outBigLog, $metaFile_log);

# check total_read
$total_read = scalar(keys %{$DATA});
LOG($outBigLog, "\nTotal Read = $LGN$total_read$N\n");

# -------------------------------------------
# 2. Processing $inputFile <filter.tlx>
# -------------------------------------------

my $outBedFile  = "$tempFolder/$fileName1.$lenMaxL.$lenMaxR.final.bed";

# header print
my $inputFile_header_print = header_print("2. PROCESSING INPUT <filter.tlx file>:","   ${LCY}$inputFile${N}\n   TO GET REF. BED COORDINATES\n   $LGN$outBedFile$N\n");
LOG($outBigLog, $inputFile_header_print);

# parse
my ($inputFile_log, $inputFile_die);
($DATA, $namewant, $inputFile_log, $inputFile_die) = parse_inputFile($inputFile, $outBedFile, $namewant);
LOG($outBigLog, $inputFile_log);
if ($inputFile_die ne 0) {
	exit(1);
}

# check total_read
$total_read = scalar(keys %{$DATA});
LOG($outBigLog, "\nTotal Read = $LGN$total_read$N\n");

# -------------------------------------------
# 3. Parsing reference seq fasta file
# -------------------------------------------

my $inputFA = "$outBedFile.fa";

# header print
my $outBedFile_header_print = header_print("3. DOING <fastaFromBed> ON REF BED FILE:","   ${LCY}$outBedFile${N}\n   TO GET REF. SEQ FASTA FILE AND PARSE IT:\n   $LGN$outBedFile.fa$N\n");
LOG($outBigLog, $outBedFile_header_print);

# parse
my ($outBedFile_log, $outBedFile_die);
($DATA, $outBedFile_log, $outBedFile_die) = parse_outBedFile($DATA, $outBedFile, $inputFA, $namewant);
LOG($outBigLog, $outBedFile_log);
if ($outBedFile_die ne 0) {
	exit(1);
}

# check total_read
$total_read = scalar(keys %{$DATA});
LOG($outBigLog, "\nTotal Read = $LGN$total_read$N\n");

# -------------------------------------------
# 4. Getting QUERY seq file from fastq
# -------------------------------------------

#my @inputFQ_R1 = <$folder1/../../0_Fastq/*$sampleIDfull\_R1.fq.gz>;
my @inputFQ_R1 = <$metaFolder/*$sampleIDfull\_R1.fq.gz>;
DIELOG($outBigLog, "Can't find input fastq R1 ($metaFolder/*$sampleIDfull\_R1.fq.gz)\n") if @inputFQ_R1 == 0;
DIELOG($outBigLog, "\nFound multiple input fastq R1 ($metaFolder/*$sampleIDfull\_R1.fq.gz)\n") if @inputFQ_R1 > 1;

#my @inputFQ_R2 = <$folder1/../../0_Fastq/*$sampleIDfull\_R2.fq.gz>;
my @inputFQ_R2 = <$metaFolder/*$sampleIDfull\_R2.fq.gz>;
DIELOG($outBigLog, "Can't find input fastq R1 ($metaFolder/*$sampleIDfull\_R1.fq.gz)\n") if @inputFQ_R2 == 0;
DIELOG($outBigLog, "Found multiple input fastq R2 ($metaFolder/*$sampleIDfull\_R2.fq.gz)\n") if @inputFQ_R2 > 1;


# header print
my $outFastqFile_header_print = header_print("4. PARSING QUERY FASTQ FILE:","   FQ1 = ${LCY}$inputFQ_R1[0]${N}\n   FQ2 = ${LCY}$inputFQ_R2[0]${N}\n   TO GET QUERY SEQUENCE\n   total read parsed = $LGN$total_read$N\n");
LOG($outBigLog, $outFastqFile_header_print);

my ($fastqFile_log, $fastqFile_die) = ("", 0);

# parse read1
($DATA, $fastqFile_log, $fastqFile_die) = parse_fastqFile($DATA, $inputFQ_R1[0], "R1", $total_read, $namewant);
LOG($outBigLog, $fastqFile_log);
if ($fastqFile_die ne 0) {
	exit(1);
}

# parse read2
($DATA, $fastqFile_log, $fastqFile_die) = parse_fastqFile($DATA, $inputFQ_R2[0], "R2", $total_read, $namewant);
LOG($outBigLog, $fastqFile_log);
if ($fastqFile_die ne 0) {
	exit(1);
}

# check total_read
$total_read = scalar(keys %{$DATA});
LOG($outBigLog, "\nTotal Read = $LGN$total_read$N\n");

# -------------------------------------------
# 5. Aligning and printing into outFile final.tsv
# -------------------------------------------

my $outoutFile_header_print = header_print("5. ALIGNING & PRINTING INTO outFile final.tsv:","   outFile = ${LCY}$outFile${N}\n");
LOG($outBigLog, $outoutFile_header_print);


my $totalall = 0;
my %totalpair;
my $good = 0;
my %okay;
my $nexted_sequndef = 0;
my $nexted_others = 0;
my $nexted_bothshort = 0;
my $nexted_begshort = 0;
my $nexted_endshort = 0;
my $nexted_juncmissing = 0;
my $print = 0;
#my $header = "name\ttype\tnuc1\tnuc2\tnumber\tperc\tpercnuc1\tpercnuc2\tmean\tsd\tchr1\tbeg1\tend1\tjunc1\tlength1\tstrand1\tchr2\tbeg2\tend2\tjunc2\tlength2\tstrand2\n";
my $header = "name\treadorient1\ttype\tnuc1\tnuc2\tnumber\tperc\ttotal_read\tigtype\tmean\tsd\tchr1\tbeg1\tend1\tjunc1\tlen1\tstrand1\tchr2\tbeg2\tend2\tjunc2\tlen2\tstrand2\tmutpos_joined\n";
print $out1 "$header";
my ($igtypez, $typez1, $typez2, $all);
foreach my $name (sort {$DATA->{$a}{order} <=> $DATA->{$b}{order}} keys %{$DATA}) {
	my $order = $DATA->{$name}{order};
	print "$order. Doing $name\n" if $order % 100 == 0;
#	my ($namez, $juncid, $chr2, $junk, $strandjunk2, $beg2, $end2, $chr1, $beg1, $end1, $strandjunk1, $begq, $endq, $lenq) = split("\t", $DATA->{$name}{coor});
	$DATA->{$name} = get_seqTR($DATA->{$name}, $name, $outBigLog, $namewant);
	my $chr1 = $DATA->{$name}{chr1};
	my $chr2 = $DATA->{$name}{chr2};
	my $beg1 = $DATA->{$name}{beg1};
	my $beg2 = $DATA->{$name}{beg2};
	my $end1 = $DATA->{$name}{end1};
	my $end2 = $DATA->{$name}{end2};
	my $strand1 = $DATA->{$name}{strand1};
	my $strand2 = $DATA->{$name}{strand2};
	my $begJ1 = $DATA->{$name}{begJ1};
	my $begJ2 = $DATA->{$name}{begJ2};
	my $seqC1L = $DATA->{$name}{seqC1L};
	my $seqC1R = $DATA->{$name}{seqC1R};
	my $seqC2R = $DATA->{$name}{seqC2R};
	my $seqC2L = $DATA->{$name}{seqC2L};
	my $seqC1 = $seqC1L . $seqC1R;
	my $seqC2 = $seqC2L . $seqC2R;
	my $lenz1 = $DATA->{$name}{lenC1L};
	my $lenz2 = $DATA->{$name}{lenC2R};
	my $seqQ1 = $DATA->{$name}{seqQ1};
	my $seqQ2 = $DATA->{$name}{seqQ2};
	my $seqCON = $DATA->{$name}{seqCON};
	my $cigarQ1 = $DATA->{$name}{cigarQ1};	
	my $cigarQ2 = $DATA->{$name}{cigarQ2};	
	my $begT1 = $DATA->{$name}{begT1};	
	my $endT1 = $DATA->{$name}{endT1};	
	my $begT2 = $DATA->{$name}{begT2};	
	my $endT2 = $DATA->{$name}{endT2};	
	my $seqTN = $DATA->{$name}{seqTN};	
	my $begTN = $DATA->{$name}{begTN};
	my $endTN = $DATA->{$name}{endTN};
	my $lenTN = $DATA->{$name}{lenTN};
	my $mhbeg = $DATA->{$name}{mhbeg};
	my $mhend = $DATA->{$name}{mhend};
	my $mhpos = $DATA->{$name}{mhpos};
	my $mhseq = $DATA->{$name}{mhseq};
	my @mhseq = split("", $mhseq);
	my $seqTR = $DATA->{$name}{seqTR};


	if ($DATA->{$name}{gap} > 0) {
		my ($seqTNL, $seqTNG, $seqTNR) = $seqTN =~ /^(.{$DATA->{$name}{endT1}})(.{$DATA->{$name}{gap}})(.*)$/;
#		my ($bcjunk, $seqTNL, $seqTNG, $seqTNR) = $seqTN =~ /^(.{$DATA->{$name}{begT1}})(.{$DATA->{$name}{endT1}})(.{$DATA->{$name}{gap}})(.*)$/;
		my $lenQ2get = $DATA->{$name}{lenTN} - $DATA->{$name}{begT2};
#		my $lenQ2get = $DATA->{$name}{endT2} - $DATA->{$name}{begT2};
#		my ($bc, $seqQ1L, $seqQ1G, $seqQ1R) = $seqQ1 =~ /^(.{$DATA->{$name}{begT1}})(.{$DATA->{$name}{endT1}})(.{$DATA->{$name}{gap}})(.*)$/;
		my ($seqQ1L, $seqQ1G, $seqQ1R, $seqQ1temp);
		($seqQ1L, $seqQ1temp) = $seqQ1 =~ /^(.{$DATA->{$name}{endT1}})(.*)$/;
		if (not defined $seqQ1temp) {
			$seqQ1L = $seqQ1;
			$seqQ1G = "";
			$seqQ1R = "";
		}
		else {
			($seqQ1L, $seqQ1G, $seqQ1temp) = $seqQ1 =~ /^(.{$DATA->{$name}{endT1}})(.{$DATA->{$name}{gap}})(.*)$/;
			if (not defined $seqQ1temp) {
				($seqQ1L, $seqQ1G) = $seqQ1 =~ /^(.{$DATA->{$name}{endT1}})(.*)$/;
				$seqQ1R = "";
			}
			else {
				($seqQ1L, $seqQ1G, $seqQ1R) = $seqQ1 =~ /^(.{$DATA->{$name}{endT1}})(.{$DATA->{$name}{gap}})(.*)$/;
				$seqQ1R = "" if not defined $seqQ1R;
			}
		}

		my ($seqQ2L, $seqQ2G, $seqQ2R, $seqQ2temp);
		if ($lenQ2get < 0) {
			$seqQ2L = "";
			$seqQ2G = "";
			$seqQ2R = $seqQ2;
		}
		else {
			($seqQ2temp, $seqQ2R) = $seqQ2 =~ /^(.*)(.{$lenQ2get})$/;
			if (not defined $seqQ2temp) {
				$seqQ2L = "";
				$seqQ2G = "";
				$seqQ2R = $seqQ2;
			}
			else {

				($seqQ2temp, $seqQ2G, $seqQ2R) = $seqQ2 =~ /^(.*)(.{$DATA->{$name}{gap}})(.{$lenQ2get})$/;
				if (not defined $seqQ2temp) {
					$seqQ2L = "";
					($seqQ2G, $seqQ2R) = $seqQ2 =~ /^(.*)(.{$lenQ2get})$/;
				}
				else {
					($seqQ2L, $seqQ2G, $seqQ2R) = $seqQ2 =~ /^(.*)(.{$DATA->{$name}{gap}})(.{$lenQ2get})$/;
				}
			}
		}
		$seqQ1L = "" if not defined $seqQ1L;
		$seqQ1G = "" if not defined $seqQ1G;
		$seqQ1R = "" if not defined $seqQ1R;
		$seqQ2L = "" if not defined $seqQ2L;
		$seqQ2G = "" if not defined $seqQ2G;
		$seqQ2R = "" if not defined $seqQ2R;
		my ($get1) = $DATA->{$name}{endT1};
		my ($getgap) = $DATA->{$name}{gap};
		LOG($outBigLog, "begT1=$DATA->{$name}{begT1}-$DATA->{$name}{endT1}, begT2=$DATA->{$name}{begT2}-$DATA->{$name}{endT2}, lengthTN=$DATA->{$name}{lenTN}, gap=$DATA->{$name}{gap}, lenQ2get=$lenQ2get\n") if $name eq $namewant;
		LOG($outBigLog, "lengthQ1=" . length($seqQ1) . ": get endT1=$DATA->{$name}{endT1} + gap=$DATA->{$name}{gap} + (rest)\n") if $name eq $namewant;
		LOG($outBigLog, "lengthQ2=" . length($seqQ2) . ": get (rest) +  gap=$DATA->{$name}{gap} + lenQ2get=$lenQ2get (lenTN=$DATA->{$name}{lenTN} - begT2=$DATA->{$name}{begT2})\n") if $name eq $namewant;

		LOG($outBigLog, "seqQ1=$seqQ1\n seqQ1L=$seqQ1L\n seqQ1G=$seqQ1G\n seqQ1R=$seqQ1R\n\n");
		LOG($outBigLog, "seqQ2=$seqQ2\n seqQ2L=$seqQ2L\n seqQ2G=$seqQ2G\n seqQ2R=$seqQ2R\n\n");
		$seqQ1 = $seqQ1L . $seqQ1R;
		$seqQ2 = $seqQ2L . $seqQ2R;
		$seqTN = $seqTNL . $seqTNR;
		$DATA->{$name}{seqQ1} = $seqQ1;
		$DATA->{$name}{seqQ2} = $seqQ2;
		$DATA->{$name}{seqTN} = $seqTN;
		$DATA->{$name}{seqQ1G} = $seqQ1G;
		$DATA->{$name}{seqQ2G} = $seqQ2G;
		$DATA->{$name}{seqTNG} = $seqTNG;
		my @res = muscle(
">seqQ1_aln\n$seqQ1
>seqQ2_aln\n$seqQ2
>seqTN_aln\n$seqTN
");

		my $res2 = ">seqTL_aln\n$seqTNL\n";
		if (length($seqQ1L) > 0) {
			$res2 .= ">seq1L_aln\n$seqQ1L\n";
		}
		if (length($seqQ2L) > 0) {
			$res2 .= ">seq2L_aln\n$seqQ2L\n";
		}
		my @res2 = muscle($res2);

		my $res3 = ">seqTR_aln\n$seqTNR\n";
		if (length($seqQ1R) > 0) {
			$res3 .= ">seq1R_aln\n$seqQ1R\n";
		}
		if (length($seqQ2R) > 0) {
			$res3 .= ">seq2R_aln\n$seqQ2R\n";
		}
		my @res3 = muscle($res3);

		my $res = parse_fasta_simple(undef, \@res);
		$res2 = parse_fasta_simple(undef, \@res2);
		if (length($seqQ1L) == 0) {
			$res2->{seq1L_aln} = join("", ("-") x length($res2->{seqTL_aln}));
		}
		if (length($seqQ2L) == 0) {
			$res2->{seq2L_aln} = join("", ("-") x length($res2->{seqTL_aln}));
		}
		$res3 = parse_fasta_simple(undef, \@res3);
		if (length($seqQ1R) == 0) {
			$res3->{seq1R_aln} = join("", ("-") x length($res3->{seqTR_aln}));
		}
		if (length($seqQ2R) == 0) {
			$res3->{seq2R_aln} = join("", ("-") x length($res3->{seqTR_aln}));
		}		
		LOG($outBigLog, "\nHERE!!! endT1 = $DATA->{$name}{endT1}\n");
		LOG($outBigLog, "seqQ1\t" . colorize($res->{seqQ1_aln}) . "\n");
		LOG($outBigLog, "seqTN\t" . colorize($res->{seqTN_aln}) . "\n");
		LOG($outBigLog, "seqQ2\t" . colorize($res->{seqQ2_aln}) . "\n\n");
#		LOG($outBigLog, "seqTL\t" . colorize($res->{seqTL_aln}) . "\n");
#		LOG($outBigLog, "seqTG\t" . colorize($res->{seqTG_aln}) . "\n");
#		LOG($outBigLog, "seqTR\t" . colorize($res->{seqTR_aln}) . "\n\n");

		LOG($outBigLog, "seqQ1\t" . colorize($res2->{seq1L_aln}) . $DATA->{$name}{seqQ1G} . colorize($res3->{seq1R_aln}) . "\n");
		LOG($outBigLog, "seqTN\t" . colorize($res2->{seqTL_aln}) . $DATA->{$name}{seqTNG} . colorize($res3->{seqTR_aln}) . "\n");
		LOG($outBigLog, "seqQ2\t" . colorize($res2->{seq2L_aln}) . $DATA->{$name}{seqQ2G} . colorize($res3->{seq2R_aln}) . "\n");
		LOG($outBigLog, "\nHERE!!!\n");


	}

	$totalall ++;
	$totalpair{$name} = 1;
#
	my $seqTNprint = $seqTN;
	if (not defined $seqCON) {
		$nexted_sequndef ++;
		LOG($outBigLog, "$LRD ERROR!!!${N} $LPR name=$name${N} is skipped as seqCON is not defined!!! Data parameters and value for name=$name:\n");
		foreach my $param (sort keys %{$DATA->{$name}}) {
			my $value = $DATA->{$name}{$param}; $value = "UNDEFINED" if not defined $value;
			LOG($outBigLog, "$param=$value\n");
		}

		if (not defined $opt_0) {
			DIELOG($outBigLog, date() . "\nPlease run this runscript:\n" . $runScript . " -n $name -0\n\n");
		}

		next;
	}
	my $strand1print = $strand1 eq "+" ? "pos" : "neg";
	my $strand2print = $strand2 eq "+" ? "pos" : "neg";
	

	my @cmdTN = muscle(
">seqCON_aln\n$seqCON
>seqQ1_aln\n$seqQ1
>seqQ2_aln\n$seqQ2
>seqC1_aln\n$seqC1L$seqC1R
>seqC2_aln\n$seqC2L$seqC2R
>seqTN_aln\n$seqTN
>seqTR_aln\n$seqTR
>seqJN_aln\n$seqC1L
>seqPR_aln\n$seqPRI
>seqAD_aln\n$seqADP
");
	
	@cmdTN = muscle(
">seqCON_aln\n$seqCON
>seqQ1_aln\n$seqQ1
>seqQ2_aln\n$seqQ2
>seqC1_aln\n$seqC1L$seqC1R
>seqC2_aln\n$seqC2L$seqC2R
>seqTN_aln\n$seqTN
>seqTR_aln\n$seqTR
");
	($DATA->{$name}) = parse_fasta_simple($DATA->{$name}, \@cmdTN);

	my $orig_TN_aln = $DATA->{$name}{'seqTN_aln'};


#	my @resTN = `$cmdTN`;

#	($DATA->{$name}) = parse_fasta_simple($DATA->{$name}, \@resTN);

	LOG($outBigLog, $DATA->{$name}{line});
#	LOG($outBigLog, "mhbeg after mod from base 1/1 to 0/1 = $DATA->{$name}{mhbeg}-$DATA->{$name}{mhend}\n");
	$mhbeg = $DATA->{$name}{mhbeg};
	$mhend = $DATA->{$name}{mhend};

	LOG($outBigLog, "\n" . date() . "\n" . "EVERYTHING BELOW IS NOW 0/0 BASED!\n" . 
"\n>seqC1L (REFERENCE SEQ, LEFT of JUNC)\n......" . colorize($seqC1L) . 
"\n>seqTLX (QUERY READ SEQUENCE)\n" . colorize($seqTN) . 
"\n>seqC2R (REFERENCE SEQ, RIGHT of JUNC)\n      " . join("", (" ") x length($seqC1L)) . colorize($seqC2R) . "\n\n") if $name eq $namewant;

	LOG($outBigLog, date() . "\n${LPR}------------\n${N}") if $name eq $namewant;

	# 1. JUNC
	LOG($outBigLog, date() . "${LGN}Find JUNC position based on READ sequence (seqTN) and LEFT REF SEQ (seqC1L) ${N}\n") if $name eq $namewant;

#	my ($junc_fix_pos) = fix_pos($DATA->{$name}{lenC1L} + 6 - 1, $seqTN, $DATA->{$name}{'seqTN_aln'}); #-1 to make from 1 based to 0 based
#	my ($newjunc) = fix_pos($junc, $DATA->{$name}{'seqTN_aln'}, $seqwants);
#	LOG($outBigLog, "\n\n prev_junc=$DATA->{$name}{lenC1L}, junc_fix_pos = $LGN$junc_fix_pos$N\n\n");


	my $testjunc1 = ""; my $testjunc2 = "";
	if (length($seqTN) > 100) {
		($testjunc1, $testjunc2) = $seqTN =~ /^(.{100})(.+)$/;
		LOG($outBigLog, "SEQTN use 100! ($testjunc1)\n");
	}
	elsif (length($seqTN) > 200) {
		($testjunc1, $testjunc2) = $seqTN =~ /^(.{200})(.+)$/;
		LOG($outBigLog, "SEQTN use 200! ($testjunc1)\n");
	}
	elsif (length($seqTN) > 300) {
		($testjunc1, $testjunc2) = $seqTN =~ /^(.{300})(.+)$/;
		LOG($outBigLog, "SEQTN use 300! ($testjunc1)\n");
	}
	else {
		($testjunc1) = $seqTN =~ /^(.+)$/;
		$testjunc2 = "";
		LOG($outBigLog, "SEQTN use ALL! ($testjunc1)\n");
	}
	my $lengthseqC1L = length($seqC1L) + 6;
	LOG($outBigLog, "testjunc1=$testjunc1, length=" . length($testjunc1) . ", lengthseqC1l = " . $lengthseqC1L . "\n");
	if (length($testjunc1) < $lengthseqC1L) {
		LOG($outBigLog, "testjunc1 < length seqC1L\n");
		my $lengthget = $lengthseqC1L;
		$lengthget = length($seqTN) > $lengthget + 20 ? $lengthget + 20 : 
						 length($seqTN) > $lengthget + 10 ? $lengthget + 10 :
						 length($seqTN) > $lengthget +  5 ? $lengthget +  5 :
						 length($seqTN) > $lengthget +  3 ? $lengthget +  3 :
						 length($seqTN) > $lengthget +  1 ? $lengthget +  1 :
						 length($seqTN) >= $lengthget +  0 ? $lengthget +  0 :
						 length($seqTN) < $lengthget +  0 ? length($seqTN) : $lengthget;

		($testjunc1, $testjunc2) = $seqTN =~ /^(.{$lengthget})(.*)$/;
		($testjunc1) = $seqTN =~ /^(.{$lengthget})/ if not defined $testjunc1 or (defined $testjunc1 and $testjunc1 eq "");
		$testjunc2 = "" if not defined $testjunc2;
	}
	else {
		LOG($outBigLog, "testjunc1=" . length($testjunc1) . " > length seqC1L=" . $lengthseqC1L . "\n");
	}

#	@cmdTN = muscle(
#">seqTN_aln\n$seqTN
#>seqJN_aln\n$seqC1L
#");

	@cmdTN = muscle(
">seqTN_aln\n$testjunc1
>seqJN_aln\n$seqC1L
");

	($DATA->{$name}) = parse_fasta_simple($DATA->{$name}, \@cmdTN);
	my $testjuncdash = join("", ("-") x length($testjunc2));
	my $resjun1 = ">1_neg_$name\n$DATA->{$name}{'seqTN_aln'}$testjunc2\n>0_JUN_$name\n$DATA->{$name}{'seqJN_aln'}$testjuncdash";

	LOG($outBigLog, "TESTJUNC:\n$resjun1\n");
	my @resjun1 = split("\n", $resjun1);
	my ($junc, $begdash, $dashes) = parse_aln(\@resjun1, "1_", "junc", undef, 1, $lenMaxL, $lenMaxR, $name, $DATA->{$name}{gap}, $outBigLog, $outLog, $namewant);

	LOG($outBigLog, "testjunc1=" . colorize($testjunc1) . "\nJUNC=$junc\n");

	if ($junc ne -1) {
		LOG($outBigLog, "Junc is good! So use all testjunc1 ($testjunc1)\n");
	#	$DATA->{$name}{'seqTN_aln'} .= $testjunc2;
	}
	else {
		LOG($outBigLog, "Junc is -1 so use all seq\n");
		@cmdTN = muscle(
">seqTN_aln\n$seqTN
>seqJN_aln\n$seqC1L
");

		($DATA->{$name}) = parse_fasta_simple($DATA->{$name}, \@cmdTN);

		$resjun1 = ">1_neg_$name\n$DATA->{$name}{'seqTN_aln'}\n>0_JUN_$name\n$DATA->{$name}{'seqJN_aln'}";
		LOG($outBigLog, "USEALL:\n$resjun1\n");
		@resjun1 = split("\n", $resjun1);
		($junc, $begdash, $dashes) = parse_aln(\@resjun1, "1_", "junc", undef, 1, $lenMaxL, $lenMaxR, $name, $DATA->{$name}{gap}, $outBigLog, $outLog, $namewant);
	}


	my $orig_junc = $junc;
	$junc -= 1;
	my $orig_junc_min1 = $junc;
	my $junc_TN_aln = $DATA->{$name}{'seqTN_aln'};
	($junc) = fix_pos($junc, $junc_TN_aln, $orig_TN_aln);
	my $ref_junc = abs($begT1 - $endT1);
	LOG($outBigLog, "\t${LGN}FINAL: $junc${N}. (QUERY junc = $orig_junc - 1 = $orig_junc_min1, fix_pos to orig TN = ${LGN}$junc${N} (used)\n");# if $name eq $namewant;# ;beg1 minus end1 abs($begT1-$endT1) = $ref_junc + $begdash + $dashes)\n\n") 
#	LOG($outBigLog, "\tFIX_POS_JUNC DIFFERENT\n") if $junc ne $junc_fix_pos;

	die "Undef junc!\n" if not defined $junc;

	$DATA->{$name}{junc} = $junc;
#	$junc = $DATA->{$name}{junc};
	my @seqTNs = split("", $orig_TN_aln);#DATA->{$name}{'seqTN_aln'});
	if ($mhpos ne -1 and defined $junc and $junc > 0) {
		$DATA->{$name}{mhend} = $junc;			
		for (my $z = $junc - ($DATA->{$name}{mhpos}-1); $z >= 0; $z--) {
			$DATA->{$name}{mhbeg} = $z if $seqTNs[$z] ne "-";
			last if $seqTNs[$z] ne "-";
		}
		LOG($outBigLog, " ->$LGN FIXED MHBEG/END from 0/1 to 0/0!$N (MHPOS = $LGN$DATA->{$name}{mhpos}$N, junc = $LGN$DATA->{$name}{junc}$N\n");
		LOG($outBigLog, "   - MHPOS = ${LGN}$DATA->{$name}{mhpos}${N}\n");
		LOG($outBigLog, "   - MHBEG-END Before : $LGN$mhbeg$N-$LCY$mhend$N\n");
		LOG($outBigLog, "   - MHBEG-END After  : MHBEG-MHEND = $LGN$DATA->{$name}{mhbeg}$N-$LCY$DATA->{$name}{mhend}$N\n\n");
	}


	# 2. PRIMER
	@cmdTN = muscle(
">seqTN_aln\n$seqTN
>seqPR_aln\n$seqPRI
");
	($DATA->{$name}) = parse_fasta_simple($DATA->{$name}, \@cmdTN);
	LOG($outBigLog, date() . "${LGN}Finding location of primer ${N}\[" . colorize($seqPRI) . "${N}]${LGN} and beg_pos of READ sequence ($seqTNprint):${N}\n") if $name eq $namewant;
	my $respri1 = ">1_neg_$name\n$DATA->{$name}{'seqTN_aln'}\n>0_PRI_$name\n$DATA->{$name}{'seqPR_aln'}";
	my @respri1 = split("\n", $respri1);
	my ($beg_pos, $end_pos, $primer_pos, $adapter_pos, $primer_pos_fix, $adapter_pos_fix) = (-1,-1,-1,-1,-1,-1);

	my $primer_TN_aln = $DATA->{$name}{'seqTN_aln'};
	my ($newjunc) = fix_pos($junc, $orig_TN_aln, $primer_TN_aln);
	($primer_pos, $primer_pos_fix, $beg_pos, $begdash, $dashes) = parse_aln(\@respri1, "1_", "primer,-1,-1,$newjunc", undef, 1, $lenMaxL, $lenMaxR, $name, $DATA->{$name}{gap}, $outBigLog, $outLog, $namewant);
	LOG($outBigLog, "primer begdash=$begdash, dashes=$dashes\n"); #if $name eq $namewant;

	($primer_pos) = fix_pos($primer_pos, $primer_TN_aln, $orig_TN_aln);
	($primer_pos_fix) = fix_pos($primer_pos_fix, $primer_TN_aln, $orig_TN_aln);
	($beg_pos) = fix_pos($beg_pos, $primer_TN_aln, $orig_TN_aln);

	LOG($outBigLog, date() . "${LGN}Finding location of adapter ${N}\[" . colorize($seqADP) . "${N}]${LGN} and end_pos of READ sequence ($seqTNprint):${N}\n"); #if $name eq $namewant;

#	my $cmdadp1 = "echo '>1_$strand1print\_$name\n$seqTN\n>0_ADP\_$name\n$seqADP\n'";
   my ($seqADPR) = $orig_TN_aln =~ /$seqADP(.*)$/; $seqADPR = "" if not defined $seqADPR;
   my ($seqADPL) = $orig_TN_aln =~ /^(.*)$seqADP/; $seqADPL = "" if not defined $seqADPL;
	my $resadp1;
	my @resadp1;
	my $adapter_TN_aln;
	if ($seqADPL ne "" and $seqADPR ne "") {
		$adapter_pos = length($orig_TN_aln) - length($seqADPR) - length($seqADP); 
		$DATA->{$name}{seqAD_aln} = join("", ("-") x $adapter_pos) . $seqADP . join("", ("-") x (length($orig_TN_aln) - $adapter_pos - length($seqADP)));
		$resadp1 = ">1_neg_$name\n$orig_TN_aln\n>0_ADP_$name\n$DATA->{$name}{'seqAD_aln'}";
		@resadp1 = split("\n", $resadp1);
		($adapter_pos, $adapter_pos_fix, $end_pos, $begdash, $dashes) = parse_aln(\@resadp1, "1_", "adapter,-1,-1,$junc", undef, 1, $lenMaxL, $lenMaxR, $name, $DATA->{$name}{gap}, $outBigLog, $outLog, $namewant);
		LOG($outBigLog, "1. endpos=$end_pos, adapter begdash=$begdash, dashes=$dashes\n"); #if $name eq $namewant;
		$adapter_TN_aln = $orig_TN_aln;
	}
	else {
	@cmdTN = muscle(
">seqTN_aln\n$seqTN
>seqAD_aln\n$seqADP
");
		($DATA->{$name}) = parse_fasta_simple($DATA->{$name}, \@cmdTN);
		$resadp1 = ">1_neg_$name\n$DATA->{$name}{'seqTN_aln'}\n>0_ADP_$name\n$DATA->{$name}{'seqAD_aln'}";
		@resadp1 = split("\n", $resadp1);
#		$resadp1 = "echo '>1_neg_$name\n$DATA->{$name}{'seqTN_aln'}\n>0_ADP_$name\n$seqADP\n'";
#		@resadp1 = `$resadp1 | muscle $muscleparam 2>/dev/null`;
#		my $currdef = "INIT";
#		my %temp;
#		for (my $k = 0; $k < @resadp1; $k++) {
#			chomp($resadp1[$k]);
#			if ($resadp1[$k] =~ /^>/) {
#				$currdef = $resadp1[$k]; $currdef =~ s/^>//;
#			}
#			else {
#				$temp{$currdef} .= $resadp1[$k];
#			}
#		}
#		my $seqwants = "";
#		foreach my $currdef2 (sort keys %temp) {
#			$seqwants = $temp{$currdef2} if $currdef2 =~ /1_neg/;
#		}
#		my ($newjunc) = fix_pos($junc, $DATA->{$name}{'seqTN_aln'}, $seqwants);
#		$newjunc += $DATA->{$name}{gap};

		$adapter_TN_aln = $DATA->{$name}{'seqTN_aln'};
		my ($newjunc) = fix_pos($junc, $orig_TN_aln, $adapter_TN_aln);
		($adapter_pos, $adapter_pos_fix, $end_pos, $begdash, $dashes) = parse_aln(\@resadp1, "1_", "adapter,-1,-1,$newjunc", undef, 1, $lenMaxL, $lenMaxR, $name, $DATA->{$name}{gap}, $outBigLog, $outLog, $namewant);
		my $orig_end_pos = $end_pos;
		($adapter_pos) = fix_pos($adapter_pos, $adapter_TN_aln, $orig_TN_aln);
		($adapter_pos_fix) = fix_pos($adapter_pos_fix, $adapter_TN_aln, $orig_TN_aln);
		($end_pos) = fix_pos($end_pos, $adapter_TN_aln, $orig_TN_aln) if $end_pos ne -99;

		LOG($outBigLog, "2. endpos=$LGN$end_pos$N (after fix_pos() from seqTN to seqTN_aln) orig endpos=$LPR$orig_end_pos$N, adapter begdash=$begdash, dashes=$dashes\n"); #if $name eq $namewant;
	}


	my ($igtype, $junc1pos, $junc2pos) = get_igtype($DATA->{$name}{junc1pos}, $DATA->{$name}{junc2pos}, $chr2, $beg2, $end2, $strand2, $outBigLog);
	$DATA->{$name}{igtype} = $igtype;

#	LOG($outBigLog, "endpos=$end_pos, adapter begdash=$begdash, dashes=$dashes\n"); #if $name eq $namewant;
	$beg_pos = -99 if not defined $beg_pos;
	$end_pos = -99 if not defined $end_pos;

	my $NA = $name eq $namewant ? undef : "NA";

	my ($seqTNL, $seqTNR) = ($seqTN, "");
	if (defined $junc) {
		($seqTNL, $seqTNR) = $seqTN =~ /^(.{$junc})(.+)$/;
	}

	my $nexted_reason = "";
	my $nexted = 0;
	if ($beg_pos eq -99 and $end_pos eq -99) {
		$nexted_bothshort ++;
		$nexted_reason .= "\n${LRD}NEXTED$N! INFO: $name strand1=$strand1 strand2=$strand2$LRD is bothshort $N(beg_pos=$beg_pos, end_pos=$end_pos, junc=$junc)\n\n";
		$nexted = 1;
	}
	elsif ($beg_pos eq -99) {
		$nexted_begshort ++;
		$nexted_reason .= "\n${LRD}NEXTED$N! INFO: $name strand1=$strand1 strand2=$strand2$LRD is begshort $N(beg_pos=$beg_pos, end_pos=$end_pos, junc=$junc)\n\n";
		$nexted = 1;
	}
	elsif ($end_pos eq -99) {
		$nexted_endshort ++;
		$nexted_reason .= "\n${LRD}NEXTED$N! INFO: $name strand1=$strand1 strand2=$strand2$LRD is endshort $N(beg_pos=$beg_pos, end_pos=$end_pos, junc=$junc)\n\n";
		$nexted = 1;
	}
	elsif ($junc < 25) {
		$nexted_juncmissing ++;
		$nexted_reason .= "\n${LRD}NEXTED$N! INFO: $name strand1=$strand1 strand2=$strand2$LRD is juncmissing $N(beg_pos=$beg_pos, end_pos=$end_pos, junc=$junc)\n\n";
		$nexted = 1;
	}
	#elsif ($igtype !~ /(IgM|IgG1|IgG3)$/) {
	#	$nexted_others ++;
	#	$nexted_reason .= "\n${LRD}NEXTED$N! INFO: $name strand1=$strand1 strand2=$strand2$LRD is not at IgM/G1/G3$N (igtype=$LCY$igtype$N, beg_pos=$beg_pos, end_pos=$end_pos, junc=$junc)\n\n";
	#	$nexted = 1;
	#}
	else {
		$okay{$name} = 1;
		$good ++;
		$nexted_reason .= "\nINFO: $name strand1=$strand1 strand2=$strand2$LGN is good$N (beg_pos=$beg_pos, end_pos=$end_pos, junc=$junc)\n\n";
	}




	my $unaligned_print = "";
	$unaligned_print .= "------------------------------------------------\n";
	$unaligned_print .= "------------------------------------------------\n";
	$unaligned_print .= "0. UNALIGNED SEQUENCES:\n";
	$unaligned_print .= "$YW$name$N\n";
	$unaligned_print .= "------------------------------------------------\n";
	$unaligned_print .= $nexted_reason;
	$unaligned_print .= "CONSENSUS  :\t      ------" . colorize($DATA->{$name}{seqC1L}) . colorize($DATA->{$name}{seqC2R}) . "\n";
	$unaligned_print .= "IgM BAIT   :\t      ------" . colorize($DATA->{$name}{seqC1L}) . "$seqC1R\n";
	$unaligned_print .= "FAR PREY   :\t      ------" . $seqC2L . colorize($DATA->{$name}{seqC2R}) . "\n\n";
	$unaligned_print .= "READ1 SEQ  :\t      " . colorize($seqQ1) . "\n";
	$unaligned_print .= "TLX SEQ    :\t      " . colorize($DATA->{$name}{seqTN}) . "\n";
#	DIELOG($outLog, "Length of seqTN is 0??\nlength seqQ2=" . length($seqQ2) . "\n" . colorize($seqQ2) . "\nlength seqTN=" . length($DATA->{$name}{seqTN}) . "\n" . colorize($DATA->{$name}{seqTN}) . "\n") if length($DATA->{$name}{seqTN}) - length($seqQ2) < 0;
	my $READSEQ2GAP = length($DATA->{$name}{seqTN}) - length($seqQ2); $READSEQ2GAP = 0 if $READSEQ2GAP < 0;
	$unaligned_print .= "READ2 SEQ  :\t      " . join("", (" ") x $READSEQ2GAP) . colorize($seqQ2) . "\n\n";
	$unaligned_print .= "TLX REF    :\t      " . colorize($DATA->{$name}{seqTR}) . "\n";
	$unaligned_print .= "MH SEQ ALN :\t      " . colorize($DATA->{$name}{mhseqlong}) . "\n";
	$unaligned_print .= "MH SEQ     :\t      " . colorize($DATA->{$name}{mhseq}) . "\n";
	$unaligned_print .= "GAP        :\t      $LGN" . $DATA->{$name}{gap} . "$N bp\n";
	my $mhlen = length($DATA->{$name}{mhseq}); $mhlen = 0 if not defined $mhlen;
	LOG($outBigLog, $unaligned_print);
	LOG($outLog, $unaligned_print, "NA");
#	LOG($outBigLog, "        THERE IS GAP ($YW$DATA->{$name}{gap}$N) SO ENDPOS GOTTEN USING MODIFIED JUNC (JUNC ($LGN$junc$N) + GAP ($YW$DATA->{$name}{gap}$N) = $LGN$end_pos$N)\n") if $DATA->{$name}{gap} > 0;
#	LOG($outLog, "        THERE IS GAP ($YW$DATA->{$name}{gap}$N) SO ENDPOS GOTTEN USING MODIFIED JUNC (JUNC ($LGN$junc$N) + GAP ($YW$DATA->{$name}{gap}$N) = $LGN$end_pos$N)\n","NA") if $DATA->{$name}{gap} > 0;

	next if $nexted eq 1;

	LOG($outBigLog, "HERE1\n");

	my $colorpos1 = $strand1 eq "+" ? "${LGN}$beg1${N}-$end1" : "$beg1-${LGN}$end1${N}";
	my $colorpos2 = $strand2 eq "+" ? "${LGN}$beg2${N}-$end2" : "$beg2-${LGN}$end2${N}";

	#my ($igtype, $junc1pos, $junc2pos) = get_igtype($DATA->{$name}{junc1pos}, $DATA->{$name}{junc2pos}, $chr2, $beg2, $end2, $strand2, $outBigLog);
#	$DATA->{$name}{igtype} = $igtype;
#	LOG($outBigLog, "\n\n" . date() . "${LGN}Aligning CONS, BAIT, READ, and PREY sequences!${N}\nIGTYPE=$LRD$igtype${N}, name=${LCY}$name${N}, beg0/junc/end0=$beg_pos/$junc/$end_pos (primer/adapter=$primer_pos/$adapter_pos, fixed=$primer_pos_fix/$adapter_pos_fix), (REF1=$chr1:$colorpos1;length1=$lenz1;junc1=${LGN}$begJ1${N}, REF2=$chr2:$colorpos2;length2=$lenz2;junc2=${LGN}$begJ2${N})\n\n",$NA);
#						my $igtypeprint = $igtype =~ /IgM$/ ? "IgM" : $igtype =~ /IgG1$/ ? "IgG1" : $igtype =~ /IgG3$/ ? "IgG3" : $igtype;
	my $lochash1;

	my $resTNs1 = 
">3_con_$name\n$DATA->{$name}{'seqCON_aln'}
>1_neg_$name\n$DATA->{$name}{'seqQ1_aln'}
>2_neg_$name\n$DATA->{$name}{'seqC1_aln'}
>3_neg_$name\n$DATA->{$name}{'seqC2_aln'}
";
	my @resTN1 = split("\n", $resTNs1);
	LOG($outBigLog, "\n\n" . date() . "${LGN}seqQ1:${N} Aligning CONS, BAIT, READ, and PREY sequences!${N} IGTYPE=$LRD$igtype${N}, name=${LCY}$name${N}, beg0/junc/end0=$beg_pos/$junc/$end_pos (primer/adapter=$primer_pos/$adapter_pos, fixed=$primer_pos_fix/$adapter_pos_fix), (REF1=$chr1:$colorpos1;length1=$lenz1;junc1=${LGN}$begJ1${N}, REF2=$chr2:$colorpos2;length2=$lenz2;junc2=${LGN}$begJ2${N})\nmhbeg=$DATA->{$name}{mhbeg}-$DATA->{$name}{mhend},mhpos=$DATA->{$name}{mhpos}\n\n",$NA);
	LOG($outLog, "\n\n" . date() . "${LGN}seqQ1:${N} Aligning CONS, BAIT, READ, and PREY sequences!${N} IGTYPE=$LRD$igtype${N}, name=${LCY}$name${N}, beg0/junc/end0=$beg_pos/$junc/$end_pos (primer/adapter=$primer_pos/$adapter_pos, fixed=$primer_pos_fix/$adapter_pos_fix), (REF1=$chr1:$colorpos1;length1=$lenz1;junc1=${LGN}$begJ1${N}, REF2=$chr2:$colorpos2;length2=$lenz2;junc2=${LGN}$begJ2${N})\nmhbeg=$DATA->{$name}{mhbeg}-$DATA->{$name}{mhend},mhpos=$DATA->{$name}{mhpos}\n\n","NA");
	my ($muthash, $typehash1, $typehash2);
	my ($beg_pos2, $end_pos2) = ($beg_pos, $end_pos);
	my ($beg_fix, $end_fix) = ($beg_pos, $end_pos);
	($muthash, $typehash1, $beg_pos2, $end_pos2) = parse_aln(\@resTN1, "1_", "none,$beg_pos,$end_pos,$DATA->{$name}{junc},$primer_pos_fix,$adapter_pos_fix,$DATA->{$name}{mhbeg},$DATA->{$name}{mhend}", $muthash, 1, $lenMaxL, $lenMaxR, $name, $DATA->{$name}{gap}, $outBigLog, $outLog, $namewant);
	$beg_fix = $beg_pos2 if $beg_fix < $beg_pos2;
	$end_fix = $end_pos2 if $end_fix > $end_pos2;
	if ($DATA->{$name}{gap} > 0) {
		$typehash1->{ins}{$DATA->{$name}{seqTNG}} = 1;
		push(@{$lochash1->{ins}{$DATA->{$name}{seqTNG}}{arr}}, $DATA->{$name}{junc});
	}
#	if ($DATA->{$name}{mh} > 0) {
#		$typehash1->{mh}{$DATA->{$name}{seqTNG}} = 1;
#		push(@{$lochash1->{mh}{$DATA->{$name}{seqTNG}}{arr}}, $DATA->{$name}{junc});
#	}

=comment
	my $resTNs2 = 
">3_con_$name\n$DATA->{$name}{'seqCON_aln'}
>1_neg_$name\n$DATA->{$name}{'seqQ2_aln'}
>2_neg_$name\n$DATA->{$name}{'seqC1_aln'}
>3_neg_$name\n$DATA->{$name}{'seqC2_aln'}
";
	my @resTN2 = split("\n", $resTNs2);
	LOG($outBigLog, "\n\n" . date() . "${LGN}seqQ2:${N} Aligning CONS, BAIT, READ, and PREY sequences!${N} IGTYPE=$LRD$igtype${N}, name=${LCY}$name${N}, beg0/junc/end0=$beg_pos/$junc/$end_pos (primer/adapter=$primer_pos/$adapter_pos, fixed=$primer_pos_fix/$adapter_pos_fix), (REF1=$chr1:$colorpos1;length1=$lenz1;junc1=${LGN}$begJ1${N}, REF2=$chr2:$colorpos2;length2=$lenz2;junc2=${LGN}$begJ2${N})\nmhbeg=$DATA->{$name}{mhbeg}-$DATA->{$name}{mhend},mhpos=$DATA->{$name}{mhpos}\n\n",$NA);
	LOG($outLog, "\n\n" . date() . "${LGN}seqQ2:${N} Aligning CONS, BAIT, READ, and PREY sequences!${N} IGTYPE=$LRD$igtype${N}, name=${LCY}$name${N}, beg0/junc/end0=$beg_pos/$junc/$end_pos (primer/adapter=$primer_pos/$adapter_pos, fixed=$primer_pos_fix/$adapter_pos_fix), (REF1=$chr1:$colorpos1;length1=$lenz1;junc1=${LGN}$begJ1${N}, REF2=$chr2:$colorpos2;length2=$lenz2;junc2=${LGN}$begJ2${N})\nmhbeg=$DATA->{$name}{mhbeg}-$DATA->{$name}{mhend},mhpos=$DATA->{$name}{mhpos}\n\n","NA");
	($muthash, $typehash2, $beg_pos2, $end_pos2) = parse_aln(\@resTN2, "1_", "none,$beg_pos,$end_pos,$DATA->{$name}{junc},$primer_pos_fix,$adapter_pos_fix,$DATA->{$name}{mhbeg},$DATA->{$name}{mhend}", $muthash, 2, $lenMaxL, $lenMaxR, $name, $DATA->{$name}{gap}, $outBigLog, $outLog, $namewant);
	$beg_fix = $beg_pos2 if $beg_fix < $beg_pos2;
	$end_fix = $end_pos2 if $end_fix > $end_pos2;

=cut

	LOG($outBigLog, "HERE2\n");

	LOG($outBigLog, "${YW}------------------------------------------------\n6a. Tabulating mutation from seqQ1 (Read pair #1)\n------------------------------------------------${N}\n" . date() . "\n",$NA);
	LOG($outLog, "${YW}------------------------------------------------\n6a. Tabulating mutation from seqQ1 (Read pair #1)\n------------------------------------------------${N}\n" . date() . "\n","NA");

	foreach my $mutpos (sort keys %{$muthash}) {
		my $mutposcurr = $mutpos;
		my $mut1 = $muthash->{$mutpos}{1}; $mut1 = "none" if not defined $mut1;
		my ($type1, $type2) = $mut1 =~ /^(ins|del|mat|mis|mh)_(.+)$/;
		if ($type1 eq "mh") {
			$type2 = $DATA->{$name}{mhseq};
			$mutposcurr = $DATA->{$name}{mhbeg};
		}
		DIELOG($outBigLog, "Can't find type1 from mutpos=$mutpos mutposcurr=$mutposcurr mut=$mut1; not ins/del/mat/mis/mh_?\n\n") if not defined $type1;
		push(@{$lochash1->{$type1}{$type2}{arr}}, $mutposcurr);
	}
	foreach my $type1 (sort keys %{$lochash1}) {
		foreach my $type2 (sort keys %{$lochash1->{$type1}}) {
			my @arr = @{$lochash1->{$type1}{$type2}{arr}};
			@arr = sort {$a <=> $b} @arr;
			$lochash1->{$type1}{$type2}{join} = join(",", @arr);
			LOG($outBigLog, "$type1\t$type2\t$lochash1->{$type1}{$type2}{join}\n");
		}
	}
	my $begorig1 = $DATA->{$name}{begorig1};
	my $endorig1 = $DATA->{$name}{endorig1};
	my $begorig2 = $DATA->{$name}{begorig2};
	my $endorig2 = $DATA->{$name}{endorig2};

	my $seqCON_orig = $DATA->{$name}{'seqCON_aln'};
	$seqCON_orig =~ s/\-//g;

	my %temp; my @temp; my $lastmutpos = -1;
	my $juncposfix = fix_pos2($DATA->{$name}{junc}, $DATA->{$name}{'seqCON_aln'}, $seqCON_orig, $outBigLog);
	$temp[$juncposfix] = "J";

#
#my $outzfixFile = $outFile . ".outzfix.bed.temp";
#my $outztempFile = $outFile . ".outztemp.bed.temp";
#my $outzfaFile = $outFile . ".outztemp.fa";

	my %coorz;
	my $begB6 = 114485808;
	my $endB6 = 117249165;
	my $IgMbeg = 114656770;
	my $IgMend = 114670045;
	my $IgG1beg = 114563452;
	my $IgG1end = 114582573;
	my $IgG3beg = 114594436;
	my $IgG3end = 114609628;
	my $printzfix;

	LOG($outBigLog, "HERE3\n");

# HERE

	my %fix;
	$fix{try} = 0;
	$fix{changed} = 0;
	$fix{total} = 0;
	$fix{add} = 0; # add to begind1
	$fix{best} = "INIT";
	$fix{bestchanged} = "INIT";
#	my @add = (0, -1, 1, -2, 2, -3, 3, -4, 4, -5, 5,-6, 6, -7, 7);
	my $fix_add = 0;

#=comment
	for (my $try = 0; $try < @add; $try ++){
		$fix{add} = $add[$try];
		foreach my $type1 (sort keys %{$lochash1}) {
			foreach my $type2 (sort keys %{$lochash1->{$type1}}) {
				my @arr = @{$lochash1->{$type1}{$type2}{arr}};
				@arr = sort {$a <=> $b} @arr;
				$lochash1->{$type1}{$type2}{join} = join(",", @arr);

#				LOG($outBigLog, "$type1\t$type2\t$lochash1->{$type1}{$type2}{join}\n");

				my $mutposcurr = fix_pos2($lochash1->{$type1}{$type2}{join}, $DATA->{$name}{'seqCON_aln'}, $seqCON_orig);#, $outBigLog);
				my $nuc = ".";
				if ($type1 eq "mat" or $type1 eq "mis") {
					my ($nuc1, $nuc2) = split("_", $type2);
					$nuc = $nuc1;
					$nuc = uc($nuc1) if $type1 eq "mat";
					$nuc = lc($nuc2) if $type1 eq "mis";
				}
				else {
					$nuc = "I" if $type1 eq "ins";
					$nuc = $type2 if $type1 eq "del";
					$nuc = "H" if $type1 eq "mh";
				}
				my @mutposcurr = split(",", $mutposcurr);

				for (my $ki = 0; $ki < @mutposcurr; $ki++) {
					my $ind = $mutposcurr[$ki];
					if ($ind > $juncposfix) { #not IGM
					LOG($outBigLog, "\tki=$LCY$ki$N, actual pos = $LPR$ind$N, $LGN try # $fix{try}:$YW Trying add=$fix{add}$N\n","NA");
						$fix{total} ++;
						my $begind2 = $strand2 eq "+" ? $begorig2 + ($ind - $juncposfix + $mhlen) - 1: $endorig2 - ($ind - $juncposfix - 1) - 1;
						my $endind2 = $strand2 eq "+" ? $begorig2 + ($ind - $juncposfix + $mhlen) - 0: $endorig2 - ($ind - $juncposfix - 1);
#						$coorz{ind2} = $ind if not defined $coorz{ind2};
#						$coorz{ind2} = $ind if $ind < $coorz{ind2};
#						$coorz{begind2} = $begind2 if not defined $coorz{begind2};
#						$coorz{begind2} = $begind2 if $begind2 < $coorz{begind2};
#						$coorz{endind2} = $endind2 if not defined $coorz{endind2};
#						$coorz{endind2} = $endind2 if $endind2 > $coorz{endind2};
						my $igtypeprint = $igtype =~ /IgM$/ ? "IgM" : $igtype =~ /IgG1$/ ? "IgG1" : $igtype =~ /IgG3$/ ? "IgG3" : $igtype;
						my $endind2print = $igtypeprint eq "IgM" ? $IgMend - $begind2 : $igtypeprint eq "IgG1" ? $IgG1end - $begind2 : $igtypeprint eq "IgG3" ? $IgG3end - $begind2 : $begind2;
						my $begind2print = $igtypeprint eq "IgM" ? $IgMend - $endind2 : $igtypeprint eq "IgG1" ? $IgG1end - $endind2 : $igtypeprint eq "IgG3" ? $IgG3end - $endind2 : $endind2;
						my $igmtouse = $igtypeprint eq "IgM" ? "IGMEND=IgMend-$begind2" : $igtypeprint eq "IgG1" ? "IGG1end=$IgG1end-$begind2" : $igtypeprint eq "IgG3" ? "IGG3end=$IgG3end-$begind2" : "UNKNOWN";
						my $constype1 = $cons{$igtypeprint}{1}{$begind2print+$fix{add}}{type1}; $constype1 = "UNDEF_constype1" if not defined $constype1;
						my $constype2 = $cons{$igtypeprint}{1}{$begind2print+$fix{add}}{type2}; $constype2 = "UNDEF_constype2" if not defined $constype2;
						my ($cons_fix, $printzfix2);
						$printzfix2 = "";
						($cons_fix, $constype2, $printzfix2) = check_cons_fix($type1, $type2, $constype1, $constype2, $printzfix2, \%cons, $begind2print, $igtypeprint, $strand2, $fix{add});
						if ($cons_fix =~ /(CHANGE|BAD)/) {
							$fix{changed} ++;
						}
						LOG($outBigLog, "\tFIX AJ851868.3\t$begind2print\t$endind2print\t$igtypeprint;$type1;$type2;$ind\t0\t$strand2\t$igtypeprint\t$constype1,$constype2\t$cons_fix\tfix try=$fix{try}, fix_add=$fix{add}\n","NA");
						LOG($outBigLog, "\tcons fix = $cons_fix\n","NA");
					}
				}
			}
		}
		if ($fix{best} eq "INIT") {
			$fix{best} = $fix{add};
			$fix{bestchange} = $fix{changed};
			LOG($outBigLog, "0. fix best = $fix{best}, fix add=$fix{add}, bestchanged= $fix{bestchange} < changed=$fix{changed}}\n");
		}
		elsif ($fix{bestchange} > $fix{changed}) {
			$fix{best} = $fix{add};
			$fix{bestchange} = $fix{changed};
			LOG($outBigLog, "1. fix best = $fix{best}, fix add=$fix{add}, bestchanged= $fix{bestchange} < changed=$fix{changed}}\n");
		}
		else {
			LOG($outBigLog, "2. try=$try, fix best = $fix{best}, fix add=$fix{add}, bestchanged= $fix{bestchange} < changed=$fix{changed}}\n") if $try % 20 == 0;
		}
		if ($try == @add - 1 or $fix{total} < 3 or ($fix{total} > 0 and $fix{changed} / $fix{total} <= 0.1)) {
			my $perc = $fix{total} == 0 ? "N/A" : int(1000 * $fix{changed}/$fix{total})/10;
			LOG($outBigLog, "\t$LGN GOOD at try # $LPR$fix{try}$N!\n","NA") if $try ne @add - 1;
			LOG($outBigLog, "\t$LGN CANT FIND at try # $LPR$fix{try}$N!\n","NA") if $try eq @add - 1;
			LOG($outBigLog, "\t  fix total = $LGN$fix{total}$N\n","NA");
			LOG($outBigLog, "\t  fix changed = $LGN$fix{changed}$N\n","NA");
			LOG($outBigLog, "\t  fix perc = $LGN$perc$N\n\n","NA");
			LOG($outBigLog, "fix best = $fix{best}, changed= $fix{bestchange}\n");
			last;
		}
		else {
			$fix{try} ++;
			$fix{changed} = 0;
			$fix{total} = 0;
		}
	}
	$fix{add} = $fix{best};
	$fix_add = $fix{best};
	my $fix_perc = $fix{total} == 0 ? "${YW}N/A$N" : int(1000*$fix{changed}/$fix{total})/10;
	LOG($outBigLog, ">FIX_ADD $YW$name$N fix_add = $fix_add, changed = $fix{changed}, try=$fix{try}, perc changed = $LGN$fix_perc$N\n");
#=cut



# dump typehash1
my $printhash1 = "";
foreach my $type1 (sort keys %{$typehash1}) {
	foreach my $type2 (sort keys %{$typehash1->{$type1}}) {
		if ($typehash1->{$type1}{$type2} =~ /HASH/) {
			my $typehash1temp = $typehash1->{$type1}{$type2};
			foreach my $key (sort keys %{$typehash1temp}) {
				$printhash1 .= "type1=$LCY$type1$N, type2=$LGN$type2$N, key=$key value=$LPR$typehash1->{$type1}{$type2}{$key}$N\n";
			}
		}
		else {
			$printhash1 .= "type1=$LCY$type1$N, type2=$LGN$type2$N, value=$LPR$typehash1->{$type1}{$type2}$N\n";
		}
	}
}
$printhash1 .= "\n\n";

	my $lochashfinal;
	my $maxind = 0;
	foreach my $type1 (sort keys %{$lochash1}) {
		foreach my $type2 (sort keys %{$lochash1->{$type1}}) {
			my @arr = @{$lochash1->{$type1}{$type2}{arr}};
			@arr = sort {$a <=> $b} @arr;
			$lochash1->{$type1}{$type2}{join} = join(",", @arr);

			#LOG($outBigLog, "FIX_POS2\n");
			my $mutposcurr = fix_pos2($lochash1->{$type1}{$type2}{join}, $DATA->{$name}{'seqCON_aln'}, $seqCON_orig);#, $outBigLog);
			#LOG($outBigLog, "$type1\t$type2\t$lochash1->{$type1}{$type2}{join}\t$mutposcurr\n");

			my $nuc = ".";
			if ($type1 eq "mat" or $type1 eq "mis") {
				my ($nuc1, $nuc2) = split("_", $type2);
				$nuc = $nuc1;
				$nuc = uc($nuc1) if $type1 eq "mat";
				$nuc = lc($nuc2) if $type1 eq "mis";
			}
			else {
				$nuc = "I" if $type1 eq "ins";
				$nuc = $type2 if $type1 eq "del";
				$nuc = "H" if $type1 eq "mh";
			}

			my @mutposorig = @arr;
			my @mutposorig2 = @arr;
			my @mutposcurr = split(",", $mutposcurr);
			my @mutposcurr2 = split(",", $mutposcurr);
			#my $fix_add = $fix{add};
			for (my $ki = 0; $ki < @mutposorig; $ki++) {
				$mutposorig2[$ki] += $fix_add;
				$mutposcurr2[$ki] += $fix_add;
			}
			#$lochash1->{$type1}{$type2}{arr} = \@mutposorig2;
			#LOG($outBigLog, "mutposorig before = " . join(",", @mutposorig) . ";$LGN fix_add=$fix_add$N after = " . join(",", @mutposorig2) . "\n");
			#LOG($outBigLog, "mutposcurr before = " . join(",", @mutposcurr) . ";$LGN fix_add=$fix_add$N after = " . join(",", @mutposcurr2) . "\n");
			for (my $ki = 0; $ki < @mutposcurr; $ki++) {
				my $ind = $mutposcurr[$ki];
				$maxind = $ind if $maxind < $ind;
			}
		}
	}

	foreach my $type1 (sort keys %{$lochash1}) {
		foreach my $type2 (sort keys %{$lochash1->{$type1}}) {
			my @arr = @{$lochash1->{$type1}{$type2}{arr}};
			@arr = sort {$a <=> $b} @arr;
			$lochash1->{$type1}{$type2}{join} = join(",", @arr);

			#LOG($outBigLog, "FIX_POS2\n");
			my $mutposcurr = fix_pos2($lochash1->{$type1}{$type2}{join}, $DATA->{$name}{'seqCON_aln'}, $seqCON_orig);#, $outBigLog);
			#LOG($outBigLog, "$type1\t$type2\t$lochash1->{$type1}{$type2}{join}\t$mutposcurr\n");

			my $nuc = ".";
			if ($type1 eq "mat" or $type1 eq "mis") {
				my ($nuc1, $nuc2) = split("_", $type2);
				$nuc = $nuc1;
				$nuc = uc($nuc1) if $type1 eq "mat";
				$nuc = lc($nuc2) if $type1 eq "mis";
			}
			else {
				$nuc = "I" if $type1 eq "ins";
				$nuc = $type2 if $type1 eq "del";
				$nuc = "H" if $type1 eq "mh";
			}

			my @mutposorig = @arr;
			my @mutposorig2 = @arr;
			my @mutposcurr = split(",", $mutposcurr);
			my @mutposcurr2 = split(",", $mutposcurr);
			my $fix_add = $fix{add};
			for (my $ki = 0; $ki < @mutposorig; $ki++) {
				$mutposorig2[$ki] += $fix_add;
				$mutposcurr2[$ki] += $fix_add;
			}
			$lochash1->{$type1}{$type2}{arr} = \@mutposorig2;
			LOG($outBigLog, "mutposorig before = " . join(",", @mutposorig) . ";$LGN fix_add=$fix_add$N after = " . join(",", @mutposorig2) . "\n");
			LOG($outBigLog, "mutposcurr before = " . join(",", @mutposcurr) . ";$LGN fix_add=$fix_add$N after = " . join(",", @mutposcurr2) . "\n");
			my $lastbegind1print = "";
#				LOG($outBigLog, ">$LCY ki = $ki = $ind$N\n");
			for (my $ki = 0; $ki < @mutposcurr; $ki++) {
				my $ind = $mutposcurr[$ki];
#				LOG($outBigLog, ">$LCY ki = $ki = $ind$N\n");
				if ($ind <= $juncposfix) { #IGM
					my $begind1 = $strand1 eq "+" ? $begorig1 + $ind : $endorig1 - $ind - 1;
					my $endind1 = $strand1 eq "+" ? $begorig1 + $ind + 1 : $endorig1 - $ind;
					$coorz{ind1} = $ind if not defined $coorz{ind1};
					$coorz{ind1} = $ind if $ind < $coorz{ind1};
					$coorz{begind1} = $begind1 if not defined $coorz{begind1};
					$coorz{begind1} = $begind1 if $begind1 < $coorz{begind1};
					$coorz{endind1} = $endind1 if not defined $coorz{endind1};
					$coorz{endind1} = $endind1 if $endind1 > $coorz{endind1};
					my $endind1print = $IgMend - $begind1;
					my $begind1print = $IgMend - $endind1;
					if ($strand1 eq "+") {
						$begind1print = $IgMbeg + $begind1;
						$endind1print = $IgMbeg + $endind1;
					}
#					my $endind1print = $endind1 - $IgMbeg;#$IgMend - $begind1;
#					my $begind1print = $begind1 - $IgMbeg;#$IgMend - $endind1;
#					my $begind1print = $IgMend - $begind1;# - $begB6;
#					my $endind1print = $IgMend - $endind1;# - $begB6;
#					my $begind1print = $endB6 - $begind1;
#					my $endind1print = $endB6 - $endind1;
#					$begind1print -= $begB6;
#					$endind1print -= $begB6;
					#if ($tyep1 ne "mat") {
					if (defined $printzfix->{$begind1print} and $printzfix->{$begind1print} !~ /;ins/) {
#						$printzfix->{$begind1print} .= "\ni=$ki, mutposcurr=$mutposcurr[$ki], ind=$ind, type1=$type1, type2=$type2, Already defined begind1print=$begind1print above!\n" if $strand1 eq "+";
#						LOG($outBigLog, "\ni=$ki, mutposcurr=$mutposcurr[$ki], ind=$ind, type1=$type1, type2=$type2, Already defined begind1print=$begind1print above!\n");
					}
					else {
#						$printzfix->{$begind1print} .= "1 beg1indprint=$begind1print, i=$ki, mutposcurr=$mutposcurr[$ki], ind=$ind, type1=$type1, type2=$type2\n" if $strand1 eq "+";
#						LOG($outBigLog, "\ni=$ki, mutposcurr=$mutposcurr[$ki], ind=$ind, type1=$type1, type2=$type2\n");
					}
#					my $igtypeprint = $igtype =~ /IgM$/ ? "IgM" : $igtype =~ /IgG1$/ ? "IgG1" : $igtype =~ /IgG3$/ ? "IgG3" : $igtype;


					my @checkmh = (0,-1,-2,1,2);
					my ($constype1, $constype2, $cons_fix, $printzfix2);
					if ($type1 eq "mh") {
						my $mh_fix_add = 0;
						for (my $kj = 0; $kj < @checkmh; $kj++) {
							$constype1 = $cons{IgM}{1}{$begind1print+$checkmh[$kj]}{type1}; $constype1 = "UNDEF_constype1" if not defined $constype1;
							$constype2 = $cons{IgM}{1}{$begind1print+$checkmh[$kj]}{type2}; $constype2 = "UNDEF_constype2" if not defined $constype2;
							next if not defined $constype1;
							$printzfix2 = $printzfix->{$begind1print};
							($cons_fix, $constype2, $printzfix2) = check_cons_fix($type1, $type2, $constype1, $constype2, $printzfix2, \%cons, $begind1print, "IgM", $strand1, $checkmh[$kj]);
							my $origcons_fix = $cons_fix; $origcons_fix =~ s/[ ]+/_/g;
							$cons_fix = $cons_fix =~ /GOOD/ ? $LGN . $cons_fix . $N : $cons_fix =~ /BAD/ ? $LRD . $cons_fix . $N : $LCY. $cons_fix . $N;
							if ($cons_fix =~ /GOOD/) {
								$mh_fix_add = $checkmh[$kj];
								$printzfix->{$begind1print} = $printzfix2;
								$printzfix->{$begind1print} .= "NG_005838.1\t$begind1print\t$endind1print\tIgM;$type1;$type2;$ind\t0\t$strand1\tIgM\t$constype1,$constype2\t$cons_fix\tfix_add=$mh_fix_add\n";
#								print $outloc "$name\tBAIT\tNG_005838.1\t$begind1print\t$endind1print\tIgM;$type1;$type2;$ind\t0\t$strand1\tIgM\t$constype1\t$constype2\t$origcons_fix\t$mh_fix_add\n";
								print $outloc "NG_005838.1_IgM\t$begind1print\t$endind1print\t$type1\_$type2\t0\t$strand1\tmutorig=$type1\_$type2;change=$origcons_fix;constype=$constype1\_$constype2;name=$name;type=BAIT;ig=IgM;fix_add=$mh_fix_add;ki=$ki;mutposorig=$mutposorig2[$ki];mutposcurr=$mutposcurr2[$ki];coor=chr12\_$begind1\_$endind1\_$name,IgM\_0\_$strand1\n";
								last;
							}
							else {
#								$printzfix->{$begind1print} .= "\tMH i=$kj, type1_type2=$type1\_$type2, constype1=$constype1, constype2=$constype2, cons_fix=$LPR$cons_fix$N,checkmh=$LPR$checkmh[$kj]$N\n";
							}
						}
#						$printzfix->{$begind1print} .= "NG_005838.1\t$begind1print\t$endind1print\tIgM;$type1;$type2;$ind\t0\t$strand1\tIgM\t$constype1,$constype2\t$cons_fix\tfix_add=$mh_fix_add\n";
					}
					else {
						$constype1 = $cons{IgM}{1}{$begind1print}{type1}; $constype1 = "UNDEF_constype1" if not defined $constype1;
						$constype2 = $cons{IgM}{1}{$begind1print}{type2}; $constype2 = "UNDEF_constype2" if not defined $constype2;
						$printzfix2 = $printzfix->{$begind1print};
						($cons_fix, $constype2, $printzfix2) = check_cons_fix($type1, $type2, $constype1, $constype2, $printzfix2, \%cons, $begind1print, "IgM", $strand1, 0);
						my ($change, $ctype1, $ctype2) = ("", $type1, $type2);
						if ($cons_fix =~ /CHANGE/ and $cons_fix !~ /DONT/) {
							($change, $ctype1, $ctype2) = $cons_fix =~ /^(CHANGE) ([A-Za-z0-9]+)_(.+)$/;
							$printzfix2 .= "\t$LRD ctype1 or 2 undefined! cons_fix=$cons_fix\n" if not defined $ctype1 or not defined $ctype2;
							$ctype1 = "NA" if not defined $ctype1; 
							$ctype2 = "NA" if not defined $ctype2;
							$printzfix2 .= "\t$YW>Changing $LGN$type1 $type2$N into $LPR$ctype1 $ctype2$N\n" if defined $ctype1 and defined $ctype2;
							$typehash1->{$type1}{$type2} --;
							$typehash1->{$ctype1}{$ctype2} ++;
							if ($typehash1->{$type1}{$type2} == 0) {
								undef $typehash1->{$type1}{$type2};
								delete $typehash1->{$type1}{$type2};
							}
						}
						my $origcons_fix = $cons_fix; $origcons_fix =~ s/[ ]+/_/g;
						$cons_fix = $cons_fix =~ /GOOD/ ? $LGN . $cons_fix . $N : $cons_fix =~ /BAD/ ? $LRD . $cons_fix . $N : $LCY. $cons_fix . $N;
						$printzfix->{$begind1print} = $printzfix2;
						$printzfix->{$begind1print} .= "NG_005838.1\t$begind1print\t$endind1print\tIgM;$type1;$type2;$ind\t0\t$strand1\tIgM\t$constype1,$constype2\t$cons_fix\t0\n";
#						print $outloc "$name\tBAIT\tNG_005838.1\t$begind1print\t$endind1print\tIgM;$type1;$type2;$ind\t0\t$strand1\tIgM\t$constype1\t$constype2\t$origcons_fix\n";
						print $outloc "NG_005838.1_IgM\t$begind1print\t$endind1print\t$ctype1\_$ctype2\t0\t$strand1\tmutorig=$type1\_$type2;change=$origcons_fix;constype=$constype1\_$constype2;name=$name;type=BAIT;ig=IgM;fix_add=0;ki=$ki;mutposorig=$mutposorig2[$ki];mutposcurr=$mutposcurr2[$ki];coor=chr12\_$begind1\_$endind1\_$name,IgM\_0\_$strand1\n";
					}
					$lastbegind1print = $begind1print;
#					LOG($outBigLog, "NG_005838.1\t$begind1print\t$endind1print\t$type1;$type2;$ind\t0\t$strand2\n");
					#}
				}
				else {
					my $begind2 = $strand2 eq "+" ? $begorig2 + ($ind - $juncposfix + $mhlen) - 1: $endorig2 - ($ind - $juncposfix - 1) - 1;
					my $endind2 = $strand2 eq "+" ? $begorig2 + ($ind - $juncposfix + $mhlen) - 0: $endorig2 - ($ind - $juncposfix - 1);
					$coorz{ind2} = $ind if not defined $coorz{ind2};
					$coorz{ind2} = $ind if $ind < $coorz{ind2};
					$coorz{begind2} = $begind2 if not defined $coorz{begind2};
					$coorz{begind2} = $begind2 if $begind2 < $coorz{begind2};
					$coorz{endind2} = $endind2 if not defined $coorz{endind2};
					$coorz{endind2} = $endind2 if $endind2 > $coorz{endind2};
					my $igtypeprint = $igtype =~ /IgM$/ ? "IgM" : $igtype =~ /IgG1$/ ? "IgG1" : $igtype =~ /IgG3$/ ? "IgG3" : $igtype;
					my $endind2print = $igtypeprint eq "IgM" ? $IgMend - $begind2 : $igtypeprint eq "IgG1" ? $IgG1end - $begind2 : $igtypeprint eq "IgG3" ? $IgG3end - $begind2 : $begind2;
					my $begind2print = $igtypeprint eq "IgM" ? $IgMend - $endind2 : $igtypeprint eq "IgG1" ? $IgG1end - $endind2 : $igtypeprint eq "IgG3" ? $IgG3end - $endind2 : $endind2;
					my $igmtouse = $igtypeprint eq "IgM" ? "IGMEND=IgMend-$begind2" : $igtypeprint eq "IgG1" ? "IGG1end=$IgG1end-$begind2" : $igtypeprint eq "IgG3" ? "IGG3end=$IgG3end-$begind2" : "UNKNOWN";
					if (defined $printzfix->{$begind2print} and $printzfix->{$begind2print} !~ /;ins/) {
#						$printzfix->{$begind2print} .= "\ni=$ki, mutposcurr=$mutposcurr[$ki], ind=$ind, type1=$type1, type2=$type2, Already defined begind2print=$begind2print above!\n" if $strand2 eq "+";
#						LOG($outBigLog, "\ni=$ki, mutposcurr=$mutposcurr[$ki], ind=$ind, type1=$type1, type2=$type2, Already defined begind2print=$begind2print above!\n");
					}
					else {
#						$printzfix->{$begind2print} .= "2 beg2=$begind2print i=$ki, mutposcurr=$mutposcurr[$ki], ind=$ind, type1=$type1, type2=$type2\n" if $strand2 eq "+";
#						LOG($outBigLog, "\ni=$ki, mutposcurr=$mutposcurr[$ki], ind=$ind, type1=$type1, type2=$type2\n");
					}

					
					if (defined $DATA->{$name}{mhseq} and length($DATA->{$name}{mhseq}) > 0  and $ind == $juncposfix + 1) {
#(($strand2 eq "+" and $ind == $juncposfix + 1) or ($strand2 eq "+" and $ind == $maxind))) {
						my $mhtype1 = "mh";
						my $mhtype2 = $DATA->{$name}{mhseq};
						if ($strand2 eq "+") {
							$mhtype2 =~ tr/ACGTacgt/TGCAtgca/;
						}
						my $lengthmhseq = length($DATA->{$name}{mhseq});
						my $begind2printmh = $strand2 eq "-" ? $begind2print - $lengthmhseq : $begind2print + $lengthmhseq;
						my $endind2printmh = $strand2 eq "-" ? $begind2print - $lengthmhseq + 1 : $begind2print + $lengthmhseq + 1;
#						$printzfix->{$begind2printmh} .= "\tMHSEQHERE! $DATA->{$name}{mhseq}, begind2printmh=$begind2printmh\n";
						my $mh_fix_add = "NA";
						my @checkmh = (0);#,-1,-2,-3,-4,-5,-6,-7,-8,-9,-10,1,2,3,4,5,6,7,8,9,10);
						my ($constype1, $constype2, $cons_fix, $printzfix2);
						for (my $kj = 0; $kj < @checkmh; $kj++) {
							$constype1 = $cons{$igtypeprint}{1}{$begind2printmh+$fix_add+$checkmh[$kj]}{type1}; $constype1 = "UNDEF_constype1" if not defined $constype1;
							$constype2 = $cons{$igtypeprint}{1}{$begind2printmh+$fix_add+$checkmh[$kj]}{type2}; $constype2 = "UNDEF_constype2" if not defined $constype2;
							if (not defined $constype1) {
								#$printzfix->{$begind2printmh} .= "undef constype1, nexted!\n";
								next;
							}
							else {
								#$printzfix->{$begind2printmh} .= "$kj: defined constype=$constype1\_$constype2, checkmh=$checkmh[$kj]\n";
							}
							$printzfix2 = $printzfix->{$begind2printmh};
							($cons_fix, $constype2, $printzfix2) = check_cons_fix($mhtype1, $mhtype2, $constype1, $constype2, $printzfix2, \%cons, $begind2printmh, $igtypeprint, $strand1, $checkmh[$kj]);
							my $origcons_fix = $cons_fix; $origcons_fix =~ s/[ ]+/_/g;
							$cons_fix = $cons_fix =~ /GOOD/ ? $LGN . $cons_fix . $N : $cons_fix =~ /BAD/ ? $LRD . $cons_fix . $N : $LCY. $cons_fix . $N;
							if ($cons_fix =~ /GOOD/) {
								$mh_fix_add = $checkmh[$kj];
								$printzfix->{$begind2printmh} = $printzfix2;
						
								$printzfix->{$begind2printmh} .= "AJ851868.3\t$begind2printmh\t$endind2printmh\t$igtypeprint;$mhtype1;$mhtype2;$ind\t0\t$strand2\t$igtypeprint\t$constype1,$constype2\t$cons_fix\tfix_add=$mh_fix_add\n";
#								print $outloc "$name\tBAIT\tNG_005838.1\t$begind2printmh\t$endind2printmh\t$igtypeprint;$mhtype1;$mhtype2;$ind\t0\t$strand2\t$igtypeprint\t$constype1\t$constype2\t$origcons_fix\t$mh_fix_add\n";
								print $outloc "AJ851868.3\t$igtypeprint\t$begind2printmh\t$endind2printmh\t$mhtype1\_$mhtype2\t0\t$strand1\tmutorig=$type1\_$type2;change=$origcons_fix;constype=$constype1\_$constype2;name=$name;type=BAIT;ig=$igtypeprint;fix_add=$mh_fix_add;ki=$ki;mutposorig=$mutposorig2[$ki];mutposcurr=$mutposcurr2[$ki];coor=chr12\_$begind2\_$endind2\_$name,$igtypeprint\_0\_$strand2\n";
								last;
							}
							else {
								$printzfix->{$begind2printmh} .= "\tMH not good! MH i=$kj, type1_type2=$type1\_$type2, constype1=$constype1, constype2=$constype2, cons_fix=$LPR$cons_fix$N,checkmh=$LPR$checkmh[$kj]$N\n";
							}
						}
						if ($mh_fix_add eq "NA") {
#							$printzfix->{$lastbegind1print} .= "\t$cons_fix\n";
							$printzfix->{$begind2printmh} .= "ki=$ki, ind=$ind, maxind=$maxind, AJ851868.3\t$begind2printmh\t$endind2printmh\t$igtypeprint;$mhtype1;$mhtype2;$ind\t0\t$strand2\t$igtypeprint\t$constype1,$constype2\t$cons_fix\tfix_add=$mh_fix_add\n";
							$printzfix->{$lastbegind1print} .= $printzfix->{$begind2printmh};
							$printzfix->{$begind2printmh} .= "";
							undef $printzfix->{$begind2printmh};
							delete $printzfix->{$begind2printmh};
#							$printzfix->{$begind2printmh} .= "AJ851868.3\t$begind2printmh\t$endind2printmh\t$igtypeprint;$mhtype1;$mhtype2;$ind\t0\t$strand2\t$igtypeprint\t$constype1,$constype2\t$cons_fix\tfix_add=$mh_fix_add\n";
							print $outloc "AJ851868.3\t$igtypeprint\t$begind2printmh\t$endind2printmh\t$mhtype1\_$mhtype2\t0\t$strand1\tmutorig=$type1\_$type2;change=$cons_fix;constype=$constype1\_$constype2;name=$name;type=BAIT;ig=$igtypeprint;fix_add=$mh_fix_add;ki=$ki;mutposorig=$mutposorig2[$ki];mutposcurr=$mutposcurr2[$ki];coor=chr12\_$begind2\_$endind2\_$name,$igtypeprint\_0\_$strand2\n";
						}
					}

					my $constype1 = $cons{$igtypeprint}{1}{$begind2print+$fix_add}{type1}; $constype1 = "UNDEF_constype1" if not defined $constype1;
					my $constype2 = $cons{$igtypeprint}{1}{$begind2print+$fix_add}{type2}; $constype2 = "UNDEF_constype2" if not defined $constype2;
					my ($cons_fix, $printzfix2);
					$printzfix2 = $printzfix->{$begind2print};
					($cons_fix, $constype2, $printzfix2) = check_cons_fix($type1, $type2, $constype1, $constype2, $printzfix2, \%cons, $begind2print, $igtypeprint, $strand2, $fix_add);
					my ($change, $ctype1, $ctype2) = ("", $type1, $type2);
					if ($cons_fix =~ /CHANGE/ and $cons_fix !~ /DONT/) {
						($change, $ctype1, $ctype2) = $cons_fix =~ /^(CHANGE) ([A-Za-z0-9]+)_(.+)$/;
						$printzfix2 .= "\t$LRD ctype1 or 2 undefined! cons_fix=$cons_fix\n" if not defined $ctype1 or not defined $ctype2;
						$ctype1 = $type1 if not defined $ctype1;
						$ctype2 = $type2 if not defined $ctype2;
						$printzfix2 .= "\t$YW>Changing $LGN$type1 $type2$N into $LPR$ctype1 $ctype2$N\n" if defined $ctype1 and defined $ctype2;
						$typehash1->{$type1}{$type2} --;
						$typehash1->{$ctype1}{$ctype2} ++;
						if ($typehash1->{$type1}{$type2} == 0) {
							undef $typehash1->{$type1}{$type2};
							delete $typehash1->{$type1}{$type2};
						}
					}
					my $origcons_fix = $cons_fix; $origcons_fix =~ s/[ ]+/_/g;

					$cons_fix = $cons_fix =~ /GOOD/ ? $LGN . $cons_fix . $N : $cons_fix =~ /BAD/ ? $LRD . $cons_fix . $N : $LCY. $cons_fix . $N;
					$printzfix->{$begind2print} = $printzfix2;


#foreach my $type1 (sort keys %type1) {
#	next if $type1 =~ /^(inspos|delpos|matpos|mispos|mhpos|tot|nuc|R2)$/;
#	foreach my $type2 (sort keys %{$type1{$type1}}) {
#		my $mutpos_joined = $lochash1->{$type1}{$type2}{join}; $mutpos_joined = "N/A" if not defined $mutpos_joined;
#		my $number = $type1{$type1}{$type2};


					$printzfix->{$begind2print} .= "AJ851868.3\t$begind2print\t$endind2print\t$igtypeprint;$type1;$type2;$ind\t0\t$strand2\t$igtypeprint\t$constype1,$constype2\t$cons_fix\tfix_add=$fix_add\n";
#					$printzfix->{$begind2print} .= "\tIgG1=$IgG1beg - $IgG1end, begorig2=$begorig2 endorig2=$endorig2, begind2=$begind2 = begorig2=$begorig2 + (ind=$ind - juncposfix=$juncposfix + mhlen=$mhlen) + 1\n";
					print $outloc "AJ851868.3_$igtypeprint\t$begind2print\t$endind2print\t$ctype1\_$ctype2\t0\t$strand2\tmutorig=$type1\_$type2;change=$origcons_fix;constype=$constype1\_$constype2;name=$name;type=PREY;ig=$igtypeprint;fix_add=$fix_add;ki=$ki;mutposorig=$mutposorig2[$ki];mutposcurr=$mutposcurr2[$ki];coor=chr12\_$begind2\_$endind2\_$name,$igtypeprint\_0\_$strand2\n";
				}
				$lastmutpos = $ind if $lastmutpos < $ind;
				if (defined $temp[$ind]) {
					$temp{$ind} = $temp[$ind] if not defined $temp{$ind};
					$temp{$ind} .= ",$nuc";
				}
				else {
					$temp[$ind] = $nuc;
				}
			}
		}
	}

	print "\n\nBEFORE:\n$printhash1\nAFTER:\n";
	# dump typehash1
	foreach my $type1 (sort keys %{$typehash1}) {
		foreach my $type2 (sort keys %{$typehash1->{$type1}}) {
			if ($typehash1->{$type1}{$type2} =~ /HASH/) {
				my $typehash1temp = $typehash1->{$type1}{$type2};
				foreach my $key (sort keys %{$typehash1temp}) {
					print "type1=$LCY$type1$N, type2=$LGN$type2$N, key=$key value=$LPR$typehash1->{$type1}{$type2}{$key}$N\n";
				}
			}
			else {
				print "type1=$LCY$type1$N, type2=$LGN$type2$N, value=$LPR$typehash1->{$type1}{$type2}$N\n";
			}
		}
	}
	print "\n\n";
	
	#	LOG($outBigLog, "AAAAAAA\n");
	open (my $outztemp, ">", $outztempFile) or die;
	open (my $outzfix, ">", $outzfixFile) or die;
	foreach my $begindprint (sort {$a <=> $b} keys %{$printzfix}) {
		print $outzfix "$printzfix->{$begindprint}";
		LOG($outBigLog, "$printzfix->{$begindprint}");
	}
	#	LOG($outBigLog, "BBBBBBB\n");

	print $outztemp "$chr1\t$coorz{begind1}\t$coorz{endind1}\t$name\t0\t$strand1\n" if defined $coorz{begind1};
	print $outztemp "$chr2\t$coorz{begind2}\t$coorz{endind2}\t$name\t0\t$strand2\n" if defined $coorz{begind2};

#my $outzfixFile = $outFile . ".outzfix.bed.temp";
#my $outztempFile = $outFile . ".outztemp.bed.temp";
#my $outzfaFile = $outFile . ".outztemp.fa";

	system("fastaFromBed -fi /home/mitochi/Bowtie2_indexes/mm9/mm9.fa -bed $outztempFile -fo $outzfaFile -s -name");
	my @fa = `cat $outzfaFile`;
	for (my $i = 0; $i < @fa; $i++) {
		chomp($fa[$i]);
		next if $fa[$i] =~ /^>/;
		my $space0 = defined $coorz{ind1} ? join("", (" ") x $coorz{ind1}) : "";
		my $space1 = defined $coorz{ind2} ? join("", (" ") x $coorz{ind2}) : "";
		LOG($outBigLog, "\n\n" . $space0 . colorize($fa[$i]) . "\n") if $i == 1;
		LOG($outBigLog, $space1 . colorize($fa[$i]) . "\n\n") if $i == 3;
	}
	
	
	for (my $i = 0; $i < $lastmutpos+1; $i++) {
		$temp[$i] = " " if not defined $temp[$i];
	}
	for (my $i = 0; $i < @temp; $i++) {
		$temp[$i] = " " if not defined $temp[$i];
	#		die "lastmutpos = $lastmutpos+1, i=$i temp=undefined!\n" if not defined $temp[$i];
	}
	
	LOG($outBigLog, "\n" . colorize($DATA->{$name}{'seqCON_aln'}) . "\n" . colorize($seqCON_orig) . "\n");
	LOG($outBigLog, colorize(join("", @temp)) . "\n\n") if defined $coorz{begind1};
	foreach my $i (sort {$a <=> $b} keys %temp) {
		LOG($outBigLog, "$i. $temp{$i}\n\n");
	}
	
	LOG($outBigLog, "\nSTRAND1=$strand1, STRAND2=$strand2\nOriginal:\n");
	LOG($outBigLog, "BAIT: $DATA->{$name}{origbeg1}-$DATA->{$name}{origend1}, junc=$DATA->{$name}{origbegJ1}, T=$DATA->{$name}{origbegT1}-$DATA->{$name}{origendT1}\n");
	LOG($outBigLog, "PREY: $DATA->{$name}{origbeg2}-$DATA->{$name}{origend2}, junc=$DATA->{$name}{origbegJ2}, T=$DATA->{$name}{origbegT2}-$DATA->{$name}{origendT2}\n");
	
	LOG($outBigLog, "\nOriginal2:\n");
	LOG($outBigLog, "BAIT: $DATA->{$name}{begorig1}-$DATA->{$name}{endorig1}, junc=$DATA->{$name}{begJ1}, T=$DATA->{$name}{begT1}-$DATA->{$name}{endT1}\n");
	LOG($outBigLog, "PREY: $DATA->{$name}{begorig2}-$DATA->{$name}{endorig2}, junc=$DATA->{$name}{begJ2}, T=$DATA->{$name}{begT2}-$DATA->{$name}{endT2}\n");
	
	LOG($outBigLog, "\nCurrent:\n");
	LOG($outBigLog, "BAIT: $DATA->{$name}{beg1}-$DATA->{$name}{end1}, junc=$DATA->{$name}{begJ1}, T=$DATA->{$name}{begT1}-$DATA->{$name}{endT1}\n");
	LOG($outBigLog, "PREY: $DATA->{$name}{beg2}-$DATA->{$name}{end2}, junc=$DATA->{$name}{begJ2}, T=$DATA->{$name}{begT2}-$DATA->{$name}{endT2}\n");
	
	
	($igtypez, $typez1, $typez2, $all) = count_final("R1", $igtypez, $typez1, $typez2, $muthash, $typehash1, $lochash1, $DATA->{$name}, $name, 1, $outBigLog, $outLog, $namewant);
	
	LOG($outBigLog, "\n---------$LGN DONE with $LCY$name$N ---------\n\n\n");
	LOG($outLog, "\n---------$LGN DONE with $LCY$name$N ---------\n\n\n","NA");
	LOG($outBigLog, "${N}\n",$NA);
	LOG($outLog, "${N}\n","NA");
	
	#	LOG($outBigLog, "${YW}------------------------------------------------\n6b. Tabulating mutation from seqQ2 (Read pair #2)\n------------------------------------------------${N}\n" . date() . "\n",$NA);
	#	LOG($outLog, "${YW}------------------------------------------------\n6b. Tabulating mutation from seqQ2 (Read pair #2)\n------------------------------------------------${N}\n" . date() . "\n","NA");
	#
	#	foreach my $mutpos (sort keys %{$muthash}) {
	#		my $mut1 = $muthash->{$mutpos}{1}; $mut1 = "none" if not defined $mut1;
	#		my $mut2 = $muthash->{$mutpos}{2}; $mut2 = "none" if not defined $mut2;
	#		if ($mutpos <= $beg_fix or $mutpos >= $end_fix and (($mut1 eq "none" and $mut2 ne "none") or ($mut2 eq "none" and $mut1 ne "none"))) {
	#			if ($mut1 eq "none" and $mut2 ne "none") {
	#				my $mut = $mut2;
	#				my ($type1, $type2) = $mut =~ /^(ins|del|mat|mis|mh)_(.+)$/;
	#				DIELOG($outBigLog, "Can't find type1 from mutpos=$mutpos mut=$mut; not ins/del/mat/mis/mh_?\n\n") if not defined $type1;
	#				$typehash1->{$type1}{$type2} ++;
	#				$typehash1->{tot}{tot} ++;
	#			}
	#		}
	#		elsif ($mut1 ne $mut2 and $mut1 ne "none" and $mut2 ne "none") {
	#			my $color = $mut1 eq $mut2 ? "" : $mutpos >= $end_fix ? "" : $mutpos <= $beg_fix ? "" : "$LRD";
	#			LOG($outBigLog, "$color" . "mutpos${N}=${LGN}$mutpos${N} mut1=${LCY}$mut1${N} mut2=${LGN}$mut2${N}\n");
	#			my $mut = $mut2;
	#			my ($type1, $type2) = $mut =~ /^(ins|del|mat|mis|mh)_(.+)$/;
	#			DIE($outBigLog, "Can't find type1 from mutpos=$mutpos mut=$mut; not ins/del/mat/mis/mh_?\n\n") if not defined $type1;
	#			$typehash1->{R2}{$type1}{$type2} ++;
	#		}
	#		$typehash1->{tot}{tot} = 100 if $typehash1->{tot}{tot} < 10;
	#	}
	#
	#	($igtypez, $typez1, $typez2, $all) = count_final("R1", $igtypez, $typez1, $typez2, $muthash, $typehash1, $DATA->{$name}, $name, 1, $outBigLog, $outLog, $namewant);
	
}
close $out1;
close $outloc;

LOG($outBigLog, "${YW}
------------------------------------------------
7. STATISTICS 
------------------------------------------------${N}\n" . date() . "\n

");
my $okay = (keys %okay);
my $totalpair = (keys %totalpair);
LOG($outBigLog, "\ntotalpair\t$totalpair\npairokay (at least 1 pair good)\t$okay\ntotalall\t$totalall\ngood\t$good\nnexted_sequndef\t$nexted_sequndef\nnexted_bothshort\t$nexted_bothshort\nnexted_begshort\t$nexted_begshort\nnexted_endshort\t$nexted_endshort\nnexted_juncmissing\t$nexted_juncmissing\nnexted_others (e.g. not IgM/G1/G3)\t$nexted_others\n");
LOG($outLog, "\ntotalpair\t$totalpair\npairokay (at least 1 pair good)\t$okay\ntotalall\t$totalall\ngood\t$good\nnexted_sequndef\t$nexted_sequndef\nnexted_bothshort\t$nexted_bothshort\nnexted_begshort\t$nexted_begshort\nnexted_endshort\t$nexted_endshort\nnexted_juncmissing\t$nexted_juncmissing\nnexted_others (e.g. not IgM/G1/G3)\t$nexted_others\n", "NA");
#close $outLog;
#close $outBed;

open (my $outALL, ">", "$outFile.all") or DIELOG($outLog, "Cannot write to $outFile.all: $!\n");
print $outALL "#sample\tigtype\ttype1\ttotal_read\taverage_length";
foreach my $typez2 (sort keys %{$typez2}) {
	print $outALL "\t$typez2";
}
print $outALL "\n";
foreach my $igtypez (sort keys %{$igtypez}) {
	my $total_read = $all->{$igtypez}{total_read};
	my $total_nuc = int(100*$all->{$igtypez}{total_nuc} / $total_read +0.5)/100;
	foreach my $typez1 (sort keys %{$typez1}) {
		print $outALL "$sampleID\t$igtypez\t$typez1\t$total_read\t$total_nuc";
		foreach my $typez2 (sort keys %{$typez2}) {
			my $value = $all->{$igtypez}{$typez1}{$typez2}; $value = 0 if not defined $value;
#			$value = $value / $total_read;
#			$value *= 100 if $typez2 =~ /^perc/;
#			if ($value eq 0) {
#3			}
#			elsif (abs($value) < 1) {
#				my ($value2) = $value =~ /^(\-?0\.0*[1-9][0-9]?)[0-9]*$/;
#				$value = $value2 if defined $value2;
#			}
#			else {
#				my $flag = $value < 0 ? -1 : 1;
#				$value = int($value*10+($flag*0.5))/10;
#			}
			print $outALL "\t$value";
		}
		print $outALL "\n";
	}
}
close $outALL;

LOG($outBigLog, "\n----------------------------\n");
LOG($outBigLog, "Example from outFile.all $LCY$outFile.all$N:\n");
LOG($outBigLog, "----------------------------\n");
system("head -n 5 $outFile.all");
LOG($outBigLog, "\n----------------------------\n\n");

sub check_cons_fix {
	my ($type1, $type2, $cons1, $cons2, $printzfix2, $cons, $begindprint, $igtypeprint, $strand, $fix_add) = @_;
	if ($strand eq "+") {
		my $newcons2 = "";
		my @cons2 = split("", $cons2);
		for (my $i = 0; $i < @cons2; $i++) {
			my $cons2temp = $cons2[$i];
			$cons2temp =~ tr/_ACTGactg/_TGACtgac/;
			$newcons2 .= $cons2temp;
		}
		$cons2 = $newcons2;
	}
	my $origtype1 = $type1;
	my $origtype2 = $type2;
	my $origcons1 = $cons1;
	my $origcons2 = $cons2;

	if ($type1 eq "mat") {
		($type2) = $type2 =~ /^\w+_(\w+)$/;
	}
	if ($cons1 eq "mat") {
		($cons2) = $cons2 =~ /^\w+_(\w+)$/;
	}

#	$printzfix2 .= "\ttype=$type1, $type2, cons=$cons1, $cons2" if $strand eq "+";
	if ($cons1 =~ /UNDEF/ or $cons2 =~ /UNDEF/) {
		return("BAD cons1 UNDEF cons2 UNDEF", $origcons2, $printzfix2);
	}
	if ($type1 eq "mat" and $cons1 =~ /^(mis|mat)$/) {
		my ($exp, $obs) = $cons2 =~ /^(\w+)_(\w+)$/ if $cons1 ne "mat";
		$obs = $cons2 if $cons1 eq "mat";
#		$printzfix2 .= ", mis/mat type1=$type1, type2=$type2, cons1=$cons1, cons2=$cons2, obs=$obs (from cons2)\n" if $strand eq "+";
		return("GOOD $type1", $origcons2, $printzfix2) if $cons2 =~ /$type2/i;
		return("DONT CHANGE mis_$obs\_$origtype2", $origcons2, $printzfix2);
	}
	elsif ($type1 eq "mat" and $cons1 =~ /^del$/) {
		my ($obs) = $cons2 =~ /^(\w)\w*$/;
#		$printzfix2 .= ", del type1=$type1, type2=$type2, cons1=$cons1, cons2=$cons2, obs=$obs (from cons2)\n" if $strand eq "+";
		return("GOOD $type1", $origcons2, $printzfix2) if $obs eq $type2;
		return("DONT CHANGE mis_$type2\_$obs", $origcons2, $printzfix2);
	}
	elsif ($type1 eq "del" and $cons1 =~ /^del$/) {
#		my ($obs) = $cons2 =~ /^(\w)\w*$/;
#		$printzfix2 .= ", del type1=$type1, type2=$type2, cons1=$cons1, cons2=$cons2\n" if $strand eq "+";
		return("CHANGE mat_$type1\_$type1", $origcons2, $printzfix2) if $type2 eq $cons2;
		return("GOOD $type1", $origcons2, $printzfix2);
	}
	elsif ($type1 eq "mis" and $cons1 =~ /^(mis|mat)$/) {
		my ($exp, $obs) = $type2 =~ /^(\w+)_(\w+)$/;
#		$printzfix2 .= ", mis/mat type1=$type1, type2=$type2, cons1=$cons1, cons2=$cons2, exp=$exp, obs=$obs (from type2)\n" if $strand eq "+";
		return("GOOD $type1", $origcons2, $printzfix2) if $cons2 !~ /$obs/i;
		return("CHANGE mat_$obs\_$obs", $origcons2, $printzfix2);
	}
	elsif ($type1 eq "mis" and $cons1 =~ /^del$/) {
		my ($exp, $obs) = $type2 =~ /^(\w)\w*_(\w)\w*$/;
#		$printzfix2 .= ", del type1=$type1, type2=$type2, cons1=$cons1, cons2=$cons2, exp=$exp, obs=$obs (from type2)\n" if $strand eq "+";
		return("CHANGE mat_$type2\_$obs", $origcons2, $printzfix2) if $obs eq $type2;
		return("GOOD $type1", $origcons2, $printzfix2);
	}
	elsif ($type1 eq "mh") {
		my $obs = $type2;
		my $obslen = length($obs);
		my $good = 0;
		my $ind = 0;
		my $newcons = "";
		for (my $i = $begindprint; $i < $begindprint + $obslen; $i++) {
			my $cons3 = $cons{$igtypeprint}{1}{$i+$fix_add}{type1}; $cons3 = "UNDEF" if not defined $cons3;
			my $cons4 = $cons{$igtypeprint}{1}{$i+$fix_add}{type2}; $cons4 = "UNDEF" if not defined $cons4;
			my $obs2 = substr($type2, $ind, 1);
			if ($cons3 =~ /^(mis|mat)$/) {
				$good ++ if $cons4 =~ /$obs2/;
				if ($cons3 eq "mat") {
					$cons3 = "";
					($cons4) =~ s/^(\w)_\w$/$1/;
					if ($newcons eq "") {
						$newcons .= "i=$cons4";
					}
					else {
						$newcons .= "$cons4";
					}
				}
				else {
					if ($newcons eq "") {
						$newcons .= "$i=$cons3\_$cons4,";
					}
					else {
						$newcons .= ",$cons3\_$cons4,";
					}
				}
			}
			$ind ++;
		}
		$newcons =~ s/,,/,/g;
		$newcons =~ s/,$//;
#		$printzfix2 .= ", type1=$type1, type2=$type2, newcons=$newcons, good=$good/$obslen\n" if $strand eq "+";
		if ($good == $obslen or ($obslen > 4 and $good >= $obslen - 2) or ($obslen > 10 and $good / $obslen > 0.9)) {
			return("GOOD good=$good obslen=$obslen $type1", $newcons, $printzfix2);
		}
		else {
			return("BAD $type1 good=$good obslen=$obslen $type1", $newcons, $printzfix2);
		}
	}
	elsif ($type1 eq "ins") {
#		$printzfix2 .= ", ins type1=$type1, type2=$type2, cons1=$cons1, cons2=$cons2\n" if $strand eq "+";

		if ($type1 eq $cons1 and $type2 eq $cons2) {
			return("CHANGE mat_$type2", $origcons2, $printzfix2);
		}
		else {
			return("GOOD $type1", $origcons2, $printzfix2);
		}
	}

#	$printzfix2 .= ", unknown\n" if $strand eq "+";
	return("UNKNOWN", $origcons2, $printzfix2);
}


sub muscle {
	my ($cmd) = @_;
	return(`echo '$cmd' | muscle $muscleparam 2> /dev/null`);

}

sub header_print {
	my ($to_print, $extra, $extra_below) = @_;
	$to_print = "" if not defined $to_print;
	$extra = "" if not defined $extra;
	$extra_below = "" if not defined $extra_below;
	my $date = date(); $date =~ s/\:\s*$//g;
	my $header_print = "";
	   $header_print .= "\n";
		$header_print .= "${YW}------------------------------------------------${N}\n";
		$header_print .= "${YW}$to_print${N}\n";
		$header_print .= "$extra\n" if $extra ne "";
		$header_print .= "${YW}------------------------------------------------${N}\n";
		$header_print .= "$date\n";
		$header_print .= "$extra_below" if $extra_below ne "";
		$header_print .= "\n";
	return($header_print);
}


sub parse_metaFile {
	my ($metaFile, $sampleIDwant) = @_;

	LOG($outBigLog, "${LGN}>Example metaFile:$N\n");


	my $metaHASH;

	my $linecount = 0;
	my $linefoundcount = 0;
	my $linefound = "";
	my $log = "";
	my $die = 0;
	open (my $in0, "<", $metaFile) or $die = 1;

	if ($die eq 1) {
		$log .= "  Failed to read from $metaFile: $!\n";
		return("", "", $log, $die);
	}

	while (my $line = <$in0>) {
		chomp($line);
		$linecount ++;
		if ($line =~ /(^Library Assembly|\tStart\t|\tEnd\t|\tPrimer\t)/) {
			$log .= "\n header: $LGN$line$N\n";
		}
		else {
			my ($sampleID, $genome, $chr, $beg, $end, $strand, $barcode, $primer, $adapter, $desc) = split("\t", $line);
			$log .= " row#$linecount : $LCY$line$N\n" if $linecount < 3;
			$log .= "         -> barcode=" . colorize($barcode) . ", primer=" . colorize($primer) . ", adapter=" . colorize($adapter) . "\n" if $linecount < 3;
			$metaHASH->{$sampleID}{coor}    = $line;
			$metaHASH->{$sampleID}{barcode} = $barcode;
			$metaHASH->{$sampleID}{primer}  = $primer;
			$metaHASH->{$sampleID}{adapter} = $adapter;
			$linefoundcount = $linecount if $sampleID eq $sampleIDwant;
			$linefound = $line if $sampleID eq $sampleIDwant;
			#Library Assembly        Chr     Start   End     Strand  MID     Primer  Adapter Description
			#W1      mm9     chr12   114664888       114664910       -       AGTCAA  CACACAAAGACTCTGGACCTC   CCACGCGTGCTCTACA        Smu_as_bait_site
		}
	}
	close $in0;
	if ((keys %{$metaHASH}) == 0) {
		$log .= "\n${LRD}ERROR!$N There is no sampleId in meta file ${LCY}$metaFile${N} as keys \$metaHASH is 0!\n";
		$die = 1;
	}
	else {
		$metaHASH->{$sampleIDwant}{bcpri} = $metaHASH->{$sampleIDwant}{barcode} . $metaHASH->{$sampleIDwant}{primer};
		$log .= "\n ----------\n ${LGN}Done parsing metaFile$N\n\n Parsed " . (keys %{$metaHASH}) . "$N sample IDs and found sampleIDwant=$LCY$sampleIDwant$N at line #$LGN$linefoundcount$N:\n $LCY$linefound\n";
		$log .= " BARCODE = " . colorize($metaHASH->{$sampleIDwant}{barcode}) . "\n";
		$log .= " PRIMER  = " . colorize($metaHASH->{$sampleIDwant}{primer}) . "\n";
		$log .= " BC+PRI  = " . colorize($metaHASH->{$sampleIDwant}{barcode}) . "" . colorize($metaHASH->{$sampleIDwant}{primer}) . "\n";
		$log .= " ADAPTER = " . colorize($metaHASH->{$sampleIDwant}{adapter}) . "\n";
	}
	return($metaHASH->{$sampleIDwant}{bcpri}, $metaHASH->{$sampleIDwant}{adapter}, $log, $die);
}

sub parse_inputFile {
	my ($inputFile, $outBedFile, $namewant) = @_;
	my ($log, $die) = ("", 0);
	my $linecount = 0;
	my ($data, @def);

	open (my $in1, "<", $inputFile) or return($data, $log . "\nCannot read from $LCY$inputFile$N: $!\n", 1);
	open (my $outBed, ">", "$outBedFile") or return($data, $log . "\nCannot write into $LCY$outBedFile$N: $!\n", 1);

	my ($total_line_inputFile) = `wc -l $inputFile` =~ /^(\d+)/;

	$log .= " -> Parsing $LGN$total_line_inputFile$N lines from inputFile=$LCY$inputFile$N\n";
	while (my $line = <$in1>) {
		chomp($line);
		next if $line =~ /^#/;
		$log .= "  - parsed ${LGN}$linecount${N} / ${LCY}$total_line_inputFile${N} on "  . date() . "\n" if $linecount % 1e3 == 0;
		my @arr = split("\t", $line);
		if ($line =~ /^Qname\t/) {
			@def = @arr;
			$linecount ++;
			next;
		}
		$linecount ++;
		my ($name, $juncid, $chr2, $junc2, $strand2, $beg2, $end2, $chr1, $beg1, $end1, $strand1, $begT1, $endT1, $begT2, $endT2, $lenTN, $cigarQ1, $cigarQ2, $seqTN, $defseqjunc, $defbarcode, $unaligned, $baitonly, $uncut, $misprimed, $freqcut, $largegap, $mapqual, $breaksite, $sequential, $repeatseq, $duplicate) = @arr;
		my ($origbeg1) = $beg1;
		my ($origend1) = $end1;
		my ($origbeg2) = $beg2;
		my ($origend2) = $end2;

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
				if ($i != 0 and $i % 5 == 0 or ($i == @def - 1)) {
					$data->{$name}{line} .= "\n";
				}
				else {
					$data->{$name}{line} .= ",";
				}
			#  if ($i != 0 and $i % 5 != 0 and $i != @def - 1);
			}
		}
		$data->{$name}{line} .= $myprintcigarseq . "\n$dashhead\n";
		# TCseq uses both 1 based for coordinate. fastaFromBed uses 0/1 (1 beg 0 end). So turn ALL TCseq end coordinate into 1 based (beg-1)
		# starts from 1, 
		$data->{$name}{origbeg1} = $beg1 - 1;
		$data->{$name}{origend1} = $end1;
		$data->{$name}{origbeg2} = $beg2 - 1;
		$data->{$name}{origend2} = $end2;
		$data->{$name}{origbegT1} = $begT1 - 1;
		$data->{$name}{origendT1} = $endT1;
		$data->{$name}{origbegT2} = $begT2 - 1;
		$data->{$name}{origendT2} = $endT2;
		$data->{$name}{origbegJ1} = $junc1 - 1;
		$data->{$name}{origbegJ2} = $junc2 - 1;
		$data->{$name}{origendJ1} = $junc1;
		$data->{$name}{origendJ2} = $junc2;
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


sub parse_outBedFile {
	my ($data, $outBedFile, $inputFA, $namewant) = @_;
	my ($log, $die) = ("", 0);
	my $cmd = "fastaFromBed -fi /home/mitochi/Bowtie2_indexes/mm9/mm9.fa -bed $outBedFile -fo $inputFA -s -name";
	$log .= date() . "$cmd\n";
	system("$cmd") == 0 or return($data, $log . "\nFailed to run $LPR$cmd${N}: $!\n\n", 1);
	
	$log .= date() . "Parsing inputFA1 ${LCY}$inputFA${N}\n";
	open (my $inFA1, "<", $inputFA) or return($data, $log . "\nCannot read from $inputFA: $!\n", 1);
	my $fasta = new FAlite($inFA1);
	my $lastdef = "INIT"; my $lastseqlong = "INIT";
	my $muhentry = "";
	while (my $entry = $fasta->nextEntry()) {
		my $seqlong = $entry->seq;
		$seqlong =~ tr/acgtn/ACGTN/;
		my $seq = $seqlong;
	   my $def = $entry->def;
		my ($name) = $def;
		$name =~ s/^>//;
		$name =~ s/\([\+\-]\)$//;
		$def = $name;
		
	
		if (not defined $name) {
			$log .= date() . "Can't parse name from fasta def $LPR$def${N}\n";
			return($data, $log, 1);
		}

		$muhentry .= "     " . $entry->def . "\n                                     " . colorize($entry->seq) . "\n" if $name eq $namewant;
		if ($lastdef eq $def) {
			my ($seqC1, $seqC2) = ($lastseqlong, $seqlong);
			($data->{$name}{seqC1L}, $data->{$name}{seqC1R}) = $seqC1 =~ /^(.{$data->{$name}{lenC1L}})(.{$data->{$name}{lenC1R}})$/; #seq ONE extra is as long as seq2 chunk (len2)
			($data->{$name}{seqC2L}, $data->{$name}{seqC2R}) = $seqC2 =~ /^(.{$data->{$name}{lenC2L}})(.{$data->{$name}{lenC2R}})$/; #seq ONE extra is as long as seq2 chunk (len2)
			if ($data->{$name}{gap} ne 0) {
#				$data->{$name}{seqC1R} = join("", ("N") x $data->{$name}{gap}) . $data->{$name}{seqC1R};
#				$data->{$name}{seqC2R} = join("", ("N") x $data->{$name}{gap}) . $data->{$name}{seqC2R};
#				$data->{$name}{lenC1R} += $data->{$name}{gap};
#				$data->{$name}{lenC2R} += $data->{$name}{gap};
			}
			$data->{$name}{seqCON} = $data->{$name}{seqC1L} . $data->{$name}{seqC2R};
			$data->{$name}{lenCON} = length($data->{$name}{seqCON});
			$data->{$name}{mhseqlong} = $data->{$name}{mhbeg} > 0 ? join("", (".") x $data->{$name}{mhbeg}) . $data->{$name}{mhseq} . join("", (".") x (length($data->{$name}{seqTN}) - length($data->{$name}{mhseq}) - $data->{$name}{mhbeg})) : "";
			$log .= "\n\n -> ${LGN}REFERENCE SEQUENCE EXAMPLE:${N}
  - READ      = ${LGN}$name${N} (used in hash)
  - ORDER     = $data->{$name}{order}
  - DEF1      = $lastdef (lastdef)
  - DEF2      = $def (currdef)
  - FULL LEN1 = ${LGN}$data->{$name}{len1}${N} bp = orig/colored 1_con  , left of junc (${LGN}$data->{$name}{lenC1L}${N} bp) + extra added to 1_con, right of junc (${LGN}$data->{$name}{lenC1R}${N} bp)
  - FULL LEN2 = ${LGN}$data->{$name}{len2}${N} bp = extra added to 2_con, left of junc (${LGN}$data->{$name}{lenC2L}${N} bp) + orig/colored 2_con  , right of junc (${LGN}$data->{$name}{lenC2R}${N} bp)
  - CONSENSUS = ${LGN}$data->{$name}{lenCON}${N} bp = REFERENCE sequence (except for microhomology, which uses the extra sequences to compare)
  - MH        = ${LGN}$data->{$name}{mhbeg}-$data->{$name}{mhend}$N ($LGN$data->{$name}{mhpos}$N bp)
  - GAP       = ${LGN}$data->{$name}{gap}$N

  - ${LGN}RAW REFERENCE SEQUENCES (top=IgM BAIT, bottom=FAR PREY, 0/1 based):${N}
$muhentry

  - ${LGN}FULL (should be identical to RAW REFERENCE SEQUENCES):${N}
    - 1_con/IGM BAIT FULL (RAW SEQ) = " . colorize($data->{$name}{seqC1L} . $data->{$name}{seqC1R}) . "
    - 2_con/FAR PREY FULL (RAW SEQ) = " . colorize($data->{$name}{seqC2L} . $data->{$name}{seqC2R}) . "

  - ${LGN}BROKEN INTO BAIT/PREY and their EXTRA${N}
    - 1_con/IGM BAIT + extra1       = " . colorize($data->{$name}{seqC1L}) . $data->{$name}{seqC1R} . "
    - extra2 + 2_con/FAR PREY       = " . $data->{$name}{seqC2L} . colorize($data->{$name}{seqC2R}) . "

  - ${LGN}consensus (combining seq1 and seq2 without extra):${N}
    - 3_con/CONSENSUS               = " . colorize($data->{$name}{seqC1L} . $data->{$name}{seqC2R}) . "
	
  - ${LGN}seqTLX and MH${N}:
    - seqTLX                        = " . colorize($data->{$name}{seqTN}) . "
    - MH                            = " . colorize($data->{$name}{mhseqlong}) . "

" if $namewant eq $name;


			$lastdef = "INIT";
			$lastseqlong = "INIT";
			$muhentry = "";
		}
		else {
			$lastdef = $def;
			$lastseqlong = $seqlong;
		}
	}
	close $inFA1;
	return($data, $log, $die);
}

sub parse_fastqFile {
	my ($data, $fastqFile, $readorient, $total_read, $namewant) = @_;
	my $seqName = $readorient eq "R1" ? "seqQ1" : "seqQ2";
	my $linecount = 0;
	my $count;
	my ($log, $die) = ("", 0);
	my $nexted = 0;
	open (my $inFQ1, "zcat $fastqFile|") or return($data, $log . "\nCannot read from $fastqFile: $!\n", 1);
	$log .= " -> Parsing Fastq of $LGN$readorient$N: $LCY$fastqFile$N\n";
	while (my $line = <$inFQ1>) {
		chomp($line);
		$log .= "   - linecount=${LGN}$linecount${N}, processed ${LGN}" . (keys %{$count}) . "${N}/${LCY}$total_read${N}\n" if $linecount % 1e4 == 0;
		$linecount ++;
		my ($name) = $line;
		$name =~ s/^\@//;
		$name =~ s/ .*$//;
		$line = <$inFQ1>; chomp($line);
		my $seq = $line;
		$line = <$inFQ1>; chomp($line);
		$line = <$inFQ1>; chomp($line);
		my $qua = $line;
		
		if (not defined $data->{$name}) {
			$log .= "     - ${LRD}$name$N nexted as doesn't exist data HASH\n" if $nexted < 5;
			$log .= "     - ...<more than 5 reads were nexted>\n" if $nexted == 5;
			$nexted ++;
			next;# if $nexted >= 5;
		}
		$count->{$name} ++;
		$seq =~ tr/acgtn/ACGTN/;
		$seq = mitochy::revcomp($seq) if $readorient eq "R2";
		$data->{$name}{$seqName} = $seq;
		if ($name eq $namewant) {
			$log .= " ->$LGN Example READ $readorient Fastq:\n";
			$log .= "   - name = $LGN$name$N\n";
			$log .= "   - seq  = $LGN$seq$N\n";

			my $defz1 = "seq$readorient";
			my $defz2 = "seqTN";
			my @res = muscle(">$defz1\n$seq\n>$defz2\n$data->{$name}{seqTN}\n");
#			my @res = `$muscleCMD`;
#			my ($res);
			my ($res) = parse_fasta_simple(undef, \@res);
			$log .= "\n";
			$log .= "$defz1\t$res->{$defz1}\n";
			$log .= "$defz2\t$res->{$defz2}\n\n";
		}
		$log .= "   - SEQ$readorient LEN = " . length($seq) . "\n" if $nexted < 5;
		last if (keys %{$count}) >= scalar(keys %{$data});
	}
	close $inFQ1;
	$log .= " -> ${LGN}Done parsing fastq $readorient$N, read processed=${LGN}" . (keys %{$count}) . "${N}/$LGN$linecount$N, read needed = ${LCY}$total_read${N}, total nexted=$LGN$nexted$N)\n\n";
	return ($data, $log, $die);
}


sub parse_fasta_simple {
	my ($data, $res) = @_;
	my $prevdef = "INIT";
	my $def = "";
#	print "\nPARSE FASTA SIMPLE--------------\n";
	for (my $i = 0; $i < @{$res}; $i++) {
		chomp($res->[$i]);
		if ($res->[$i] =~ /^>/) {
			$res->[$i] =~ s/^>//;
			$def = $res->[$i];
			$data->{$def} = "";
			if ($prevdef ne "INIT") {
#				print "$prevdef\t" . colorize($data->{$prevdef}) . "\n";
			}
			$prevdef = $def;
		}
		else {
			$data->{$def} .= $res->[$i];
		}
	}
	if ($prevdef ne "INIT") {
#		print "$prevdef\t" . colorize($data->{$prevdef}) . "\n";
	}
#	print "\nPARSE FASTA SIMPLE-------------- DONE!\n";
	return($data);

}

sub count_final {
	my ($readorient1, $igtypez, $typez1, $typez2, $muthash, $typehash1, $lochash1, $data, $name, $printout, $outBigLog, $outLog, $namewant) = @_;
	my $NA = $name eq $namewant ? undef : "NA";
	my %type1 = %{$typehash1};
	$type1{tot}{tot} = 0;
	foreach my $type1 (sort keys %type1) {
		next if $type1 =~ /^(inspos|delpos|matpos|mispos|mhpos|nuc|tot|R2)$/;
		foreach my $type2 (sort keys %{$type1{$type1}}) {
			my $number = $type1{$type1}{$type2};
			LOG($outBigLog, "$type1\_$type2,$number (total=$type1{tot}{tot}+$number=", $NA);
			LOG($outLog, "$type1\_$type2,$number (total=$type1{tot}{tot}+$number=", "NA");
			$type1{tot}{tot} += $number;
			LOG($outBigLog, "$type1{tot}{tot})\n", $NA);
			LOG($outLog, "$type1{tot}{tot})\n", "NA");
		}
	}

	my %cur;
	LOG($outBigLog, "tot,$type1{tot}{tot}\n\n",$NA);
	LOG($outLog, "tot,$type1{tot}{tot}\n\n", "NA");
	foreach my $type1 (sort keys %type1) {
		next if $type1 =~ /^(inspos|delpos|matpos|mispos|mhpos|tot|nuc|R2)$/;
		foreach my $type2 (sort keys %{$type1{$type1}}) {
			my $mutpos_joined = $lochash1->{$type1}{$type2}{join}; $mutpos_joined = "N/A" if not defined $mutpos_joined;
			my $number = $type1{$type1}{$type2};
			
			my $perc = $type1{tot}{tot} == 0 ? 0 : $number/$type1{tot}{tot}*100;
			if ($perc == 0) {
				$perc = $perc;
			}
			elsif ($perc < 1) {
				my ($perc2) = $perc =~ /^(0\.0*[1-9][0-9]?)[0-9]*$/;
				$perc = $perc2 if defined $perc2;
			}
			else {
				$perc = int($perc*10+0.5)/10;
			}
			my ($nuc1, $nuc2);
			if ($type1 =~ /^(mis|mat)$/) {# and $type2 =~ /\w+_\w+/) {
				($nuc1, $nuc2) = split("_", $type2);
				DIELOG($outBigLog, "type2 isnt X_X...\n") if $type2 !~ /^.+_.+$/;
				my $orignuc1 = $nuc1;
				$nuc1 = "mis_$orignuc1->$nuc2" if $type1 eq "mis"; #nuc1= "mis_A->C"
				$nuc1 = "mat_$nuc2" if $type1 eq "mat"; #nuc1=mat_A
				$nuc2 = "mis_$orignuc1" if $type1 eq "mis"; #nuc2 = "mis_A"
				$cur{count_atleast1}{$type1} = 1; #single, mat/mis
				$cur{count_atleast1}{$nuc1} = 1; #single, mat_A/mis_A->C
				$cur{count_atleast1}{$nuc2} = 1 if $type1 eq "mis"; #single, mis_A
				$cur{count_total}{$type1} += $number; #single, mat/mis
				$cur{count_total}{$nuc1} += $number; #single, mat_A/mis_A->C
				$cur{count_total}{$nuc2} += $number if $type1 eq "mis"; #single, mis_A
				$cur{percperbp}{$type1} += $perc; #single, mat/mis
				$cur{percperbp}{$nuc1} += $perc; #single, mat_A/mis_A->C
				$cur{percperbp}{$nuc2} += $perc if $type1 eq "mis"; #single, mis_A
			}
			else {
				my $lenTNype = length($type2);
				$lenTNype = 4 if $type1 =~ /^(ins|del|mh)$/ and $lenTNype > 3 and $lenTNype < 10;
				$lenTNype = 4 if $type1 =~ /^(ins|del|mh)$/ and $lenTNype >= 10;
				$nuc1 = $type1 . "_" . $lenTNype;
				$nuc2 = $type2;
				$cur{count_atleast1}{$type1} = 1; #single, del
				$cur{count_atleast1}{$nuc1} = 1; #single, del_4
				$cur{count_total}{$type1} += $number; #del
				$cur{count_total}{$nuc1} += $number; #del_4
				$cur{percperbp}{$type1} += $perc; #del
				$cur{percperbp}{$nuc1} += $perc; #del_4
			}
			my $percnuc = 0;
			if ($nuc2 =~ /^[ACGT]$/) {
				my $totalnuc = $type1{nuc}{$nuc2}; $totalnuc = 0 if not defined $totalnuc;
				$percnuc = $totalnuc == 0 ? 0 : $number/$type1{nuc}{$nuc2}*100;
			}
			my $print1 = "$name\t$readorient1\t$type1\t$nuc1\t$nuc2\t$number\t$perc\t$type1{tot}{tot}\t$data->{igtype}";
			#my $print1 = "$name,$readorient1,$type1,$nuc1,$nuc2,$number,$perc,$type1{tot}{tot},igtype=$data->{igtype},";
			LOG($outBigLog, $print1, $NA);
			LOG($outLog, $print1, "NA");
			print $out1 $print1 if $printout eq 1;
			my $type3 = $type1 eq "ins" ? "inspos" : $type1 eq "del" ? "delpos" : $type1 eq "mat" ? "matpos" : $type1 eq "mis" ? "mispos" : $type1 eq "mh" ? "mhpos" : "";
			if ($type3 ne "") {
#				LOG($outBigLog, " (mean,sd (bp) from junction break: $type1{$type3}{$type2})\n", $NA);
				#LOG($outBigLog, "Undefined mean sd at type1=$type1 type2=$type2 type3=$type3\n") if not defined $type1{$type3}{$type2};
				my ($mean, $sd);
				if (not defined $type1{$type3}{$type2}) {
					($mean, $sd) = ("NA", "NA");
				}
				else {
					($mean, $sd) = split(",", $type1{$type3}{$type2});
				}
				my $print2 = "\t$mean\t$sd";
				LOG($outBigLog, $print2, $NA);
				LOG($outLog, $print2, "NA");
				print $out1 $print2 if $printout eq 1;
			}
			else {
				DIELOG($outBigLog, "shouldn't happen with type1=$type1, type2=$type2, number=$number, type3=$type3\n");
			}
			my $print3 = "\t$data->{chr1}\t$data->{beg1}\t$data->{end1}\t$data->{begJ1}\t$data->{len1}\t$data->{strand1}\t$data->{chr2}\t$data->{beg2}\t$data->{end2}\t$data->{begJ2}\t$data->{len2}\t$data->{strand2}\t$data->{origbeg1}\t$data->{origend1}\t$data->{origbeg2}\t$data->{origend2}\t$mutpos_joined\n";
			LOG($outBigLog, "MUTPOSJOINED UNDEFINED! AT:\n$print3\n\n") if not defined $mutpos_joined;
#			my $print3 = ",$data->{chr1},$data->{begorig1},$data->{endorig1},$data->{begJ1},$data->{len1},$data->{strand1},$data->{chr2},$data->{begorig2},$data->{endorig2},$data->{begJ2},$data->{len2},$data->{strand2}\n";
			LOG($outBigLog, $print3, $NA);
			LOG($outLog, $print3, "NA");
			print $out1 $print3 if $printout eq 1;
		}
	}
	if (defined $type1{R2}) {
		foreach my $type1 (sort keys %{$type1{R2}}) {
			next if $type1 =~ /^(inspos|delpos|matpos|mispos|mhpos|tot|nuc|R2)$/;
			foreach my $type2 (sort keys %{$type1{R2}{$type1}}) {
				my $number = $type1{R2}{$type1}{$type2};
				
				my $perc = $type1{tot}{tot} == 0 ? 0 : $number/$type1{tot}{tot}*100;
				if ($perc == 0) {
					$perc = $perc;
				}
				elsif ($perc < 1) {
					my ($perc2) = $perc =~ /^(0\.0*[1-9][0-9]?)[0-9]*$/;
					$perc = $perc2 if defined $perc2;
				}
				else {
					$perc = int($perc*10+0.5)/10;
				}
				my ($nuc1, $nuc2);
				if ($type1 =~ /^(mis|mat)$/) {# and $type2 =~ /\w+_\w+/) {
					($nuc1, $nuc2) = split("_", $type2);
					DIELOG($outBigLog, "type2 isnt X_X...\n") if $type2 !~ /^.+_.+$/;
					my $orignuc1 = $nuc1;
					$nuc1 = "mis_$orignuc1->$nuc2" if $type1 eq "mis"; #nuc1= "mis_A->C"
					$nuc1 = "mat_$nuc2" if $type1 eq "mat"; #nuc1=mat_A
					$nuc2 = "mis_$orignuc1" if $type1 eq "mis"; #nuc2 = "mis_A"
	#				$cur{count_atleast1}{$type1} = 1; #single, mat/mis
	#				$cur{count_atleast1}{$nuc1} = 1; #single, mat_A/mis_A->C
	#				$cur{count_atleast1}{$nuc2} = 1 if $type1 eq "mis"; #single, mis_A
	#				$cur{count_total}{$type1} += $number; #single, mat/mis
	#				$cur{count_total}{$nuc1} += $number; #single, mat_A/mis_A->C
	#				$cur{count_total}{$nuc2} += $number if $type1 eq "mis"; #single, mis_A
	#				$cur{percperbp}{$type1} += $perc; #single, mat/mis
	#				$cur{percperbp}{$nuc1} += $perc; #single, mat_A/mis_A->C
	#				$cur{percperbp}{$nuc2} += $perc if $type1 eq "mis"; #single, mis_A
				}
				else {
					my $lenTNype = length($type2);
					$lenTNype = 4 if $type1 =~ /^(ins|del|mh)$/ and $lenTNype > 3 and $lenTNype < 10;
					$lenTNype = 4 if $type1 =~ /^(ins|del|mh)$/ and $lenTNype >= 10;
					$nuc1 = $type1 . "_" . $lenTNype;
					$nuc2 = $type2;
	#				$cur{count_atleast1}{$type1} = 1; #single, del
	#				$cur{count_atleast1}{$nuc1} = 1; #single, del_4
	#				$cur{count_total}{$type1} += $number; #del
	#				$cur{count_total}{$nuc1} += $number; #del_4
	#				$cur{percperbp}{$type1} += $perc; #del
	#				$cur{percperbp}{$nuc1} += $perc; #del_4
				}
				my $percnuc = 0;
#				if ($nuc2 =~ /^[ACGT]$/) {
#					my $totalnuc = $type1{nuc}{$nuc2}; $totalnuc = 0 if not defined $totalnuc;
#					$percnuc = $totalnuc == 0 ? 0 : $number/$type1{nuc}{$nuc2}*100;
#				}
				my $print1 = "$name\tR2\t$type1\t$nuc1\t$nuc2\t$number\t$perc\t$type1{tot}{tot}\t$data->{igtype}";
				#my $print1 = "$name,$readorient1,$type1,$nuc1,$nuc2,$number,$perc,$type1{tot}{tot},igtype=$data->{igtype},";
				LOG($outBigLog, $print1, $NA);
				LOG($outLog, $print1, "NA");
				print $out1 $print1 if $printout eq 1;
#				my $type3 = $type1 eq "ins" ? "inspos" : $type1 eq "del" ? "delpos" : $type1 eq "mat" ? "matpos" : $type1 eq "mis" ? "mispos" : $type1 eq "mh" ? "mhpos" : "";
#				if ($type3 ne "") {
##				LOG($outBigLog, " (mean,sd (bp) from junction break: $type1{$type3}{$type2})\n", $NA);
#					my ($mean, $sd) = split(",", $type1{$type3}{$type2});
#					my $print2 = "\t$mean\t$sd";
#					LOG($outBigLog, $print2, $NA);
#					print $out1 $print2 if $printout eq 1;
#				}
#				else {
#					DIELOG($outBigLog, "shouldn't happen with type1=$type1, type2=$type2, number=$number, type3=$type3\n");
#				}
				my $print3 = "\tNA\tNA\t$data->{chr1}\t$data->{beg1}\t$data->{end1}\t$data->{begJ1}\t$data->{len1}\t$data->{strand1}\t$data->{chr2}\t$data->{beg2}\t$data->{end2}\t$data->{begJ2}\t$data->{len2}\t$data->{strand2}\n";
#			my $print3 = ",$data->{chr1},$data->{begorig1},$data->{endorig1},$data->{begJ1},$data->{len1},$data->{strand1},$data->{chr2},$data->{begorig2},$data->{endorig2},$data->{begJ2},$data->{len2},$data->{strand2}\n";
				LOG($outBigLog, $print3, $NA);
				LOG($outLog, $print3, "NA");
				print $out1 $print3 if $printout eq 1;
			}
		}
	}

	$igtypez->{$data->{igtype}} = 1;
	foreach my $typezz1 (sort keys %cur) {
		$typez2->{$typezz1} = 1;
		foreach my $typezz2 (sort keys %{$cur{$typezz1}}) {
			$typez1->{$typezz2} = 1;
			$all->{$data->{igtype}}{$typezz2}{$typezz1} += $cur{$typezz1}{$typezz2};
		}
	}
	LOG($outBigLog, "$name has no A!\n","NA") if not defined $cur{count_atleast1}{"mat_A"};
	LOG($outBigLog, "$name has no G!\n","NA") if not defined $cur{count_atleast1}{"mat_G"};
	LOG($outBigLog, "$name has no C!\n","NA") if not defined $cur{count_atleast1}{"mat_C"};
	LOG($outBigLog, "$name has no T!\n","NA") if not defined $cur{count_atleast1}{"mat_T"};
	LOG($outBigLog, "$name has no total A!\n","NA") if not defined $cur{count_total}{"mat_A"};
	LOG($outBigLog, "$name has no total G!\n","NA") if not defined $cur{count_total}{"mat_G"};
	LOG($outBigLog, "$name has no total C!\n","NA") if not defined $cur{count_total}{"mat_C"};
	LOG($outBigLog, "$name has no total T!\n","NA") if not defined $cur{count_total}{"mat_T"};

	LOG($outLog, "$name has no A!\n","NA") if not defined $cur{count_atleast1}{"mat_A"};
	LOG($outLog, "$name has no G!\n","NA") if not defined $cur{count_atleast1}{"mat_G"};
	LOG($outLog, "$name has no C!\n","NA") if not defined $cur{count_atleast1}{"mat_C"};
	LOG($outLog, "$name has no T!\n","NA") if not defined $cur{count_atleast1}{"mat_T"};
	LOG($outLog, "$name has no total A!\n","NA") if not defined $cur{count_total}{"mat_A"};
	LOG($outLog, "$name has no total G!\n","NA") if not defined $cur{count_total}{"mat_G"};
	LOG($outLog, "$name has no total C!\n","NA") if not defined $cur{count_total}{"mat_C"};
	LOG($outLog, "$name has no total T!\n","NA") if not defined $cur{count_total}{"mat_T"};

	$all->{$data->{igtype}}{total_nuc} += $type1{tot}{tot};
	$all->{$data->{igtype}}{total_read} ++;
	LOG($outBigLog, "\n-------------------------\n",$NA);
	LOG($outLog, "\n-------------------------\n", "NA");
	return($igtypez, $typez1, $typez2, $all);
}


sub get_seqTR {
	my ($data, $name, $outLog, $namewant) = @_;

	my $print_seqTR .= "\n\n>${LGN}get_seqTR${N}:\n";

	# 1. add begT1 (begT1=7 with adapter, so begT1-1 . I)
	my $begAdd = $data->{begT1} - $data->{begTN};
	$begAdd = $begAdd > 0 ? $begAdd . "I" : "";
	$data->{cigarTN} = $begAdd . $data->{cigarQ1};
	$print_seqTR .= "\n1. add begT1 and cigarQ1\n";
	$print_seqTR .= "   begAdd = ${LGN}$begAdd${N}\n";
	$print_seqTR .= "   cigarQ1 = ${LGN}$data->{cigarQ1}${N}\n";
	$print_seqTR .= "   cigarTN: ${LCY}$data->{cigarTN}${N} (added begAdd I + data->cigarQ1 = ${LGN}${begAdd}I${N} . ${LGN}$data->{cigarQ1}${N})\n";
	
	#2. add begT2-endT1. if < 0, then deletion. if > 0, then insertion
	if ($data->{mhend} - $data->{mhbeg} > 0) {
		$print_seqTR .= "\n2. MH exist so add MHlen MH\n";
		$print_seqTR .= "   prev cigarQ2 = ${LGN}$data->{cigarQ2}${N} - MH (${LGN}" . ($data->{mhend} - $data->{mhbeg}) . " bp)${N}\n";
		($data->{cigarQ2}) = fix_cigarQ2(abs($data->{mhend} - $data->{mhbeg}), $data->{cigarQ2});
	}
	if ($data->{gap} > 0) {
		$print_seqTR .= "\n2. GAP exist so add $data->{gap} bp gap\n";
		$print_seqTR .= "   prev cigarQ2 = ${LGN}$data->{cigarQ2}${N} - GAP (${LGN}" . ($data->{mhend} - $data->{mhbeg}) . " bp)${N}\n";
		($data->{cigarQ2}) = "$data->{gap}I" . $data->{cigarQ2};
	}
	$print_seqTR .= "   curr cigarQ2 = ${LGN}$data->{cigarQ2}${N}\n";
	$data->{cigarTN} .= $data->{cigarQ2};
	$print_seqTR .= "   cigarTN: ${LCY}$data->{cigarTN}${N} (added cigarQ2=${LGN}$data->{cigarQ2}${N})\n";

	#3. add lenT - endT1 - 1 . "I"
	$print_seqTR .= "\n3. add end of seq (-) from lenTN (${LGN}$data->{lenTN})${N} - endT2 (${LGN}$data->{endT2}${N}) = ${LGN}" . ($data->{lenTN} - $data->{endT2}) . "${N})\n";
	my $endAdd = ($data->{lenTN} - $data->{endT2});
   	$endAdd = $endAdd > 0 ? $endAdd . "I" : $endAdd < 0 ? $endAdd . "D" : "";
	$print_seqTR .= "   endAdd = ${LGN}$endAdd${N} (I if lenTN-endT2 > 0, or D if < 0)\n";
	$data->{cigarTN} .= $endAdd;
	$print_seqTR .= "   cigarTN = ${LCY}$data->{cigarTN}${N} (added endAdd ${LGN}$endAdd${N})\n";
	$data->{seqTR} = print_cigar($data->{seqTN}, $data->{cigarTN}, $name, $outLog, $namewant);
	$print_seqTR .= "\n${LGN} FINAL cigarTN = ${N} ${LCY}$data->{cigarTN}${N}\n\n";
	LOG($outBigLog, $print_seqTR) if $name eq $namewant;
	return ($data);
}


sub junc2check {
		my ($beg1, $end1, $beg2, $end2, $strand1, $strand2, $junc1, $junc2, $linecount, $line, $outLog) = @_;
	my $junc2check = $strand2 eq "-" ? $end2 : $beg2;
	DIELOG($outBigLog, "Died ${LCY}$inputFile${N} at linecount=$linecount: junc2check $junc2check isn't the same as junc2 $junc2
			junc1 ($junc1) = strand1 ($strand1) eq + ? end1 ($end1) : beg1 ($beg1)
			junc2check ($junc2check) = strand2 ($strand2) eq + ? beg2 ($beg2) : end2 ($end2)
			junc2 = $junc2\n\nline:\n$line\n\n") if $junc2check ne $junc2;
}

sub fix_cigarQ2 {
	my ($del, $cigarQ2) = @_;
	#LOG($outBigLog, "fix cigarQ2 del=$del, cigarQ2=$cigarQ2\n");
	my ($nums2, $alps2, $ciglen2) = parse_cigar($cigarQ2);
	my $newcigar = "";
	for (my $i = 0; $i < @{$alps2}; $i++) {
		my $alp2 = $alps2->[$i];
		my $num2 = $nums2->[$i];
		if ($num2 - $del <= 0) {
			$del -= $num2;
			$del = 0 if $del < 0;
			next;
		}
		else {
			$num2 -= $del;
			$del -= $del;
			$del = 0 if $del < 0;
			$newcigar .= "$num2$alp2";
		}
	}
	#LOG($outBigLog, "del=$del\nbefore = $cigarQ2\nafter  = $newcigar\n\n");
	return ($newcigar);
}

sub print_cigar {
   my ($seq1orig, $cigar, $name, $namewant) = @_;
   my ($nums, $alps, $ciglen) = parse_cigar($cigar);
	my ($seq1, $seq2, $flag) = ("", "", "");
	my @seq1orig = split("", $seq1orig);
	my ($ind1, $ind2, $ind) = (0,0,0);
	my $die = "";
	my $printz = "";
	my $undeftoend = 0;
   for (my $i = 0; $i < @{$nums}; $i++) {
		my $alp = $alps->[$i];
		my $num = $nums->[$i];
		for (my $j = $ind; $j < $ind + $num; $j++) {
			if (not defined ($seq1orig[$j])) {
				$die .= "," if $die ne "";
				$die .= $j-1;
				$undeftoend = 1;
			}
			else {
				$undeftoend = 0;
			}
			my $seq1origchunk = $seq1orig[$j];
			$seq1origchunk = "-" if not defined $seq1origchunk;
			if ($alp eq "M") {
				$seq1 .= $seq1origchunk;
				$seq2 .= $seq1origchunk;
				$flag .= "|";
				$ind1 ++;
				$ind2 ++;
				$printz .= "\t$seq1origchunk M $seq1origchunk\tind=$ind, ind1=$ind1, ind2=$ind2, alp=$alp, num=$num\n"; #if $name eq $namewant;
			}
			elsif ($alp eq "X" or $alp eq "N") {
				$seq1 .= $seq1origchunk;
				$seq2 .= $seq1origchunk;
				$flag .= "m";
				$ind1 ++;
				$ind2 ++;
				$printz .= "\t$seq1origchunk m $seq1origchunk\tind=$ind, ind1=$ind1, ind2=$ind2, alp=$alp, num=$num\n"; #if $name eq $namewant;
			}
			elsif ($alp eq "D") {
				$seq1 .= $seq1origchunk;
				$seq2 .= "-";			
				$ind1 ++;
				$flag .= "D";
				$printz .= "\t$seq1origchunk D -\tind=$ind, ind1=$ind1, ind2=$ind2, alp=$alp, num=$num\n"; #if $name eq $namewant;
			}
			elsif ($alp eq "I") {
				$seq1 .= "-";
				$seq2 .= $seq1origchunk;
				$ind2 ++;
				$flag .= "I";
				$printz .= "\t- I $seq1origchunk\tind=$ind, ind1=$ind1, ind2=$ind2, alp=$alp, num=$num\n"; #if $name eq $namewant;
			}
		}
		$ind += $num;
   }
	LOG($outBigLog, "seqTN  " . colorize($seq1orig) . "\nREF    " . colorize($seq1) . "\nseqTR  " . colorize($seq2) . "\nflag   " . $flag . "\n\n") if $name eq $namewant;
#	LOG($outBigLog, "\n\nDEAD\n\n$printz\n\nname=$name\nind1=$ind1, ind2=$ind2, ind=$ind\n$cigar\nseqTN  " . colorize($seq1orig) . "\nREF    " . colorize($seq1) . "\nseqTR  " . colorize($seq2) . "\nflag   " . $flag . "\n\n") if $die eq 1;
	LOG($outBigLog, "\n\ndied at line 2500+ 3_Mutation.pl:\nUndef=seq1orig[j]\nj=$die\nname=$LCY$name$N\n\n$printz\n\nname=$name\nind1=$ind1, ind2=$ind2, ind=$ind\n$cigar\nseqTN  " . colorize($seq1orig) . "\nREF    " . colorize($seq1) . "\nseqTR  " . colorize($seq2) . "\nflag   " . $flag . "\n\n") if $die ne "" and $undeftoend eq 0;
	return($seq1);
}
__END__
