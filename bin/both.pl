#!/usr/bin/perl

use strict; use warnings; use mitochy; use Getopt::Std; use FAlite;
use vars qw($opt_v $opt_i $opt_n);
getopts("vi:n:");

my $inputFolder = $opt_i;
die "\nusage: $YW$0$N ${LGN}[-n namewant ]$N -i $CY<folder containing.50.50.final.tsv>$N\n\n" unless defined $opt_i and -e $opt_i;
my $namewant = $opt_n;

open (my $out1, ">", "both.csv") or die "Cannot write to both.csv: $!\n";
my @mut = qw(mat mis del ins mh);# mhinsdel mhins mhdel insdel mismhinsdel mismh misins misdel);
my $muts = \@mut;
my $mutslen = @mut;
my @nucs = qw(A C G T);
my $nucs = \@nucs;
my %head;
foreach my $nuc1 (@{$nucs}[0..@{$nucs}-1]) {
	$head{"mat_perkb_$nuc1"} = 1;
	$head{"mat_$nuc1"} = 1;
	foreach my $nuc2 (@{$nucs}[0..@{$nucs}-1]) {
		next if $nuc1 eq $nuc2;
		$head{"mis_perkb_$nuc1->$nuc2"} = 1;
		$head{"mis_$nuc1->$nuc2"} = 1;
	}
}
for (my $i = 1; $i < 5; $i++) {
	$head{"ins_$i"} = 1;

	$head{"del_$i"} = 1;

	$head{"mh_$i"} = 1;
}

##############
### HEADER ###
#############

##
# any

my $done;

my $header1 = "";
$header1 .= "name,sampleID,sampleType,sampleNo,sample,treat,cut,isotype,myfill2,len,blunt,any1,any2,any3";
print $out1 "name,sampleID,sampleType,sampleNo,sample,treat,cut,isotype,myfill2,len,blunt,any1,any2,any3";
#$header1 .= "name,sampleID,sample,sampleNo,sample,cut,isotype,len";
#print $out1 "name,sampleID,sample,sampleNo,sample,cut,isotype,len";
# single mutation detail
foreach my $mut (sort keys %head) {
	next if grep(/^$mut$/, @mut);
	print $out1 ",$mut\_any";
	$header1 .= ",$mut";
}
# single mutation any
undef $done;
foreach my $mut1 (@{$muts}[0..$mutslen-1]) {
	my $mut = $mut1;
	my @currmuts = ($mut1);
	$done = add2done(\@currmuts, $done, 1, 1);

	print $out1 ",$mut\_any";
	$header1 .= ",$mut\_any";

}

# double mutation any
undef $done;
foreach my $mut1 (@{$muts}[0..$mutslen-1]) {
	foreach my $mut2 (@{$muts}[0..$mutslen-1]) {
		next if $mut2 eq $mut1;
		my $mut = $mut1 . $mut2;
		next if defined $done->{$mut};
		my @currmuts = ($mut1, $mut2);
		$done = add2done(\@currmuts, $done, 2, 1);

		print $out1 ",$mut\_any";
		$header1 .= ",$mut\_any";
	}
}
# triple mutation any
undef $done;
foreach my $mut1 (@{$muts}[0..$mutslen-1]) {
	foreach my $mut2 (@{$muts}[0..$mutslen-1]) {
		next if $mut2 eq $mut1;
		foreach my $mut3 (@{$muts}[0..$mutslen-1]) {
			next if $mut3 eq $mut1 or $mut3 eq $mut2;
			my $mut = $mut1 . $mut2 . $mut3;
			next if defined $done->{$mut};
			my @currmuts = ($mut1, $mut2, $mut3);
			$done = add2done(\@currmuts, $done, 3, 1);

			print $out1 ",$mut\_any";
			$header1 .= ",$mut\_any";
		}
	}
}

##
# only
# single mutation only
foreach my $mut1 (@{$muts}[0..$mutslen-1]) {
	my $mut = $mut1;
	next if defined $done->{$mut};
	my @currmuts = ($mut1);
	$done = add2done(\@currmuts, $done, 1, 1);

	print $out1 ",$mut\_only";
	$header1 .= ",$mut\_only";
}

# double mutation only
undef $done;
foreach my $mut1 (@{$muts}[0..$mutslen-1]) {
	foreach my $mut2 (@{$muts}[0..$mutslen-1]) {
		next if $mut2 eq $mut1;
		my $mut = $mut1 . $mut2;
		next if defined $done->{$mut};
		my @currmuts = ($mut1, $mut2);
		$done = add2done(\@currmuts, $done, 2, 1);
		print $out1 ",$mut\_only";
		$header1 .= ",$mut\_only";
	}
}

# triple mutation only
undef $done;
foreach my $mut1 (@{$muts}[0..$mutslen-1]) {
	foreach my $mut2 (@{$muts}[0..$mutslen-1]) {
		next if $mut2 eq $mut1;
		foreach my $mut3 (@{$muts}[0..$mutslen-1]) {
			next if $mut3 eq $mut1 or $mut3 eq $mut2;
			my $mut = $mut1 . $mut2 . $mut3;
			next if defined $done->{$mut};
			my @currmuts = ($mut1, $mut2, $mut3);
			$done = add2done(\@currmuts, $done, 3, 1);
			print $out1 ",$mut\_only";
			$header1 .= ",$mut\_only";
		}
	}
}

print $out1 "\n";
print "$header1\n";
my @inputFiles = <$inputFolder/*50.50.final.tsv>;
for (my $m = 0; $m < @inputFiles; $m++) {

undef $namewant if not defined $opt_n;
my $input1 = $inputFiles[$m];
my ($folder1, $fileName1) = mitochy::getFilename($input1, "folderfull");
die "input must be a filter.tlx.50.50.final.tsv!\n" if $input1 !~ /filter.tlx.50.50.final.tsv$/;


#$header1 .= "name,sampleID,sampleNo,sample,treat,sample,cut,isotype,myfill2,len";
my ($sampleID) = $fileName1 =~ /^.*([WRSD][0-9](cut)?)/; die "can't get samp;leID from $LCY$fileName1$N\n" if not defined $sampleID;
my ($sampleType, $sampleNo) = $fileName1 =~ /^.*([WRSD])([0-9])/; die "can't get samp;le from $LCY$fileName1$N\n" if not defined $sampleType or not defined $sampleNo;
my ($cut) = $fileName1 =~ /$sampleType${sampleNo}cut/ ? "CUT" : "UNCUT";
$cut = "CUT" if $fileName1 =~ /[WDSR]2/;
my $treat = $sampleType eq "W" ? "1. WT" : $sampleType eq "S" ? "2. Setx KO" : $sampleType eq "R" ? "3. RNH2b KO" : $sampleType eq "D" ? "4. Setx RNH2b DKO" : "5. Unknown";
my $sample = $sampleNo eq 1 ? "IgG1 undigested germline" : $sampleNo eq 3 ? "IgG3 undigested enriched" : $sampleNo eq 2 ? "IgG1 digested enriched" : "Unknown sampledescription";
$sample =~ s/undigested/digested/g if $cut eq "CUT";
my $myfill = "$sampleNo\\n$sample\\n$cut\\nBAIT=IgM\\n";
my %name;
my %len;
my %data;
my $printwant = "";
my %count;
my %info;
print "\n$LGN$m/" . scalar(@inputFiles) . "$N. Processing $LCY$input1$N\n";
print "--------------------\n";
print " -> sampleID   = $LCY$sampleID$N\n";
print " -> sampleNo   = $LCY$sampleNo$N\n";
print " -> sampleType = $LCY$sampleType$N\n";
print " -> treat      = $LCY$treat$N\n";
print " -> sample     = $LCY$sample$N\n";
print " -> cut        = $LCY$cut$N\n";
print " -> myfill     = $LCY$myfill$N\n";
open (my $in1, "<", $input1) or die "Cannot read from $input1: $!\n";
while (my $line = <$in1>) {
	chomp($line);
	next if $line =~ /^#/;
	next if $line =~ /^name\t/;
	my @arr = split("\t", $line);
	my $name = $arr[0];
	next if $arr[1] eq "R2";
	$namewant = $name if not defined $namewant;
	$data{mat}{$name} = 0 if not defined	$data{mat}{$name};
	$data{mis}{$name} = 0 if not defined	$data{mis}{$name};
	$data{del}{$name} = 0 if not defined	$data{del}{$name};
	$data{ins}{$name} = 0 if not defined	$data{ins}{$name};
	$data{mh}{$name}  = 0 if not defined	$data{mh}{$name};
	my $isotype = $arr[8]; $isotype =~ s/^\d+\. (.+)$/$1/;
	$info{$name}{isotype} = $isotype;
	$info{$name}{len} = $arr[7] if not defined $info{$name}{len};
	my $nucdone;
	foreach my $nuc1 (@{$nucs}[0..@{$nucs}-1]) {
		$count{"mat_perkb_$nuc1"}{$name} = 0 if not defined $count{"mat_$nuc1"}{$name};
		$count{"mat_$nuc1"}{$name} = 0 if not defined $count{"mat_$nuc1"}{$name};
		$count{"mat_perkb_$nuc1"}{$name} += $arr[6]*10 if $arr[3] eq "mat_$nuc1";
		$count{"mat_$nuc1"}{$name} += $arr[5] if $arr[3] eq "mat_$nuc1";
		foreach my $nuc2 (@{$nucs}[0..@{$nucs}-1]) {
			next if $nuc1 eq $nuc2;
			$count{"mis_perkb_$nuc1->$nuc2"}{$name} = 0 if not defined $count{"mis_$nuc1->$nuc2"}{$name};
			$count{"mis_$nuc1->$nuc2"}{$name} = 0 if not defined $count{"mis_$nuc1->$nuc2"}{$name};
			$count{"mis_perkb_$nuc1->$nuc2"}{$name} += $arr[6]*10 if $arr[3] eq "mis_$nuc1->$nuc2";
			$count{"mis_$nuc1->$nuc2"}{$name} += $arr[5] if $arr[3] eq "mis_$nuc1->$nuc2";
		}
	}
	for (my $i = 1; $i < 5; $i++) {
		$count{"ins_$i"}{$name} = 0 if not defined $count{"ins_$i"}{$name};
		$count{"ins_$i"}{$name} = 1 if $arr[3] eq "ins_$i";

		$count{"del_$i"}{$name} = 0 if not defined $count{"del_$i"}{$name};
		$count{"del_$i"}{$name} = 1 if $arr[3] eq "del_$i";

		$count{"mh_$i"}{$name} = 0 if not defined $count{"mh_$i"}{$name};
		$count{"mh_$i"}{$name} = 1 if $arr[3] eq "mh_$i";
	}
	$data{mat}{$name} = 1 if $arr[2] eq "mat";
	$data{mis}{$name} = 1 if $arr[2] eq "mis";
	$data{del}{$name} = 1 if $arr[2] eq "del";
	$data{ins}{$name} = 1 if $arr[2] eq "ins";
	$data{mh}{$name}  = 1 if $arr[2] eq "mh";
	$name{$name} = 1;
	$printwant .= "   $YW$name:$N\n" if $name eq $namewant and $printwant eq "";
	$printwant .= "    -> mut=$LCY$arr[2]$N ($LCY$arr[3]$N) datavalue=$LGN$data{$arr[2]}{$name}$N count=$LGN$LGN$count{$arr[3]}{$name}$N/len=$LCY$info{$name}{len}$N\n" if $name eq $namewant;
=comment
	my @arr2;
	@arr2 = @arr; $arr2[2] = "mat";
	$name{$name}{mat} = join("\t", @arr2);
	@arr2 = @arr; $arr2[2] = "mis";
	$name{$name}{mis} = join("\t", @arr2);
	@arr2 = @arr; $arr2[2] = "ins";
	$name{$name}{ins} = join("\t", @arr2);
	@arr2 = @arr; $arr2[2] = "del";
	$name{$name}{del} = join("\t", @arr2);
	@arr2 = @arr; $arr2[2] = "mh";
	$name{$name}{mh} = join("\t", @arr2);
#	@arr2 = @arr; $arr2[2] = "mhins";
#	$name{$arr[0]}{mhins} = join("\t", @arr2);
#	@arr2 = @arr; $arr2[2] = "mhdel";
#	$name{$arr[0]}{mhdel} = join("\t", @arr2);
#	@arr2 = @arr; $arr2[2] = "mhinsdel";
#	$name{$arr[0]}{mhinsdel} = join("\t", @arr2);
#	@arr2 = @arr; $arr2[2] = "mhinsdel";
#	$name{$arr[0]}{mhinsdel} = join("\t", @arr2);
#	@arr2 = @arr; $arr2[2] = "mhinsdel";
#	$name{$arr[0]}{mhinsdel} = join("\t", @arr2);

=cut
}
close $in1;
print "   Example:\n   $LCY" . $printwant . "$N\n";
#my @mut = qw(insdel mhins mhdel mhinsdel);

##############
### VALUES ###
##############

# any, values
#print $out1 "name,sampleID,sample,sampleNo,sample,cut,isotype,len";
foreach my $name (sort keys %name) {
	my $len = $info{$name}{len};
	my $isotype = $info{$name}{isotype};
	my $myfill2 = $myfill . "PREY at $isotype";
	my $printout .=  "$name,$sampleID,$sampleType,$sampleNo,$sample,$treat,$cut,$isotype,$myfill2,$len";
	my $printhead2 = $printout;

	#blunt
	if ($data{ins}{$name} eq 0 and $data{del}{$name} eq 0 and $data{mh}{$name} eq 0) {
		$printout .= ",1";
	}
	else {
		$printout .= ",0";
	}

	my $total_mut_read = 0;
	foreach my $mut (sort @mut) {
		$total_mut_read ++ if $data{$mut}{$name} eq 1 and $mut =~ /^(ins|mh|del)$/;
	}

	#at least 1
	if ($total_mut_read >= 1) {
		$printout .= ",1";
	}
	else {
		$printout .= ",0";
	}

	#at least 2
	if ($total_mut_read >= 2) {
		$printout .= ",1";
	}
	else {
		$printout .= ",0";
	}

	#at least 3
	if ($total_mut_read >= 3) {
		$printout .= ",1";
	}
	else {
		$printout .= ",0";
	}


	# single mutation detail
	undef $done;
	foreach my $mut (sort keys %count) {
		next if grep(/^$mut$/, @mut);
		if ($mut !~ /perkb/) {
			$printout .= ",1" if $count{$mut}{$name} ne 0;
			$printout .= ",0" if $count{$mut}{$name} eq 0;
		}
		else {
			$printout .= ",$count{$mut}{$name}";# if $count{$mut}{$name} ne 0;
		}
		$printhead2 .= ",$mut";
	}
	# single mutation any
	undef $done;
	foreach my $mut1 (@{$muts}[0..$mutslen-1]) {
		my $mut = $mut1;
		next if defined $done->{$mut};
		my @currmuts = ($mut1);
		$done = add2done(\@currmuts, $done, 1);

		if ($data{$mut1}{$name} eq 1) {
			$printout .=  ",1";
		}
		else {
			$printout .=  ",0";
		}
		$printhead2 .= ",$mut";
	}

	# double mutation any
	undef $done;
	foreach my $mut1 (@{$muts}[0..$mutslen-1]) {
		foreach my $mut2 (@{$muts}[0..$mutslen-1]) {

			next if $mut2 eq $mut1;
			my $mut = $mut1 . $mut2;
			next if defined $done->{$mut};
			my @currmuts = ($mut1, $mut2);
			$done = add2done(\@currmuts, $done, 2);

			if ($data{$mut1}{$name} eq 1 and $data{$mut2}{$name} eq 1) {
				$printout .=  ",1";
			}
			else {
				$printout .=  ",0";
			}
			$printhead2 .= ",$mut";
		}
	}
	# triple mutation any
	undef $done;
	foreach my $mut1 (@{$muts}[0..$mutslen-1]) {
		foreach my $mut2 (@{$muts}[0..$mutslen-1]) {
			next if $mut2 eq $mut1;
			foreach my $mut3 (@{$muts}[0..$mutslen-1]) {
				next if $mut3 eq $mut1 or $mut3 eq $mut2;
				my $mut = $mut1 . $mut2 . $mut3;
				next if defined $done->{$mut};
				my @currmuts = ($mut1, $mut2, $mut3);
				$done = add2done(\@currmuts, $done, 3);
				
				if ($data{$mut1}{$name} eq 1 and $data{$mut2}{$name} eq 1 and $data{$mut3}{$name} eq 1) {
					$printout .=  ",1";
				}
				else {
					$printout .=  ",0";
				}
				$printhead2 .= ",$mut";
			}
		}
	}

	#only, values
	# single mutation only
	undef $done;
	foreach my $mut1 (@{$muts}[0..$mutslen-1]) {
		my $mut = $mut1;
		next if defined $done->{$mut};
		my @currmuts = ($mut1);
		$done = add2done(\@currmuts, $done, 1);

		my $unique = 0;
		if ($data{$mut1}{$name} ne 1) {
			$printout .=  ",0";
			$printhead2 .= ",$mut";
			next;
		}
		foreach my $mutcheck (@{$muts}[0..$mutslen-1]) {
			next if $mut1 eq $mutcheck;
			$unique ++ if $data{$mutcheck}{$name} eq 1; #if other mutation exist, next
			last if $unique ne 0;
		}
		if ($data{$mut1}{$name} eq 1 and $unique eq 0) {
			$printout .=  ",1";
		}
		else {
			$printout .=  ",0";
		}
		$printhead2 .= ",$mut";
	}
	# double mutation only
	undef $done;
	foreach my $mut1 (@{$muts}[0..$mutslen-1]) {
		foreach my $mut2 (@{$muts}[0..$mutslen-1]) {
			next if $mut1 eq $mut2;
			my $mut = $mut1 . $mut2;
			next if defined $done->{$mut};
			my @currmuts = ($mut1, $mut2);
			$done = add2done(\@currmuts, $done, 2);

			my $unique = 0;
			if ($data{$mut1}{$name} ne 1 or $data{$mut2}{$name} ne 1) {
				$printout .=  ",0";
				$printhead2 .= ",$mut";
				next;
			}
			foreach my $mutcheck (@{$muts}[0..$mutslen-1]) {
				next if $mut1 eq $mutcheck or $mut2 eq $mutcheck;
				$unique ++ if $data{$mutcheck}{$name} eq 1; #if other mutation exist, next
				last if $unique ne 0;
			}
			if ($data{$mut1}{$name} eq 1 and $data{$mut2}{$name} eq 1 and $unique eq 0) {
				$printout .=  ",1";
			}
			else {
				$printout .=  ",0";
			}
			$printhead2 .= ",$mut";
		}
	}

	# triple mutation only
	undef $done;
	foreach my $mut1 (@{$muts}[0..$mutslen-1]) {
		foreach my $mut2 (@{$muts}[0..$mutslen-1]) {
			next if $mut1 eq $mut2;
			foreach my $mut3 (@{$muts}[0..$mutslen-1]) {
				next if $mut1 eq $mut3 or $mut2 eq $mut3;
				my $mut = $mut1 . $mut2 . $mut3;
				next if defined $done->{$mut};
				my @currmuts = ($mut1, $mut2, $mut3);
				$done = add2done(\@currmuts, $done, 3);

				my $unique = 0;
				if ($data{$mut1}{$name} ne 1 or $data{$mut2}{$name} ne 1 or $data{$mut3}{$name} ne 1) {
					$printout .=  ",0";
					$printhead2 .= ",$mut";
					next;
				}
				foreach my $mutcheck (@{$muts}[0..$mutslen-1]) {
					next if $mut1 eq $mutcheck or $mut2 eq $mutcheck or $mut3 eq $mutcheck;
					$unique ++ if $data{$mutcheck}{$name} ne 0; #if other mutation exist, print 0
					last if $unique ne 0;
				}
				if ($data{$mut1}{$name} eq 1 and $data{$mut2}{$name} eq 1 and $data{$mut3}{$name} eq 1 and $unique eq 0) {
					$printout .=  ",1";
				}
				else {
					$printout .=  ",0";
				}
				$printhead2 .= ",$mut";
			}
		}
	}
	$printhead2 .= "\n";
	$printout .=  "\n";
	
	print "   $YW$printout$N" if $name eq $namewant;
#	print $out1 $printhead2 if $name eq $namewant;
	print $out1 $printout;
}
	print " -> Processed $LGN" . (keys %name) . "$N reads\n\n";
	print "--------------------\n";
}
close $out1;

sub add2done {
	my ($muts, $done, $level, $print) = @_;
	undef $print;
	my $mutslen = @{$muts};
	print "\nPrinting add2done\n" if defined $print;
	foreach my $mut1 (@{$muts}[0..$mutslen-1]) {
		next if $mut1 eq "mat" and $level eq 0;
		if ($level eq 1) {
			my $mut = $mut1;
			$done->{$mut} = 1;
			print " -> $LCY$mut1$N\n" if defined $print;
			next;
		}
		foreach my $mut2 (@{$muts}[0..$mutslen-1]) {
			next if $mut2 eq "mat" and $level eq 0;
			next if $mut1 eq $mut2;
			if ($level eq 2) {
				my $mut = $mut1 . $mut2;
				$done->{$mut} = 1;
				print " -> $LCY$mut1$LGN$mut2$N\n" if defined $print;
				next;
			}
			foreach my $mut3 (@{$muts}[0..$mutslen-1]) {
				next if $mut3 eq "mat" and $level eq 0;
				next if $mut1 eq $mut3 or $mut2 eq $mut3;
#				$done->{"$mut1$mut2$mut3"} = 1;
				if ($level eq 3) {
					my $mut = $mut1 . $mut2 . $mut3;
					$done->{$mut} = 1;
					print " -> $LCY$mut1$LGN$mut2$LPR$mut3$N\n" if defined $print;
					next;
				}
				foreach my $mut4 (@{$muts}[0..$mutslen-1]) {
					next if $mut4 eq "mat" and $level eq 0;
					next if $mut1 eq $mut4 or $mut2 eq $mut4 or $mut3 eq $mut4;
					if ($level eq 4) {
						my $mut = $mut1 . $mut2 . $mut3 . $mut4;
						$done->{$mut} = 1;
						print " -> $LCY$mut1$LGN$mut2$LPR$mut3$N\n" if defined $print;
						next;
					}
					if ($level eq 0) {
						my $mut = $mut1 . $mut2 . $mut3 . $mut4;
						$done->{$mut} = 1;
						print " -> $LCY$mut1$LGN$mut2$LPR$mut3$N\n" if defined $print;
						next;
					}
				}
			}
		}
	}
	return $done;
}

__END__
	foreach my $mut1 (@{$muts}[0..$mutslen-1]) {
		foreach my $mut2 (@mut[0..$mutslen-1]) {
			next if $mut2 eq $mut1;
			my $mut = $mut1 . $mut2;
			next if defined $done2a{$mut};
			$done2a{"$mut1$mut2"} = 1;
			$done2a{"$mut2$mut1"} = 1;
			if ($name{$name}{$mut1} ne 1 or $name{$name}{$mut2} ne 1) {
				print $out1 ",0";
			}
			else {
				my $unique = 0;
				foreach my $mutcheck (@mut[0..$mutslen-1]) {
					next if $mut1 eq $mutcheck or $mut2 eq $mutcheck;
					$unique ++ if $name{$name}{$mutcheck} ne 0; #if other mutation exist, next
					last if $unique ne 0;
				}
				if ($name{$name}{$mut1} eq 1 and $name{$name}{$mut2} eq 1 and $unique eq 0) {
					print $out1 ",1";
				}
				else {
					print $out1 ",0";
				}
				foreach my $mut3 (@mut[0..$mutslen-1]) {
					next if $mut3 eq $mut1 or $mut3 eq $mut2;
					my $mut = $mut1 . $mut2 . $mut3;
					next if defined $done2a{$mut};
					$done2a{"$mut1$mut2$mut3"} = 1;
					$done2a{"$mut1$mut3$mut2"} = 1;
					$done2a{"$mut2$mut1$mut3"} = 1;
					$done2a{"$mut2$mut3$mut1"} = 1;
					$done2a{"$mut3$mut2$mut1"} = 1;
					$done2a{"$mut3$mut1$mut2"} = 1;
					if ($name{$name}{$mut1} ne 1 or $name{$name}{$mut2} ne 1 or $name{$name}{$mut3} ne 1) {
						print $out1 ",0";
					}
					else {
						my $unique2 = 0;
						foreach my $mutcheck (@mut[0..@mut-1]) {
							next if $mut1 eq $mutcheck or $mut2 eq $mutcheck or $mut3 eq $mutcheck;
							$unique2 ++ if $name{$name}{$mutcheck} eq 1; #if other mutation exist, next
							last if $unique2 ne 0;
						}
						if ($name{$name}{$mut1} eq 1 and $name{$name}{$mut2} eq 1 and $name{$name}{$mut3} eq 1 and $unique2 eq 0) {
							print $out1 ",1";
						}
						else {
							print $out1 ",0";
						}
					}
				}
			}
		}
	}
	print $out1 "\n";
}
	print $out1 $name{$name}{mis} . "\t";
	print $out1 $name{$name}{del} . "\t";
	print $out1 $name{$name}{ins} . "\t";
	print $out1 $name{$name}{mh} . "\t";
	print $out1 "1\t" if $data{ins}{$name} == 1 and $data{del}{$name} == 1 and $data{mh}{$name} == 1;
	print $out1 "1\t" if $data{ins}{$name} == 1 and $data{del}{$name} == 1 and $data{mh}{$name} != 1;
#	print $out1 $name{$name}{mhins} . "\n" if $data{del}{$name} == 1 and $data{mh}{$name} == 1;
#	print $out1 $name{$name}{mhdel} . "\n" if $data{ins}{$name} == 1 and $data{del}{$name} == 1;
#	print $out1 $name{$name}{insdel} . "\n" if $data{ins}{$name} == 1 and $data{del}{$name} == 1;
}
close $out1;

__END__
name    type    nuc1    nuc2    number  perc    percnuc1        percnuc2        mean    sd      chr1    beg1    end1    junc1   length1 strand1 chr2    beg2     end2    junc2   length2 strand2
M02034:514:000000000-J6BC6:1:1101:23105:3064    R1      mat     mat_A   A       23      23      100     03. IgG1        2.73    0.31    chr12   114664625114664909       114664809       284     -       chr12   114574599       114574883       114574782       284     -
M02034:514:000000000-J6BC6:1:1101:23105:3064    R1      mat     mat_C   C       22      22      100     03. IgG1        -7.13   0.32    chr12   114664625114664909       114664809       284     -       chr12   114574599       114574883       114574782       284     -
M02034:514:000000000-J6BC6:1:1101:23105:3064    R1      mat     mat_G   G       35      35      100     03. IgG1        1.97    0.27    chr12   114664625114664909       114664809       284     -       chr12   114574599       114574883       114574782       284     -

