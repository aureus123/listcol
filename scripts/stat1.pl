#!/usr/bin/perl -w
use strict;

# Generate instances

print "stat1 - Write statistics (over Set 1)\n\n";

sub manual {
	print "Usage: ./stat.pl algo exp\n";
	print "  algo is the algorithm (vc,st,bp,...)\n";
	print "  exp is the number of experiment: 0 (test), 1 to 5 (real).\n"; 
	exit 1;
}

if ($#ARGV < 1) { manual; }

my $algo = $ARGV[0];

my $experiment = $ARGV[1];
if ($experiment < 0 || $experiment > 5) { die "exp must be between 0 and 5."; }

my $n;
my $k;
my $r;
my $q;
my $p;
my $result;

print "n,r,q,p,solved,gap,nodes,time\n";
for ($n = 40; $n <= 100; $n+=20) {
	$k = $n / 2;
	for ($r = 1; $r <= 5; $r+=2) {
		for ($q = 25; $q <= 75; $q+=25) {
			for ($p = 25; $p <= 75; $p+=25) {
				my $graph = "R$n" . "_" . $r . "_" . $q . "_" . $p;
				my $name = "Set1/$graph";
				my $galgo = "exp$experiment/$graph" . "." . $algo ;

				my $flag = 0;
				stat("$galgo.out") or $flag = 1;
				if ($flag == 0) {
					open HANDLE, '<', "$galgo.out" or die "Error reading output file.";
					$_ = <HANDLE>;
					chomp;
					my @datas = split ':', $_;
					my $solved = $datas[0];
					my $LB = $datas[1];
					my $UB = $datas[2];
					my $nodes = $datas[3];
					my $time = $datas[4];
					close HANDLE;

					print "$n,$r," . ($q/100) . "," . ($p/100) . ",";
					if ($solved == 0) {
						print "NO,";
						if ($LB <= 0 || $UB >= 99999999) { print "-,"; }
						else { my $gap = 100*($UB-$LB)/$UB; print "$gap,"; }
						print "-,-\n";
					}
					if ($solved == 1) { print "YES,0.0,$nodes,$time\n"; }
					if ($solved == 2) { print "INF,-,-,$time\n"; }
				}
			}
		}
	}
}

exit 0;

