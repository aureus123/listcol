#!/usr/bin/perl -w
use strict;

# Check again the output files

print "recheck - Check again the output files\n\n";

sub manual {
	print "Usage: ./recheck.pl algo exp\n";
	print "  algo is the algorithm (vc,st,bp,...)\n";
	print "  exp is the number of experiment: 0 (test).\n"; 
	exit 1;
}

if ($#ARGV < 1) { manual; }

my $algo = $ARGV[0];

my $experiment = $ARGV[1];
if ($experiment < 0 || $experiment > 0) { die "exp must be 0."; }

my $thread = 1;
my $n;
my $k;
my $r;
my $q;
my $p;
my $result;

for ($n = 10; $n <= 20; $n+=10) {
	$k = $n / 2;
	for ($r = 1; $r <= 5; $r+=2) {
		for ($q = 25; $q <= 75; $q+=25) {
			for ($p = 25; $p <= 75; $p+=25) {
				my $graph = "R$n" . "_" . $r . "_" . $q . "_" . $p;
				my $name = "Set0/$graph";
				my $galgo = "exp$experiment/$graph" . "." . $algo ;

				print "Checking instance $name with $algo...\n";

				my $flag = 0;
				stat("$galgo.out") or $flag = 1;
				if ($flag == 0) {
					system("cp $name.graph t$thread.graph");
					system("cp $name.list t$thread.list");
					system("cp $name.cost t$thread.cost");
					system("rm t$thread.sol 2>/dev/null");
					system("rm t$thread.out 2>/dev/null");

					my $flag2 = 0;
					stat("$name.sol") or $flag2 = 1;
					if ($flag2 == 0) { system("cp $name.sol t$thread.sol"); }

					$result = system("cp $galgo.out t$thread.out");

					$result = system("listcol/checker t$thread")>>8;
					if ($result) { die "Error in checker."; }
					$flag2 = 0;
					stat("t$thread.sol") or $flag2 = 1;
					if ($flag2 == 0) { system("mv t$thread.sol $name.sol"); }
				}
			}
		}
	}
}

my $t;
for ($n = 10; $n <= 20; $n+=10) {
	$k = $n / 2;
	for ($p = 25; $p <= 75; $p+=25) {
		for ($r = 1; $r <= 5; $r+=2) {
			for ($t = 5; $t <= 20; $t+=5) {
				if ($t != 15) {
					my $graph = "m$n" . "_" . $p . "_" . $r . "_" . $t;
					my $name = "Set0/$graph";
					my $galgo = "exp$experiment/$graph" . "." . $algo ;

					print "Checking instance $name with $algo...\n";

					my $flag = 0;
					stat("$galgo.out") or $flag = 1;
					if ($flag == 0) {
						system("cp $name.graph t$thread.graph");
						system("cp $name.list t$thread.list");
						system("cp $name.cost t$thread.cost");
						system("rm t$thread.sol 2>/dev/null");
						system("rm t$thread.out 2>/dev/null");

						my $flag2 = 0;
						stat("$name.sol") or $flag2 = 1;
						if ($flag2 == 0) { system("cp $name.sol t$thread.sol"); }

						$result = system("cp $galgo.out t$thread.out");

						$result = system("listcol/checker t$thread")>>8;
						if ($result) { die "Error in checker."; }
						$flag2 = 0;
						stat("t$thread.sol") or $flag2 = 1;
						if ($flag2 == 0) { system("mv t$thread.sol $name.sol"); }
					}
				}
			}
		}
	}
}


exit 0;

