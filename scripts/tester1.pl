#!/usr/bin/perl -w
use strict;

# Generate instances

print "tester - Test an algorithm\n\n";

sub manual {
	print "Usage: ./tester.pl algo exp\n";
	print "  algo is the algorithm (vc,st,bp,...)\n";
	print "  exp is the number of experiment: 0 (test), 1 to 5 (real).\n"; 
	exit 1;
}

if ($#ARGV < 1) { manual; }

my $algo = $ARGV[0];

my $experiment = $ARGV[1];
if ($experiment < 0 || $experiment > 5) { die "exp must be between 0 and 5."; }

my $thread = 1;
my $n;
my $k;
my $r;
my $q;
my $p;
my $result;

for ($n = 40; $n <= 100; $n+=20) {
	$k = $n / 2;
	for ($r = 1; $r <= 5; $r+=2) {
		for ($q = 25; $q <= 75; $q+=25) {
			for ($p = 25; $p <= 75; $p+=25) {
				my $graph = "R$n" . "_" . $r . "_" . $q . "_" . $p;
				my $name = "Set1/$graph";
				my $galgo = "exp$experiment/$graph" . "." . $algo ;

				print "Solving instance $name with $algo...\n";

				my $flag = 0;
				stat("$galgo.out") or $flag = 1;
				if ($flag == 1) {
					system("cp $name.graph t$thread.graph");
					system("cp $name.list t$thread.list");
					system("cp $name.cost t$thread.cost");
					system("rm t$thread.sol 2>/dev/null");

					my $flag2 = 0;
					stat("$name.sol") or $flag2 = 1;
					if ($flag2 == 0) { system("cp $name.sol t$thread.sol"); }

					$result = system("listcol/$algo t$thread >t$thread.log 2>/dev/null")>>8;

					$flag2 = 0;
					stat("t$thread.out") or $flag2 = 1;
					if ($flag2 == 1) {
						print "Error in resolution of $name with $algo\n";
						system("echo 0:-99999999:99999999:-1:0:0 >t$thread.out");
					}
					else {
						$result = system("listcol/checker t$thread")>>8;
						if ($result) { die "Error in checker."; }
						$flag2 = 0;
						stat("t$thread.sol") or $flag2 = 1;
						if ($flag2 == 1) { system("mv t$thread.sol $name.sol"); }
					}

					system("mv t$thread.out $galgo.out");
					system("mv t$thread.log $galgo.log");
				}
			}
		}
	}
}

exit 0;

