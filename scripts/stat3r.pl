#!/usr/bin/perl -w
use strict;

# Generate instances

print "stat - Write statistics (over Set 0)\n\n";

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

my $n = 200;
my $k;
my $r;
my $q;
my $p;
my $result;

print "k,r,q,p,solved,rel,time\n";
for ($k = 50; $k <= 150; $k+=50) {
	for ($r = 1; $r <= 5; $r+=2) {
		for ($q = 25; $q <= 75; $q+=25) {
			for ($p = 25; $p <= 75; $p+=25) {
				my $graph = "W$k" . "_" . $r . "_" . $q . "_" . $p;
				my $name = "Set3/$graph";
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
					my $time = $datas[4];
					close HANDLE;

					print "$k,$r," . ($q/100) . "," . ($p/100) . ",";
					if ($solved == 0) { print "NO,-,-\n"; }
					if ($solved == 2) { print "INF,-,$time\n"; }
					if ($solved == 3) {
						print "YES,$LB,$time\n";
						#						system("cat $galgo.log | grep \"Number of columns\" | cut -d\" \" -f5");
					}
				}
			}
		}
	}
}

exit 0;
