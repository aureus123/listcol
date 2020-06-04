#!/usr/bin/perl -w
use strict;

# Generate instances

print "stat - Write statistics (over Set 0)\n\n";

sub manual {
	print "Usage: ./stat.pl algo exp\n";
	print "  algo is the algorithm (vc,st,bp,...)\n";
	print "  exp is the number of experiment: 0 (test).\n"; 
	exit 1;
}

if ($#ARGV < 1) { manual; }

my $algo = $ARGV[0];

my $experiment = $ARGV[1];
if ($experiment < 0 || $experiment > 0) { die "exp must be 0."; }

my $n;
my $k;
my $r;
my $q;
my $p;
my $result;

print "n,r,q,p,solved,rel,time\n";
for ($n = 10; $n <= 30; $n+=10) {
	$k = $n / 2;
	for ($r = 1; $r <= 5; $r+=2) {
		for ($q = 25; $q <= 75; $q+=25) {
			for ($p = 25; $p <= 75; $p+=25) {
				my $graph = "R$n" . "_" . $r . "_" . $q . "_" . $p;
				my $name = "Set0/$graph";
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

					print "$n,$r," . ($q/100) . "," . ($p/100) . ",";
					if ($solved == 0) { print "NO,-,-\n"; }
					if ($solved == 2) { print "INF,-,$time\n"; }
					if ($solved == 3) { print "YES,$LB,$time\n"; }
				}
			}
		}
	}
}

my $t;
print "\nn,p,r,t,solved,rel,time\n";
for ($n = 10; $n <= 30; $n+=10) {
	$k = $n / 2;
	for ($p = 25; $p <= 75; $p+=25) {
		for ($r = 1; $r <= 5; $r+=2) {
			for ($t = 5; $t <= 20; $t+=5) {
				if ($t != 15) {
					my $graph = "m$n" . "_" . $p . "_" . $r . "_" . $t;
					my $name = "Set0/$graph";
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

						print "$n," . ($p/100) . ",$r," . ($t/100) . ",";
						if ($solved == 0) { print "NO,-,-\n"; }
						if ($solved == 2) { print "INF,-,$time\n"; }
						if ($solved == 3) { print "YES,$LB,$time\n"; }
					}
				}
			}
		}
	}
}

exit 0;
