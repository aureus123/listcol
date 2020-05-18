#!/usr/bin/perl -w
use strict;

# Generate instances of Set 1

print "gen1 - Generate instances of Set 1\n\n";

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
				my $name = "Set1/R$n" . "_" . $r . "_" . $q . "_" . $p;

				$result = system("listcol/gengraph $name.graph $n $p >$name.log 2>/dev/null")>>8;
				if ($result) { die "Error in gengraph."; }

				$result = system("listcol/genrandominst $name $k 1 3 $r $q >>$name.log 2>/dev/null")>>8;
				if ($result) { die "Error in genrandominst."; }
			}
		}
	}
}

exit 0;
