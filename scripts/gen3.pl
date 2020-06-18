#!/usr/bin/perl -w
use strict;

# Generate instances of Set 3

print "gen3 - Generate instances of Set 3\n\n";

my $n = 200;
my $k;
my $r;
my $q;
my $p;
my $result;

for ($k = 50; $k <= 150; $k+=50) {
	for ($r = 1; $r <= 5; $r+=2) {
		for ($q = 25; $q <= 75; $q+=25) {
			for ($p = 25; $p <= 75; $p+=25) {
				my $name = "Set3/W$k" . "_" . $r . "_" . $q . "_" . $p;

				print "Generating $name ...\n";

				do {

					$result = system("listcol/gengraph $name.graph $n $p >$name.log 2>/dev/null")>>8;
					if ($result) { die "Error in gengraph."; }

					$result = system("listcol/genrandominst $name $k 1 3 $r $q >>$name.log 2>/dev/null")>>8;
					print "Res: $result\n";
					if ($result) { print "An error has occurred. Trying again...\n"; }
				} while ($result);
			}
		}
	}
}

exit 0;
