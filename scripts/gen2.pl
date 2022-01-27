#!/usr/bin/perl -w
use strict;

# Generate instances of Set 2

print "gen2 - Generate instances of Set 2\n\n";

my $n;
my $k;
my $p;
my $r;
my $t;
my $result;

for ($n = 40; $n <= 100; $n+=20) {
	$k = $n / 2;
	for ($p = 25; $p <= 75; $p+=25) {
		for ($r = 1; $r <= 5; $r+=2) {
			for ($t = 5; $t <= 20; $t+=5) {
				if ($t != 15) {

					my $name = "Set2/m$n" . "_" . $p . "_" . $r . "_" . $t;

					print "Generating $name ...\n";

					do {

						$result = system("listcol/gengraph $name.graph $n $p >$name.log 2>/dev/null")>>8;
						if ($result) { die "Error in gengraph."; }

						$result = system("listcol/genmuinst $name $k $r $t >>$name.log 2>/dev/null")>>8;
						print "Res: $result\n";
						if ($result) { print "An error has occurred. Trying again...\n"; }
					} while ($result);
				}
			}
		}
	}
}

exit 0;
