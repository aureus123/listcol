#!/usr/bin/perl -w
use strict;

# Generate instances of Set 0 (test)

print "gen0 - Generate instances of Set 0 (test)\n\n";

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
				my $name = "Set0/R$n" . "_" . $r . "_" . $q . "_" . $p;

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

my $t;
for ($n = 10; $n <= 20; $n+=10) {
	$k = $n / 2;
	for ($p = 25; $p <= 75; $p+=25) {
		for ($r = 1; $r <= 5; $r+=2) {
			for ($t = 5; $t <= 20; $t+=5) {
				if ($t != 15) {

					my $name = "Set0/m$n" . "_" . $p . "_" . $r . "_" . $t;

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
