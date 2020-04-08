#!/usr/bin/perl -w
use strict;

# Conversor de grafos de formato DIMACS a formato propio.
#
# (c) 2007 Daniel E. Severin.
#
# Véase en UTF-8 unicode.

print "col2grafo - Conversor de formatos (DIMACS -> grafo)\n\n";

sub manual {
	print "Uso: ./col2grafo archivo.col archivo.grafo\n";
	exit 1;
}

if ($#ARGV < 1) { manual; }

my $archivo_in = $ARGV[0];
my $archivo_out = $ARGV[1];

print "Archivo de entrada: $archivo_in\n";
open IN, '<', $archivo_in or die "Error al leer este archivo: $!";
print "Archivo de salida: $archivo_out\n";
open OUT, '>', $archivo_out or die "Error al crear este archivo: $!";

my $vertices = -1;
my $aristas = -1;
my $cuenta = 0;

while (<IN>) {
	chomp;
	my @datos = split ' ', $_;
	if ($datos[0] eq "p" && $datos[1] eq "edge") {
		# Imprime la cantidad de vértices y aristas.
		$vertices = $datos[2];
		$aristas = $datos[3];
		print OUT "$vertices:$aristas\n";
	}
	if ($datos[0] eq "e") {
		my $v1 = $datos[1]-1;
		my $v2 = $datos[2]-1;
		if ($v1 < $v2) {
			# Imprime cada arista.
			print OUT "$v1,$v2\n";
			$cuenta++;
		}
	}
}

close OUT;
close IN;

if ($vertices == -1 || $aristas == -1) {
	print "No se ha encontrado cantidad de vértices o aristas!\n";
	exit 1;
}
if ($cuenta != $aristas) {
	print "No coincide la cantidad de aristas con la informacion del archivo ";
	print "(cuenta=$cuenta, aristas=$aristas).\n";
	print "POR FAVOR COLOQUE $cuenta EN EL REGISTRO DE ARISTAS!\n";
	exit 1;
}

print "Cantidad de vértices: $vertices\n";
print "Cantidad de aristas: $aristas\n";
print "Listop!\n";
exit 0;
