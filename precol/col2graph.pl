#!/usr/bin/perl -w
use strict;

# Conversor de grafos de formato DIMACS a formato propio.
#
# (c) 2007 Daniel E. Severin.
#
# Véase en UTF-8 unicode.

print "col2grafo - Conversor de formatos (DIMACS -> grafo)\n\n";

sub manual {
	print "Uso: ./col2grafo archivo.col archivo.graph archivo.list\n";
	exit 1;
}

if ($#ARGV < 2) { manual; }

my $archivo_in = $ARGV[0];
my $archivo_out = $ARGV[1];
my $archivo_list = $ARGV[2];

print "Archivo de entrada: $archivo_in\n";
open IN, '<', $archivo_in or die "Error al leer este archivo: $!";
print "Archivo de salida: $archivo_out\n";
open OUT, '>', $archivo_out or die "Error al crear este archivo: $!";
print "Archivo de lista: $archivo_list\n";
open LST, '>', $archivo_list or die "Error al crear este archivo: $!";

my $vertices = -1;
my $aristas = -1;
my $cuenta = 0;
my $i;

while (<IN>) {
	chomp;
	my @datos = split ' ', $_;
	if ($datos[0] eq "p" && $datos[1] eq "edges") {
		# Imprime la cantidad de vértices y aristas.
		$vertices = $datos[2];
		$aristas = $datos[3];
		print OUT "$vertices:$aristas\n";
		print LST "$vertices:*****\n";
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
	if ($datos[0] eq "f") {
		my $length = (scalar @datos) - 2;
		print LST "$length";
		for ($i = 1; $i <= $length; $i++) {
			my $co = $datos[$i+1]-1;
			print LST " $co";
		}
		print LST "\n";
	}	
}

close OUT;
close IN;

print "POR FAVOR COLOQUE LA CANTIDAD DE COLORES EN EL ARCHIVO LIST SEGUNDO PARAMETRO!\n";

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
