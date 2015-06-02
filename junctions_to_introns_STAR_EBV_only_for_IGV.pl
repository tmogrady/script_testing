#!/usr/bin/perl
#junctions_to_introns_STAR.pl by TO'G. Takes junctions.bed files from STAR and converts them into intron bed files separated into EBV and cellular files
use warnings;
use strict;

foreach my $file (@ARGV) {
	open(INF, "<$file" ) or die "couldn't open file";
	open(OUT, ">$file.introns_EBV_only.bed");

	my $line;
	while( $line = <INF> ) {
		
		chomp($line);
		my @cols = split( "\t", $line );	
			
		tr/12/+-/ foreach ($cols[3]);	#change the numeric strand indicators to + or -
			
		next if $cols[0] ne "chrEBV(Akata_107955to171322_1to107954)"; #skip lines that aren't EBV
		
		my $chrStart = $cols[1]-1; #fixes the start coordinate, which is off by 1
				
		print OUT "$cols[0]\t$chrStart\t$cols[2]\t$cols[4]\t$cols[6]\t$cols[3]\n";

		}
		
		close(OUT);
		close(INF);
}