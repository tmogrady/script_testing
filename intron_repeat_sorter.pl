#!/usr/bin/perl
#intron_repeat_sorted.pl by TO'G. Takes intron bed files and sorts them into a file of introns that have a donor and/or acceptor in an Akata repeat region and those that don't.
use warnings;
use strict;

foreach my $file (@ARGV) {
	open(INF, "<$file" ) or die "couldn't open file";
	open(OUT1, ">$file.repeat_junctions.bed");
	open(OUT2, ">$file.no_repeat_junctions.bed");
				
	while(my $line = <INF>) {
		chomp($line);
		my @cols = split( "\t", $line );
		next if ($cols[0] ne "chrEBV(Akata_107955to171322_1to107954)");		#ignore lines from chromosomes other than EBV
		if (($cols[1]>30543 && $cols[1]<30921) || ($cols[1]>32356 && $cols[1]<34874) || ($cols[1]>34864 && $cols[1]<35921) || ($cols[1]>59496 && $cols[1]<59658) || ($cols[1]>61183 && $cols[1]<63319) || ($cols[1]>70797 && $cols[1]<71307) || ($cols[1]>72285 && $cols[1]<72400) || ($cols[1]>75265 && $cols[1]<98628) || ($cols[1]>99662 && $cols[1]<99782) || ($cols[1]>100506 && $cols[1]<100560) || ($cols[1]>101558 && $cols[1]<103096) || ($cols[1]>103633 && $cols[1]<104677) || ($cols[1]>120666 && $cols[1]<120913) || ($cols[1]>121369 && $cols[1]<121414) || ($cols[1]>133220 && $cols[1]<133363) || ($cols[1]>134541 && $cols[1]<134631) || ($cols[1]>141051 && $cols[1]<141388) || ($cols[1]>144956 && $cols[1]<145026) || ($cols[1]>145029 && $cols[1]<145055) || ($cols[1]>145134 && $cols[1]<145160) || ($cols[1]>145160 && $cols[1]<145238) || ($cols[1]>145238 && $cols[1]<145309) || ($cols[1]>145319 && $cols[1]<145397) || ($cols[1]>145397 && $cols[1]<145468) || ($cols[1]>145644 && $cols[1]<145731) || ($cols[1]>145731 && $cols[1]<145818) || ($cols[1]>148270 && $cols[1]<148387) || ($cols[1]>150810 && $cols[1]<150993) || ($cols[1]>151248 && $cols[1]<151521) || ($cols[1]>153321 && $cols[1]<153420) || ($cols[1]>158984 && $cols[1]<159620) || ($cols[2]>30543 && $cols[2]<30921) || ($cols[2]>32356 && $cols[2]<34874) || ($cols[2]>34864 && $cols[2]<35921) || ($cols[2]>59496 && $cols[2]<59658) || ($cols[2]>61183 && $cols[2]<63319) || ($cols[2]>70797 && $cols[2]<71307) || ($cols[2]>72285 && $cols[2]<72400) || ($cols[2]>75265 && $cols[2]<98628) || ($cols[2]>99662 && $cols[2]<99782) || ($cols[2]>100506 && $cols[2]<100560) || ($cols[2]>101558 && $cols[2]<103096) || ($cols[2]>103633 && $cols[2]<104677) || ($cols[2]>120666 && $cols[2]<120913) || ($cols[2]>121369 && $cols[2]<121414) || ($cols[2]>133220 && $cols[2]<133363) || ($cols[2]>134541 && $cols[2]<134631) || ($cols[2]>141051 && $cols[2]<141388) || ($cols[2]>144956 && $cols[2]<145026) || ($cols[2]>145029 && $cols[2]<145055) || ($cols[2]>145134 && $cols[2]<145160) || ($cols[2]>145160 && $cols[2]<145238) || ($cols[2]>145238 && $cols[2]<145309) || ($cols[2]>145319 && $cols[2]<145397) || ($cols[2]>145397 && $cols[2]<145468) || ($cols[2]>145644 && $cols[2]<145731) || ($cols[2]>145731 && $cols[2]<145818) || ($cols[2]>148270 && $cols[2]<148387) || ($cols[2]>150810 && $cols[2]<150993) || ($cols[2]>151248 && $cols[2]<151521) || ($cols[2]>153321 && $cols[2]<153420) || ($cols[2]>158984 && $cols[2]<159620)) {
			print OUT1 $line, "\n";
		}
		else {
			print OUT2 $line, "\n";
		}
	}

	close(INF);
	close(OUT1);
	close(OUT2);
}
	

