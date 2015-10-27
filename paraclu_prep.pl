#!/usr/bin/perl

#Takes a SAM file of deepCAGE (or other Illumina) data and counts reads with 5' ends at each genomic position to prepare an input file for Paraclu
#Strand assignment for paired-end sequencing assumes TruSeq stranded prep. Reports the 5' end of each READ, not each fragment. Includes only reads flagged as primary alignments.

#USAGE:
# perl PATH/paraclu_prep.pl /PATH/inputfile(s).sam

#TO'G 9/9/15

use warnings;
use strict;

die "USAGE: 'perl <PATH/paraclu_filter.pl> </PATH/sam_file.sam>'" unless @ARGV == 1;

foreach my $file(@ARGV) {
    
    print "preparing input file...\n";
    
    system("sort -k 3,3 -k 4,4n \Q$file\E > \Q$file\E.sorted.temp");
    system("awk '\$2==0 \|\| \$2==81 \|\| \$2==83 \|\| \$2==89 \|\| \$2==137 \|\| \$2==161 \|\| \$2==163' \Q$file\E.sorted.temp > \Q$file\E.sorted.plus.sam.temp");
    system("awk '\$2==16 \|\| \$2==73 \|\| \$2==97 \|\| \$2==99 \|\| \$2==145 \|\| \$2==147 \|\| \$2==153' \Q$file\E.sorted.temp > \Q$file\E.sorted.minus.sam.temp");
    
    #processing of plus sam file
    
    print "processing plus strand...\n";
    
    open(INF, "<$file.sorted.plus.sam.temp") or die "couldn't open file";
    open(OUT, ">$file.plus_read_starts.txt") or die "couldn't open file";
    
    my $previous_coordinate=1;
	my $count=0;
	my $previous_chr = "start";
	
	while (my $line = <INF>) {
		chomp($line);
        next if ($line =~ m/^@/); #skips header lines
        my @cols = split("\t", $line);
			
        if (($cols[2] eq $previous_chr) and ($cols[3] == $previous_coordinate)) {
            $count++; #increases the count by 1
			}
            
        else {
            if ($previous_chr eq "start") { #doesn't print out the placeholder first line.
                $previous_chr = $cols[2];	#sets the previous chromosome, previous coordinate and count values
                $previous_coordinate = $cols[3];
                $count = 1;
            }
				
            else {
                print OUT $previous_chr, "\t+\t", $previous_coordinate, "\t", $count, "\n"; #prints to output file
                $previous_chr = $cols[2];
                $previous_coordinate = $cols[3];
                $count = 1;
            }
        }
	}
	print OUT $previous_chr, "\t+\t", $previous_coordinate, "\t", $count, "\n"; #prints the last start coordinates to output file
	close(INF);
	close(OUT);

	system("rm \Q$file\E.sorted.plus.sam.temp");
	
    
    #processing of MINUS sam file
    
    print "processing minus strand...\n";
    
	open(INF, "<$file.sorted.minus.sam.temp") or die "couldn't open file";
	open(OUT, ">$file.sorted.minus.sam.read_starts.txt.temp") or die "couldn't open file";
	
	my @dist;
	my $sum;
	my %minus_starts;
	
	while (my $line = <INF>) {
		chomp($line);
        my @cols = split("\t", $line);
		
        while ($cols[5] =~ /(\d+)[DMNX=]/g) { #these lines use the CIGAR string to determine the downstream coordinate
            push (@dist, $1);
        }
        $sum += $_ for @dist;
        my $start_coord = $cols[3] + $sum - 1;
			
        my $chr_start_coord = "$cols[2]\:$start_coord"; #combines the chromosome and 5' end coordinate into a key to use for the hash
			
        $sum = 0;
        @dist = ();
			
        if (exists $minus_starts{$chr_start_coord}) { #if the key is already in the hash, increases the value (count) by the read depth for that putative isoform
            $minus_starts{$chr_start_coord} = $minus_starts{$chr_start_coord} + 1;
        }
            
        else {
            $minus_starts{$chr_start_coord} = 1; #if the key is not already in the hash, adds it with a value (count) of the read depth for that putative isoform
        }
	}
    
    foreach my $chr_start_coord (sort keys %minus_starts) { #prints out a(n inadequately) sorted file
		my @split_keys = split("\:", $chr_start_coord);
		print OUT "$split_keys[0]\t-\t$split_keys[1]\t$minus_starts{$chr_start_coord}\n";
	}
	close(INF);
	close(OUT);
    
    system("sort -k 1,1 -k 3,3n \Q$file\E.sorted.minus.sam.read_starts.txt.temp > \Q$file\E.minus_read_starts.txt");
	system("rm \Q$file\E.sorted.minus.sam.read_starts.txt.temp");
	system("rm \Q$file\E.sorted.minus.sam.temp");
    system("rm \Q$file\E.sorted.temp");
}