#!usr/bin/perl
#Takes a bed file of PacBio SMRT isoforms and compares them to a list of validated 5' ends, 3' ends and introns to create a list of validated transcript structures.
#TO'G 5/14/15

use warnings;
use strict;

my ($test_file, $valid_starts_file, $valid_ends_file, $valid_introns_file) = (@ARGV);

#Create an array of validated start sites from the start sites input file:
open(INF, "<$valid_starts_file") or die "couldn't open file";

my @valid_start;

while (my $line = <INF>) {
	chomp($line);
    next if ($line =~ /^track/); #skips the track definition line
	push (@valid_start, $line); #puts each line of the start sites file into an array to be checked later
}
close(INF);

#Check each start site in the SMRT reads file against the array of validated start sites:
open(INF, "<$test_file") or die "couldn't open file";
open(OUT, ">$test_file.valid_start.bed.temp");

my @good_start;
my $new_start_line;

while (my $line = <INF>) {
	chomp($line);
    next if ($line =~ /^track/); #skips the track definition line
	my ($chrom, $chromStart, $chromEnd, $name, $score, $strand, $thickStart, $thickEnd, $itemRgb, $blockCount, $blockSizes, $blockStarts)  = split("\t", $line);
	if ($strand eq "+") {
		foreach my $valid_start (@valid_start) { #checks to see if the 5' end of the (plus strand) SMRT transcript matches a range of possible start site values from the list of validated start sites
			my @start_cols = split("\t", $valid_start);
            my ($range_start, $range_end, $SMRT_depth) = split(":", $start_cols[6]);
			if (($chrom eq $start_cols[0]) and ($strand eq $start_cols[5]) and ($chromStart >= $range_start) and ($chromStart <= $range_end)) {
				$new_start_line = "$line\t$start_cols[1]"; #creates a line for the read, changing the start site to the consensus start site and adding an extra field with the original start site
				push (@good_start, $new_start_line); #if the start site matches, pushes the line into a new array of SMRT transcripts with validated 5' ends
				print OUT $new_start_line, "\n"; #prints out a file of SMRT reads with validated 5' ends
                last;
			}
		}
	}
	elsif ($strand eq "-"){
		foreach my $valid_start (@valid_start) { #checks to see if the 5' end of the (minus strand) SMRT transcript matches a range of possible start site values from the list of validated start sites
			my @start_cols = split("\t", $valid_start);
            my ($range_start, $range_end, $SMRT_depth) = split(":", $start_cols[6]);
			if (($chrom eq $start_cols[0]) and ($strand eq $start_cols[5]) and ($chromEnd >= $range_start) and ($chromEnd <= $range_end)) {
				$new_start_line = "$line\t$start_cols[2]"; #creates a line for the read, changing the start site to the consensus start site and adding an extra field with the original start site
				push (@good_start, $new_start_line); #if the start site matches, pushes the line into a new array of SMRT transcripts with validated 5' ends
				print OUT $new_start_line, "\n"; #prints out a file of SMRT reads with validated 5' ends
                last;
			}
		}
	}
}

my $good_start_number = scalar @good_start;

print "Validated $good_start_number start sites.\n"; #at this point, have output a file of reads (in their original form) that have validated 5' ends, and have an array in memory of reads that have validated 5' ends, and their newly estimated 5' ends

close(OUT);
close(INF);

#Create an array of validated end sites from the end sites input file:
open(INF, "<$valid_ends_file") or die "couldn't open file";

my @valid_end;

while (my $line = <INF>) {
	chomp($line);
    next if ($line =~ /^track/); #skips the track definition line
	push (@valid_end, $line); #puts each line of the end sites file into an array to be checked later
}

close(INF);

open(OUT, ">$test_file.valid_start_and_end.bed.temp");

my @good_start_and_end;
my $new_end_line;

foreach my $good_start (@good_start) { #starts with the array of SMRT transcripts with validated 5' ends
	my ($chrom, $chromStart, $chromEnd, $name, $score, $strand, $thickStart, $thickEnd, $itemRgb, $blockCount, $blockSizes, $blockStarts, $new_coord) = split("\t", $good_start);
	if ($strand eq "+") { #determines 3' end of the SMRT transcript
		foreach my $valid_end (@valid_end) { #checks to see if the 3' end of the SMRT transcript matches a range of possible 3' end values from the list of validated start sites
            my @end_cols = split("\t", $valid_end);
            my ($range_start, $range_end, $SMRT_depth) = split(":", $end_cols[6]);
            if (($chrom eq $end_cols[0]) and ($strand eq $end_cols[5]) and ($chromEnd >= $range_start) and ($chromEnd <= $range_end)) {
                $new_end_line = "$good_start\t$end_cols[2]";
                push (@good_start_and_end, $new_end_line);
                print OUT $new_end_line, "\n"; #prints out a file of SMRT reads with validated 5' and 3' ends
                last;
            }
        }
    }
	elsif ($strand eq "-") {
        foreach my $valid_end (@valid_end) { #checks to see if the 3' end of the SMRT transcript matches a range of possible 3' end values from the list of validated start sites
            my @end_cols = split("\t", $valid_end);
            my ($range_start, $range_end, $SMRT_depth) = split(":", $end_cols[6]);
            if (($chrom eq $end_cols[0]) and ($strand eq $end_cols[5]) and ($chromStart >= $range_start) and ($chromStart <= $range_end)) {
                $new_end_line = "$chrom\t$chromStart\t$chromEnd\t$name\t$score\t$strand\t$thickStart\t$thickEnd\t$itemRgb\t$blockCount\t$blockSizes\t$blockStarts\t$end_cols[1]\t$new_coord";
                push (@good_start_and_end, $new_end_line);
                print OUT $new_end_line, "\n"; #prints out a file of SMRT reads with validated 5' and 3' ends
                last;
            }
        }
    }
}

my $good_start_end_number = scalar @good_start_and_end;

print "Validated $good_start_end_number end sites.\n";

close(OUT);

open(INF, "<$valid_introns_file") or die "couldn't open file";

my @valid_intron;

while (my $line = <INF>) {
	chomp($line);
    next if ($line =~ /^track/); #skips the track definition line
	push (@valid_intron, $line); #creates an array of valid splice junctions
}

close(INF);

open(OUT, ">$test_file.validated.bed.temp");

my $start;
my $end;
my @intron_start;
my @intron_end;
my @intron_coord_pair;
my @good_intron_counter;

foreach my $good_start_and_end (@good_start_and_end) { #starts with the array of SMRT transcripts with validated 5' and 3' ends
	my @cols = split("\t", $good_start_and_end);
	my $intron_strand = $cols[5];
	my $intron_number = $cols[9] - 1;
	if ($intron_number == 0) { #if a SMRT transcript has validated 5' and 3' ends and no introns, it is fully validated
		print OUT $good_start_and_end, "\n";
	}
	else {
		my @block_sizes = split(",", $cols[10]);
		my @block_starts = split(",", $cols[11]);
		for (my $i = 0; $i < $intron_number; $i = $i + 1) { #for the transcript currently in the "while" loop, creates an array of intron start sites relative to the genome
			$start = $cols[1] + $block_sizes[$i] + $block_starts[$i];
			push(@intron_start, $start);
		}
		for (my $i2 = 1; $i2 < $cols[9]; $i2 = $i2 + 1) { #for the transcript currently in the "while" loop, creates an array of intron end sites relative to the genome
			$end = $cols[1] + $block_starts[$i2];
			push(@intron_end, $end);
		}
		for (my $i3 = 0; $i3 < $intron_number; $i3 = $i3 + 1) { #for the transcript currently in the "while" loop, matches up intron start and end sites to create an array of complete intron coordinates relative to the genome
			my $intron_coords = "$intron_start[$i3]:$intron_end[$i3]";
			push (@intron_coord_pair, $intron_coords);
		} 
		@intron_start = ();
		@intron_end = (); #intron starts and ends have been assigned to the @intron_coords array; empty them for the next transcript
		foreach my $intron_coord_pair (@intron_coord_pair) { #goes through each intron in the SMRT transcript
			my @coords = split(":", $intron_coord_pair); #allows extraction of the start and end coordinates from each intron in the SMRT transcript
			foreach my $valid_intron (@valid_intron) { #goes through each intron in the array of validated introns
				my @valid_coords = split("\t", $valid_intron); #allows extraction of the start and end coordinates from each validated intron
				next if $cols[0] ne $valid_coords[0]; #enforces chromosome matching
                next if $intron_strand ne $valid_coords[5]; #enforces strand matching
				if (($coords[0] == $valid_coords[1]) and ($coords[1] == $valid_coords[2])) {
					push(@good_intron_counter, $intron_coord_pair);	#puts introns that are validated for this transcript into an array (this really just functions as a counter)
				}			
			}
		}
		@intron_coord_pair = (); #once each intron in the SMRT transcript has been examined, empty the array for the next transcript
		if (@good_intron_counter == $intron_number) { #check to see if all of the introns in the transcript are validated
			print OUT $good_start_and_end, "\n";
		}
		@good_intron_counter = (); #after checking to see if all the introns in the transcript are validated, empties this array for the next transcript

	}
}

print "Finished validating introns.\n";

close(OUT);

system("sort -k 2,2n -k 3,3n \Q$test_file\E.valid_start.bed.temp > \Q$test_file\E.valid_start.bed");
system("rm \Q$test_file\E.valid_start.bed.temp");

system("sort -k 2,2n -k 3,3n \Q$test_file\E.valid_start_and_end.bed.temp > \Q$test_file\E.valid_start_and_end.bed");
system("rm \Q$test_file\E.valid_start_and_end.bed.temp");

system("sort -k2,2n -k3,3n \Q$test_file\E.validated.bed.temp > \Q$test_file\E.validated.bed"); #This uniq command shouldn't be necessary: need to figure out why I'm getting duplicates (get as many duplicates as there are features in the valid_ends file, but things only match one)
system("rm \Q$test_file\E.validated.bed.temp");

open(INF, "<$test_file.validated.bed");
open(OUT, ">$test_file.validated_corrected.temp");

my $new_block_size;
my @exon_start;
my @exon_end;
my @exon_coord_pair;

while (my $line = <INF>) {
    chomp($line);
    my @cols = split("\t", $line);
    my $exon_number = $cols[9];
    if ($exon_number == 1) { #if a SMRT transcript has validated 5' and 3' ends and no introns, it is fully validated. Just adjust the start and end to the consensus sites and fix the BlockSize accordingly
        $new_block_size = $cols[13]-$cols[12];
		print OUT "$cols[0]\t$cols[12]\t$cols[13]\t$cols[3]\t$cols[4]\t$cols[5]\t$cols[12]\t$cols[13]\t$cols[8]\t$cols[9]\t$new_block_size\t$cols[11]\n";
	}
    else { #need to adjust the start and end sites and also the blockStarts and blockSizes
        my @block_sizes = split(",", $cols[10]);
		my @block_starts = split(",", $cols[11]);
		for (my $i = 0; $i < $exon_number; $i = $i + 1) { #for the transcript currently in the "while" loop, creates an array of exon start sites relative to the genome
			$start = $cols[1] + $block_starts[$i];
			push(@exon_start, $start);
		}
		for (my $i2 = 0; $i2 < $exon_number; $i2 = $i2 + 1) { #for the transcript currently in the "while" loop, creates an array of intron end sites relative to the genome
			$end = $cols[1] + $block_starts[$i2] + $block_sizes[$i2];
			push(@exon_end, $end);
		}
        shift(@exon_start); #removes the first exon start
        unshift(@exon_start, $cols[12]); #replaces the first exon start with the adjust chrStart value
        pop(@exon_end); #removes the last exon end
        push(@exon_end, $cols[13]); #replaces the last exon end with the adjusted chrEnd value
		for (my $i3 = 0; $i3 < $exon_number; $i3 = $i3 + 1) { #for the transcript currently in the "while" loop, matches up intron start and end sites to create an array of complete intron coordinates relative to the genome
			my $exon_coords = "$exon_start[$i3]:$exon_end[$i3]";
			push (@exon_coord_pair, $exon_coords);
		} 
		@exon_start = ();
		@exon_end = (); #intron starts and ends have been assigned to the @intron_coords array; empty them for the next transcript
		my @new_starts;
		my @new_sizes;
        #print $cols[3], "\t", @exon_coord_pair, "\n";
		foreach my $exon_coord_pair (@exon_coord_pair) { #goes through each exon in the SMRT transcript
			my @coords = split(":", $exon_coord_pair);
			my $blockStart = $coords[0] - $cols[12];
			my $blockSize = $coords[1] - $coords[0];
			push (@new_starts, $blockStart);
			push (@new_sizes, $blockSize);
		}
        @exon_coord_pair = ();
        shift(@new_starts); #removes the first value of the new_starts array, so we can replace it with 0 (it won't be 0 already if chrStart has been updated)
		my $assembled_starts = join(",", 0, @new_starts);
		my $assembled_sizes = join(",", @new_sizes);
		print OUT "$cols[0]\t$cols[12]\t$cols[13]\t$cols[3]\t$cols[4]\t$cols[5]\t$cols[12]\t$cols[13]\t$cols[8]\t$cols[9]\t$assembled_sizes\t$assembled_starts\n";
    }
}

#then need to add code to collapse the transcripts with identical structure into a single feature

close(INF);
close(OUT);

system("sort -k 2,2n -k 3,3n -k11,11 -k12,12 -k5,5n \Q$test_file\E.validated_corrected.temp > \Q$test_file\E.validated_corrected.bed");

open(INF, "<$test_file.validated_corrected.bed");
open(OUT, ">$test_file.isoforms.bed");

#To collapse transcripts with matching structures into isoforms. Need to deal with count for first line, make sure last feature prints out and separate plus and minus strands to deal with the slim possibility that antisense transcripts share the same coordinates.

my $prev_chr_plus = "start";
my $prev_chrStart_plus = 0;
my $prev_chrEnd_plus = 0;
my $prev_name_plus;
my $count_plus = 0;
my $prev_rgb_plus;
my $prev_blocks_plus;
my $prev_blockSizes_plus = "1,1";
my $prev_blockStarts_plus = "0,0";
my $prev_chr_minus = "start";
my $prev_chrStart_minus = 0;
my $prev_chrEnd_minus = 0;
my $prev_name_minus;
my $count_minus = 0;
my $prev_rgb_minus;
my $prev_blocks_minus;
my $prev_blockSizes_minus = "1,1";
my $prev_blockStarts_minus = "0,0";

while(my $line = <INF>) {
	chomp($line);
	my @cols = split("\t", $line);
    if ($cols[5] eq "+") {
        if (($cols[0] eq $prev_chr_plus) and ($cols[1] == $prev_chrStart_plus) and ($cols[2] == $prev_chrEnd_plus) and ($cols[10] eq $prev_blockSizes_plus) and ($cols[11] eq $prev_blockStarts_plus)) {
            $count_plus = $count_plus + $cols[4];
            $prev_chr_plus = $cols[0];
            $prev_name_plus = $cols[3];
            $prev_rgb_plus = $cols[8];
            $prev_blocks_plus = $cols[9];
            
        }
        else {
            if ($count_plus == 0) {
                $prev_chrStart_plus = $cols[1];
                $prev_chrEnd_plus = $cols[2];
                $prev_blockSizes_plus = $cols[10];
                $prev_blockStarts_plus = $cols[11];
                $count_plus = $cols[4];
                $prev_chr_plus = $cols[0];
                $prev_name_plus = $cols[3];
                $prev_rgb_plus = $cols[8];
                $prev_blocks_plus = $cols[9];
            }
            else {
                print OUT "$prev_chr_plus\t$prev_chrStart_plus\t$prev_chrEnd_plus\t$prev_name_plus\t$count_plus\t\+\t$prev_chrStart_plus\t$prev_chrEnd_plus\t$prev_rgb_plus\t$prev_blocks_plus\t$prev_blockSizes_plus\t$prev_blockStarts_plus\n";
                $prev_chrStart_plus = $cols[1];
                $prev_chrEnd_plus = $cols[2];
                $prev_blockSizes_plus = $cols[10];
                $prev_blockStarts_plus = $cols[11];
                $count_plus = $cols[4];
                $prev_chr_plus = $cols[0];
                $prev_name_plus = $cols[3];
                $prev_rgb_plus = $cols[8];
                $prev_blocks_plus = $cols[9];
            }
        }
    }
    elsif ($cols[5] eq "-") {
        if (($cols[0] eq $prev_chr_minus) and ($cols[1] == $prev_chrStart_minus) and ($cols[2] == $prev_chrEnd_minus) and ($cols[10] eq $prev_blockSizes_minus) and ($cols[11] eq $prev_blockStarts_minus)) {
            $count_minus = $count_minus + $cols[4];
            $prev_chr_minus = $cols[0];
            $prev_name_minus = $cols[3];
            $prev_rgb_minus = $cols[8];
            $prev_blocks_minus = $cols[9];
            
        }
        else {
            if ($count_minus == 0) {
                $prev_chrStart_minus = $cols[1];
                $prev_chrEnd_minus = $cols[2];
                $prev_blockSizes_minus = $cols[10];
                $prev_blockStarts_minus = $cols[11];
                $count_minus = $cols[4];
                $prev_chr_minus = $cols[0];
                $prev_name_minus = $cols[3];
                $prev_rgb_minus = $cols[8];
                $prev_blocks_minus = $cols[9];
            }
            else {
                print OUT "$prev_chr_minus\t$prev_chrStart_minus\t$prev_chrEnd_minus\t$prev_name_minus\t$count_minus\t\-\t$prev_chrStart_minus\t$prev_chrEnd_minus\t$prev_rgb_minus\t$prev_blocks_minus\t$prev_blockSizes_minus\t$prev_blockStarts_minus\n";
                $prev_chrStart_minus = $cols[1];
                $prev_chrEnd_minus = $cols[2];
                $prev_blockSizes_minus = $cols[10];
                $prev_blockStarts_minus = $cols[11];
                $count_minus = $cols[4];
                $prev_chr_minus = $cols[0];
                $prev_name_minus = $cols[3];
                $prev_rgb_minus = $cols[8];
                $prev_blocks_minus = $cols[9];
            }
        }
    }
}
if ($count_plus > 0) {#prints out the last feature (plus strand)
    print OUT "$prev_chr_plus\t$prev_chrStart_plus\t$prev_chrEnd_plus\t$prev_name_plus\t$count_plus\t\+\t$prev_chrStart_plus\t$prev_chrEnd_plus\t$prev_rgb_plus\t$prev_blocks_plus\t$prev_blockSizes_plus\t$prev_blockStarts_plus\n";
}

if ($count_minus > 0) {#prints out the last feature (minus strand)
    print OUT "$prev_chr_minus\t$prev_chrStart_minus\t$prev_chrEnd_minus\t$prev_name_minus\t$count_minus\t\-\t$prev_chrStart_minus\t$prev_chrEnd_minus\t$prev_rgb_minus\t$prev_blocks_minus\t$prev_blockSizes_minus\t$prev_blockStarts_minus\n";
}

#print OUT "$cols[0]\t$cols[1]\t$cols[2]\t$cols[3]\t$count\t$cols[5]\t$cols[6]\t$cols[7]\t$cols[8]\t$cols[9]\t$cols[10]\t$cols[11]\n"; Need to print out last feature.

close(INF);
close(OUT);

system("rm \Q$test_file\E.validated_corrected.temp");