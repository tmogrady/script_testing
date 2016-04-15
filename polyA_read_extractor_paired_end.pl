#USAGE:
# perl /PATH/polyA_read_extractor.pl number_of_end_As_in_a_row "name_of_viral_chromosome" /PATHtofiles/*
	#output will be in home directory.
	
#This script works with sam files of aligned paired-end RNA-seq data from libraries prepared with the TruSeq stranded method (first read in pair at the 3' end of the fragment).
	
#EKF (Flemington lab) 11/19/2014. Modified by TO'G 12/3/2014
use warnings;use strict;        

my ($min_As, $virus_chr_term, @files) = @ARGV;
                       
foreach my $file(@files) {
	
	open(INF, "<$file") or die "couldn't open input file";
	open(OUT1, ">$file.$virus_chr_term.polyA_ends") or die "couldn't open output file";
	open(OUT2, ">$file.cellular_polyA_ends") or die "couldn't open output file";
	
	while (my $line = <INF>) {
		chomp($line);
		
		my @split_line = split("\t", $line);
		
		next if ($split_line[0] eq "\@HD" || $split_line[0] eq "\@PG" || $split_line[0] eq "\@SQ"); #skips SAM file header lines
		
		if ($split_line[1] == 81 || $split_line[1] == 83 || $split_line[1] == 89) {  #selects reads with FLAG codes indicating they are first in pair on the plus strand
					
			if (($split_line[5] =~ m/[$min_As-9]S$/) || ($split_line[5] =~ m/\d\dS$/) and ($split_line[9] =~ m/A{$min_As}$/)) { # selects reads with a run of As at the end
				if ($split_line[2] =~ m/\Q$virus_chr_term\E/) {
				
					print OUT1 $line, "\n";
				}
			
				else {
					print OUT2 $line, "\n";
				}

			}
			else { next; }
		}
		
		elsif ($split_line[1] == 73 || $split_line[1] == 97 || $split_line[1] == 99) {  #selects reads with FLAG codes indicating they are first in pair on the minus strand
					
			if (($split_line[5] =~ m/^[$min_As-9]S/) || ($split_line[5] =~ m/^\d\dS/) and ($split_line[9] =~ m/^T{$min_As}/)) { #selects reads with a run of Ts at the beginning
				if ($split_line[2] =~ m/\Q$virus_chr_term\E/) {
				
					print OUT1 $line, "\n";
				}
			
				else {
					print OUT2 $line, "\n";
				}

			}
			else { next;}
		}
			
		else { next; }
		
	}
	
	close(INF);
	close(OUT1);
	close(OUT2);

	system ("sort -k 3,3 -k 4,4n \Q$file\E.\Q$virus_chr_term\E.polyA_ends > \Q$file\E.\Q$virus_chr_term\E.polyA_ends.sorted.sam");
	system ("sort -k 3,3 -k 4,4n \Q$file\E.cellular_polyA_ends > \Q$file\E.cellular_polyA_ends.sorted.sam");
	
	system ("rm \Q$file\E.\Q$virus_chr_term\E.polyA_ends");
	system ("rm \Q$file\E.cellular_polyA_ends");
}