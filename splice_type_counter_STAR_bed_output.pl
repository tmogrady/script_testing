#!/usr/bin/perl
#W_splice_type_counter.pl by TO'G. Takes SJ.out.tab intron files from STAR alignments to chrEBV(Akata_107955to171322_1to107954) and counts depth of coverage for W1-W2, W2-W1 and W1-W1 junctions.
use warnings;
use strict;

foreach my $file (@ARGV) {
    open(INF, "<$file" ) or die "couldn't open file";
    open(OUT, ">$file.junctions_of_interest_depths.txt");
    open(OUT1, ">$file.W1W2_junctions.bed");
    open(OUT2, ">$file.W2W1_junctions.bed");
    open(OUT3, ">$file.W1W1_junctions.bed");
    open(OUT4, ">$file.W2W2_junctions.bed");
    open(OUT5, ">$file.C2W1_junctions.bed");
    open(OUT6, ">$file.LP12_junctions.bed");
    open(OUT7, ">$file.W2Y1_junctions.bed");
    open(OUT8, ">$file.W1_BHRF_junctions.bed");
    open(OUT9, ">$file.BZLF2_junctions.bed");
    open(OUT10, ">$file.LMP2_junctions.bed");

    my $out_line;
    my $W1W2count = 0;
    my $W2W1count = 0;
    my $W1W1count = 0;
    my $W2W2count = 0;
    my $C2W1count = 0;
    my $LP12count = 0;
    my $W2Y1count = 0;
    my $W1BHRF1count = 0;
            
    while(my $line = <INF>) {
        chomp($line);
        my @cols = split( "\t", $line );
        if ($cols[3] == 1) {
            $out_line = "$cols[0]\t$cols[1]\t$cols[2]\t$cols[4]:$cols[6]:$cols[7]:$cols[8]\t$cols[6]\t+"
        }
        elsif ($cols[3] == 2) {
            $out_line = "$cols[0]\t$cols[1]\t$cols[2]\t$cols[4]:$cols[6]:$cols[7]:$cols[8]\t$cols[6]\t-"
        }
        else { next }
        next if ($cols[0] ne "chrEBV(Akata_107955to171322_1to107954)");		#ignore lines from chromosomes other than EBV
        if (($cols[1] == 77886||$cols[1] == 80959||$cols[1] == 84032||$cols[1] == 87105||$cols[1] == 90178||$cols[1] == 93251||$cols[1] == 96324) && ($cols[2] == 77966||$cols[2] == 81039||$cols[2] == 84112||$cols[2] == 87185||$cols[2] == 90258||$cols[2] == 93331||$cols[2] == 96404))
        {
            $W1W2count = $W1W2count + $cols[7];
            print OUT1 $out_line, "\n";
        }
        elsif (($cols[1] == 78099||$cols[1] == 81172||$cols[1] == 84245||$cols[1] == 87318 ||$cols[1] == 90391||$cols[1] == 93464||$cols[1] == 96537) && ($cols[2] == 77819||$cols[2] == 80892||$cols[1] == 83965||$cols[2] == 87038||$cols[2] == 90111||$cols[2] == 93184||$cols[2] == 96257))
        {
            $W2W1count = $W2W1count + $cols[7];
            print OUT2 $out_line, "\n";
        }
        elsif (($cols[1] == 77886||$cols[1] == 80959||$cols[1] == 84032||$cols[1] == 87105||$cols[1] == 90178||$cols[1] == 93251||$cols[1] == 96324) && ($cols[2] == 77819||$cols[2] == 80892||$cols[2] == 83965||$cols[2] == 87038||$cols[2] == 90111||$cols[2] == 93184||$cols[2] == 96257))
        {
            $W1W1count = $W1W1count + $cols[7];
            print OUT3 $out_line, "\n";
        }
        elsif (($cols[1] == 78099||$cols[1] == 81172||$cols[1] == 84245||$cols[1] == 87318||$cols[1] == 90391||$cols[1] == 93464||$cols[1] == 96537) && ($cols[2] == 77966||$cols[2] == 81039||$cols[2] == 84112||$cols[2] == 87185||$cols[2] == 90258||$cols[2] == 93331||$cols[2] == 96404))
        {
            $W2W2count = $W2W2count + $cols[7];
            print OUT4 $out_line, "\n";
        }
        elsif ($cols[1] == 74923 && ($cols[2] == 77819||$cols[2] == 80892||$cols[2] == 83965||$cols[2] == 87038||$cols[2] == 90111||$cols[2] == 93184||$cols[2] == 96257))
        {
            $C2W1count = $C2W1count + $cols[7];
            print OUT5 $out_line, "\n";
        }
        elsif (($cols[1] == 77676||$cols[1] == 80749||$cols[1] == 83822||$cols[1] == 86895||$cols[1] == 89968||$cols[1] == 93041||$cols[1] == 96114) && ($cols[2] == 77824||$cols[2] == 80897||$cols[2] == 83970||$cols[2] == 87043||$cols[2] == 90116||$cols[2] == 93189||$cols[2] == 96262))
        {
            $LP12count = $LP12count + $cols[7];
            print OUT6 $out_line, "\n";
        }
        elsif (($cols[1] == 78099||$cols[1] == 81172||$cols[1] == 84245||$cols[1] == 87318||$cols[1] == 90391||$cols[1] == 93464||$cols[1] == 96537) && $cols[2] == 98745)
        {
            $W2Y1count = $W2Y1count + $cols[7];
            print OUT7 $out_line, "\n";
        }
        elsif (($cols[1] == 77886||$cols[1] == 80959||$cols[1] == 84032||$cols[1] == 87105||$cols[1] == 90178||$cols[1] == 93251||$cols[1] == 96324) && ($cols[2] == 105314))
        {
            $W1BHRF1count = $W1BHRF1count + $cols[7];
            print OUT8 $out_line, "\n";
        }
        elsif (
        ($cols[1] == 140964 && $cols[2] == 152403) ||
        ($cols[1] == 146996 && $cols[2] == 152403) ||
        ($cols[1] == 135024 && $cols[2] == 152589) ||
        ($cols[1] == 140964 && $cols[2] == 152589) ||
        ($cols[1] == 135024 && $cols[2] == 152631) ||
        ($cols[1] == 135277 && $cols[2] == 152631) ||
        ($cols[1] == 140964 && $cols[2] == 152631) ||
        ($cols[1] == 141564 && $cols[2] == 152631) ||
        ($cols[1] == 142030 && $cols[2] == 152631) ||
        ($cols[1] == 142205 && $cols[2] == 152631) ||
        ($cols[1] == 140964 && $cols[2] == 152776) ||
        ($cols[1] == 135024 && $cols[2] == 152838) ||
        ($cols[1] == 135357 && $cols[2] == 152838) ||
        ($cols[1] == 140964 && $cols[2] == 152838) ||
        ($cols[1] == 141979 && $cols[2] == 152838) ||
        ($cols[1] == 142030 && $cols[2] == 152838) ||
        ($cols[1] == 142770 && $cols[2] == 152838)
        )
        {
            print OUT9 "$out_line\t$cols[1]\t$cols[2]\t51,0,102\n";
        }
        elsif (
        ($cols[1] == 146732 && $cols[2] == 152589) ||
        ($cols[1] == 146996 && $cols[2] == 152589) ||
        ($cols[1] == 147815 && $cols[2] == 152589) ||
        ($cols[1] == 150368 && $cols[2] == 152589)
        )
        {
            print OUT9 "$out_line\t$cols[1]\t$cols[2]\t76,0,153\n";
        }
        elsif (
        ($cols[1] == 146732 && $cols[2] == 152631) ||
        ($cols[1] == 146996 && $cols[2] == 152631) ||
        ($cols[1] == 147237 && $cols[2] == 152631) ||
        ($cols[1] == 147815 && $cols[2] == 152631) ||
        ($cols[1] == 147815 && $cols[2] == 152838) ||
        ($cols[1] == 150368 && $cols[2] == 152838)
        )
        {
            print OUT9 "$out_line\t$cols[1]\t$cols[2]\t0,153,0\n";
        }
        elsif 
        ($cols[1] == 146996 && $cols[2] == 152838)
        {
            print OUT9 "$out_line\t$cols[1]\t$cols[2]\t0,102,0\n";
        }
        elsif (
        ($cols[1] == 12367 && $cols[2] == 15626) ||
        ($cols[1] == 40495 && $cols[2] == 41171) ||
        (($cols[1] == 41940) && ($cols[2] == 51368||$cols[2] == 51368||$cols[2] == 57650)) ||
        ($cols[1] == 48326 && $cols[2] == 48433) ||
        ($cols[1] == 48517 && $cols[2] == 50212) ||
        ($cols[1] == 50472 && $cols[2] == 51368) ||
        ($cols[1] == 51507 && $cols[2] == 57650) ||
        ($cols[1] == 56156 && $cols[2] == 57650) ||
        ($cols[1] == 56432 && $cols[2] == 57650) ||
        (($cols[1] == 58046) && ($cols[2] == 63426||$cols[2] == 63728||$cols[2] == 64648||$cols[2] == 64942)) ||
        ($cols[1] == 58945 && $cols[2] == 63426) ||
        (($cols[1] == 63642) && ($cols[2] == 63728||$cols[2] == 64394)) ||
        (($cols[1] == 63828) && ($cols[2] == 63908||$cols[2] == 64394)) ||
        ($cols[1] == 64158 && $cols[2] == 64239) ||
        (($cols[1] == 64321) && ($cols[2] == 64394||$cols[2] == 64648)) ||
        (($cols[1] == 64566) && ($cols[2] == 64648||$cols[2] == 64840)) ||
        (($cols[1] == 64865) && ($cols[2] == 64942||$cols[2] == 68776)) ||
        (($cols[1] == 65052) && ($cols[2] == 66638||$cols[2] == 68776)) ||
        (($cols[1] == 65168) && ($cols[2] == 66638||$cols[2] == 68776)) ||
        ($cols[1] == 65530 && $cols[2] == 68776) ||
        (($cols[1] == 66265) && ($cols[2] == 66638||$cols[2] == 68776)) ||
        ($cols[1] == 66872 && $cols[2] == 68776) ||
        ($cols[1] == 67046 && $cols[2] == 68776) ||
        ($cols[1] == 67212 && $cols[2] == 68776) ||
        ($cols[1] == 67592 && $cols[2] == 68776) ||
        ($cols[1] == 68321 && $cols[2] == 68776)
        )
        {
            print OUT10 $out_line, "\n";
        }
    }

my $output = "W1-W2 depth: $W1W2count\nW2-W1 depth: $W2W1count\nW1-W1 depth: $W1W1count\nW2-W2 depth: $W2W2count\nC2-W1 depth: $C2W1count\nLP 1-2 depth: $LP12count\nW2-Y1 depth: $W2Y1count\nW1-BHRF1 depth: $W1BHRF1count\nBZLF2 and LMP2 junctions in separate file.\n";
    
print $output;
    
print OUT $output;
    
    close(INF);
    close(OUT);
    close(OUT1);
    close(OUT2);
    close(OUT3);
    close(OUT4);
    close(OUT5);
    close(OUT6);
    close(OUT7);
    close(OUT8);
    close(OUT9);
    close(OUT10);
}


