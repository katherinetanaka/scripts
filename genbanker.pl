#! /Users/Katherine/perl5/perlbrew/perls/perl-5.16.2/bin/perl

use strict;
use warnings;
use Getopt::Long;

############################################################
# A Genbanker for visualization in Artemis, for instance...#
############################################################

# Global
my $fileopt;
my $type;
my $color = "no";
#####################
#### Subroutines ####
#####################

sub genbankWriter{
# Retrieve and define variables
    my $file = shift;
    my $filename;
    for($file =~ m/(.*\/)*(.*)/){
        $filename = join(".", ($2, "gb")); 
    }
    my $type = shift;
    my $color = shift;
    
    my $key;
    my $awk;
    
    if ($type eq "IS") {
        $key = "mobile_element";
        $awk = "awk '(\$3 >= 60 || \$3 ~ /(100.00)/) && \$4 >= 30 {print}'";
    } elsif ($type eq "read") {
        $key = "misc_feature";
        $awk = "awk '(\$3 >= 90 || \$3 ~ /(100.00)/) && \$4 >= 10 {print}'";
    } else {
        die ("genbankWriter usage: genbankWriter('file', 'type') with type in 'IS', 'read'.\n");
    }
# Sort and filter blast results, write genbank file   
    open(my $genbank, ">", $filename) || die("Could not create a new file :(\n");
    open(my $blastres, '-|', "cat $file | $awk | sort -n -k9,9") || die("Invalid file or could not be opened.\n");
# Begin printing    
    
#    print $genbank "FEATURES", "\t\t\t", "Locations/Qualifiers\n";
    while(my $line = <$blastres>){
        print $line;
        $line =~ m/^(.+)\s(.+)\s(.+)\s(\d+)\s(\d+)\s(\d+)\s(\d+)\s(\d+)\s(\d+)\s(\d+)\s(.+)\s(.+)$/;
        # Direct or reverse complement????
        
        if ($key eq "mobile_element") {
            if ($7 < $8) {
                print $genbank "     ", "$key", "  ", "$9..$10\n";
            } else{
                print $genbank "     ", "$key", "  ", "complement($9..$10)\n";
            }
            print $genbank "                     ", "/mobile_element_type=\"$1\"\n";
            print $genbank "                     ", "/note=\"Identity percentage=$3\"\n";
        unless($color eq "no"){
            print $genbank "                     ", "/color=$color\n";
        }
        } else {
            if ($7 < $8) {
                print $genbank "     ", "$key", "    ", "$9..$10\n";
            } else{
                print $genbank "     ", "$key", "    ", "complement($9..$10)\n";
            }
            print $genbank "                     ", "/note=\"cluster=$1, percentage=$3\"\n";
        unless($color eq "no"){
            print $genbank "                     ", "/color=$color\n";
        } 
    }
    }
#    print $genbank "//";
    close($blastres);
    close($genbank);
}

###############
#### Main #####
###############

GetOptions('file=s' => \$fileopt, 'type=s' => \$type, 'color=s' => \$color) || die("./samread.pl -file [FILE] -cores [CORES]\n");
genbankWriter($fileopt, $type, $color); # write genbank

print "End of method\n";