#! /usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;

my $fileopt;
my %master;
my %master2; 
my @secondkeys = ("ref", "strand", "score", "id", "tstart", "tstop", "tlength",
                      "qstart", "qstop", "qlength");
my %report;

# 2016/01/14: New parsing format to accomodate tabulate version, which is enough for
#  the application I need

# Get tabulated file
GetOptions('file=s', \$fileopt);
open(my $file, "<", $fileopt) || die("Could not find such a file.\n");

# First thing is to sort relevant information in a hash
while (my $line = <$file>) {
    chomp $line;
    my @fields = split(/\s/, $line);
    
# Create name for key: remove movie name
    my @subnames = split(/\//, $fields[0]);
    my $key = join("/", $subnames[1], $subnames[2], $subnames[3]);
    
# Store: push is the key to keep multiple alignment under the same name
    
    my @updatedF = @fields[1,3,4,5,6,7,8,9,10,11];
    for (my $i = 0; $i < scalar(@updatedF); $i ++) {
        push(@{$master{$key}{$secondkeys[$i]}}, $updatedF[$i]);
    }
}
close($file);

=head1
Conversion: From multiple parameters under the same read name to multiple alignments
under the same read name
 Nested Map are used to iterate over the parameter directly outside the bracket: together,
 they build an array for each of the previous array index, with the index as a key 
=cut
for my $keys3 (keys(%master)) {
    %{$master2{$keys3}} = map { 
        my $num = $_;
        $num => [ map { $master{$keys3}{$_}[$num] } @secondkeys ]
        }     0..scalar(@{$master{$keys3}{"score"}})-1;
}

# Inversion: alignment arrays on strand 1 need to have their start/stop switched and
#  positions needs to be recalculated (against tlength)
for my $keys4 (keys(%master2)) {
    for my $keysPrime (keys($master2{$keys4})) {
        if ($master2{$keys4}{$keysPrime}[1] eq '1') {
          my $refarray = \@{$master2{$keys4}{$keysPrime}};
          #  print join("\t", @{$master2{$keys4}{$keysPrime}}[1,4,5]), "\t";
          my $nStart = @$refarray[6] - @$refarray[5];
          my $nStop  = @$refarray[6] - @$refarray[4];
          @$refarray[4,5] = ($nStart, $nStop);
          # print join("\t", @{$master2{$keys4}{$keysPrime}}[1,4,5]), "\n";
        }
    } 
}

# Filter reads to keep only partials and split aligned
  
# Getting the loop done
for my $keys5 (keys(%master2)) {
     my $refAlign = \@{$master2{$keys5}{"0"}};  # 0 should be main alignment
    my $outerStart = @$refAlign[4] - (@$refAlign[5] - @$refAlign[4]);
    my $outerStop = @$refAlign[5] + (@$refAlign[5] - @$refAlign[4]);
    my $length = @$refAlign[9] - (@$refAlign[8] - @$refAlign[7]);  # Trailing length on either side
    
    # Setting 0 first
    if ($length <= 20) {  # Remove completely aligned
        $report{$keys5}{"0"} = "complete";   
    } else {
        $report{$keys5}{"0"} = "partial";
    }
    # Setting the secondary alignments
    for ( my $j = 1; $j < keys($master2{$keys5}); $j ++) {
                my $refSub = \@{$master2{$keys5}{$j}};
                
        if (@$refSub[0] eq @$refAlign[0] && (@$refSub[4] >= $outerStart && $outerStop >= @$refSub[4])) {  
# Second alignment centroid overlap with main: probably adapter non detected
                $report{$keys5}{$j} = "misadapter";
        } elsif (@$refSub[7] >= @$refAlign[7] && @$refSub[7] <= @$refAlign[8]) {
# Secondary alignment cover same query nt than main: aligned over repeat region
                $report{$keys5}{$j} = "repeat";
        } else {
                $report{$keys5}{$j} = "discordant";
        }
    }
}

# Express results
for my $keys6 (keys(%report)) {
    for ( my $k = 1; $k < keys($report{$keys6}); $k ++) {
        if ($report{$keys6}{$k} eq "discordant") {
           if(@{$master2{$keys6}{"0"}}[0] eq "pAsa5_propre" && (@{$master2{$keys6}{$k}}[0] eq "pAsa5_propre")) {
                print $keys6, "\n0:", join("\t", @{$master2{$keys6}{"0"}}), "\t", 
                $report{$keys6}{"0"}, "\n", $k, ":", join("\t", @{$master2{$keys6}{$k}}), "\n";
            }
        }
    }
}



# print Dumper \%report;
print "All done now!\n";
