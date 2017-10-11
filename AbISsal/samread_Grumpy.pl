#! /usr/bin/perl

=head1
CALL: ./samread_Grumpy.pl -file [FILE] -cores [CORES]
a replacement to samread2.pl
part of the parallel Blast fork
todo:   extract left-anchored and right-anchored reads
        rev comp left-anchored reads (so its anchored at 1)
        sort from longest to shortest
        hash by hierarchy: name and ref seq main
Katherine Tanaka, July 27th, 2015
=cut

use strict;
use Benchmark;
use Getopt::Long;
use Data::Dumper;
use Parallel::ForkManager;

# Some variables have to stay global!
#### Name counting variables ####
my $samfile;
my $cores;
my @names;
my %ambi;

my @twos; # Finally I will make them arrays, but more probably for future uses
my @threes;
my @fours;
my @fives;
my @append_name;
#################################

##### File parsing variables ####
my $samfile;
my %LengthHeader;
my %Lefties;
my %Rightoes;

my $leftcounter = 0;
my $rightcounter = 0;
my $badcounter = 0;
my $fullcounter = 0;
my $nocounter = 0;
my $miscounter = 0;
################################


######################
#### Benchmarking ####
######################
=head3
what are times?
t0 = beginning, now after variables declaration
t1 = pairing
t2 = after parsing and left-right creating, before preparation for blast
t3 = after parallelBlast methods
t4 = end
=cut
my $t0 = Benchmark->new;
my $t1;
my $t2;
my $t3;
my $t4;

#####################
#### Subroutines ####
#####################

sub printableTime{
    my @times = localtime(time);
    @times = reverse(@times[0..5]);
    my $date = join("-", @times[0..2]);
    my $hour = join("m", @times[3..5]);
    my $dirname = join("_", ($hour, $date));
    return $dirname;
}
sub revcompAll {
# Converter key == drop a dependance!
my %converter = (   A => "T",
                    T => "A",
                    G => "C",
                    C => "G");
# Core subroutine
    die("Bad subroutine revcompAll\n") unless @_;
    my %options = %{ shift @_ };
    my @sequences = @{$options{"sequence"}};
    my @newseq;
    foreach my $seq (@sequences){
        my @decompose = split(//, $seq);
        my @reverse = reverse(@decompose);
        my @revcomp;
        foreach my $base (@reverse) {
            push(@revcomp, $converter{$base});
        }
        push(@newseq, join("", @revcomp));
    }
    @{$options{"sequence"}} = @newseq;
    return(\%options);
}
sub cutter{
    die("Bad subroutine cutter\n") unless @_;
    my %options = %{ shift @_ };
    my @sequences = @{$options{"sequence"}};
    my @cutter = @{$options{"truncLength"}};
    my @newseq;
    my $i = 0;
    
        while($i < scalar(@sequences)){
            push(@newseq, substr($sequences[$i], -$cutter[$i]));
            $i ++;
        }
    @{$options{"sequence"}} = @newseq;
    return(\%options);
}
sub thesorter {
    die("Bad subroutine\n") unless @_;
    my %options = %{ shift @_ };
    #print Dumper \%options;
    my ($sorting, $sorted1, $sorted2, $sorted3) = (@options{"truncLength",
                                                            "read", "refseq",
                                                            "sequence"});
    my @derefsort = @$sorting;
    
    my @sorted_index = sort { $derefsort[$b] <=> $derefsort[$a] } 0 .. $#derefsort;
    my @lengthS = @$sorting[@sorted_index];
    my @readS = @$sorted1[@sorted_index];
    my @refS = @$sorted2[@sorted_index];
    my @seqS = @$sorted3[@sorted_index];
  
    @{$options{"truncLength"}} = @lengthS;
    @{$options{"read"}} = @readS;
    @{$options{"refseq"}} = @refS;
    @{$options{"sequence"}} = @seqS;
    return (\%options);

}
sub parallelBlast{
# Define variables
    my $directory = shift;
    my $side = shift;
    my $dbname = shift;
    my $cores = shift;
    my %hashsub = %{ shift @_ };
    
    my $manager = new Parallel::ForkManager($cores);
    my $stopindex = keys(%hashsub)-1;
# Begin parallel task
print "...Begin parallel blast -- $side --- in $directory...\n";    
   for(my $m = 0; $m <= $stopindex; $m ++){
    
    if ($m % 100 == 0 && $m != 0){
        $manager->wait_all_children;
        print $m, " reads done to the $side. Reevaluating task to be done\n";
        open(my $apipe, '-|', " cat $directory/*_$side | awk '\$1 !~ /#/ && \$3 >= 95 && \$7 <= 3 && \$9 <= 3 && \$4 >= \$14-3 {print \$2}' | sort -u");
        while (my $line = <$apipe>){
            chomp $line;
       #     print "\t deleting $line\n";
            delete($hashsub{$line});
        }
        close($apipe); 
        print "More reads deleted.\n";        
    }
    # Begin a task
    
     if (exists($hashsub{$m}) && length($hashsub{$m})>=30){
     $manager->start and next;
         open(my $pipe, "| blastn -db $dbname -out $directory/$m\_$side -outfmt '7 std qlen slen' -max_target_seqs $stopindex");
             print $pipe ">$m\n", $hashsub{$m}, "\n";
     close($pipe);
     $manager->finish;
    } 
    
}
$manager->wait_all_children;
return 0;
}

sub orphan{
  my @orphans;
  my $ref_all = shift;
  my $ref_counted = shift;
  my @all = @$ref_all;
  my @counted = @$ref_counted;
  
  my %hash_all     = map { $all[$_] => $all[$_] }         0..$#all;
  my %hash_counted = map { $counted[$_] => $counted[$_] } 0..$#counted;
  
  for my $keys (keys(%hash_counted)){
    unless(defined $hash_all{$keys}){
        push(@orphans, $keys);
    }
  }
  return @orphans;
}

sub thebin {
    my $dir  = shift;
    my $side = shift;
    my $max  = shift;
    my %hashQuery = %{ shift @_ };
    my @printer;

#### Bining part ####    
open(my $results, "-|", "cat $dir/*_$side |sort -n -k1,1 -k2,2|sed -e '/#/d'");
my $n = 0;
my $previous = -1;
while(my $line = <$results>){
$line =~ m/^(\d+)\s(\d+)\s(.+)\s(\d+)\s(\d+)\s(\d+)\s(\d+)\s(\d+)\s(\d+)\s(\d+)\s(.+)\s(.+)\s(\d+)\s(\d+)$/;
    ######## Redundancy skip #########
    if( (grep {/^$1$/} @printer) && (grep {/^$2$/} @printer) ){
        $previous = $1;
        next;
    } elsif($1 != $previous && $1 == $2){  # First of number, first match itself == bin
    #    print "In -$side: New bin is created! $n bin! Bin of the $1! ", scalar(keys(%hashQuery)), " are there!\n";
        #### Bin creator
        $hashQuery{$1}{bin} = $1;
        $hashQuery{"Chevron_$n"}{query} = ">Bin $n";
        push(@printer, "Chevron_$n");
        $n ++;
        push(@printer, $1);
        ##### End of bin creator                        # Get next bin done
    } elsif ($1 != $previous && $1 != $2){  # Orphan testing
        push(my @quest, $1);
        my @isthere = orphan(\@printer, \@quest);
        if(@isthere != 0){
            print "$1 is an orphan. It should have a bin on its own :) So...\n";
            print "In -$side: New bin is created! $n bin! Bin of the $1! ", scalar(keys(%hashQuery)), " are there!\n";
        #### Bin creator
            $hashQuery{$1}{bin} = $1;
            $hashQuery{"Chevron_$n"}{query} = ">Bin $n";
            push(@printer, "Chevron_$n");
            $n ++;
            push(@printer, $1);
        ##### End of bin creator
       }
       
    } else{                                  # Which should be $1 == $previous
      #  print $2, "\t";
        push(my @otherquest, $2);
        my @search = orphan(\@printer, \@otherquest);
       # print @search, "\n";
        if(@search != 0){     # Not indexed, but no new bin created, append $1
            if($3 >= 95 && $7 <= 3 && $9 <= 3 && $4 >= ($14-3)){
                $hashQuery{$2}{bin} = $1;
                push(@printer, $2);
            }
        }
    }
    $previous = $1;
    
}
close($results);

#### Searching for orphans part via $max ####
my @pos_children = (0..$max);
my @orphans = orphan(\@printer, \@pos_children);
@orphans = sort { $a <=> $b } @orphans;

if(@orphans > 0){
    print "There (are) orphan(s):\n";
    print join("\t", @orphans), "\n";
} else {
    print "No orphan in there. \n";
}
#### End of orphans searching ####
    print "End of bining for -- $side --\n";
    return(\@printer, \@orphans, \%hashQuery);

}

############# End of subroutines ##############

##############
#### Main ####
##############

#### File I/O and arguments ####
GetOptions('file=s' => \$samfile, 'cores=s' => \$cores) || die("./samread.pl -file [FILE] -cores [CORES]\n");



######## Name counting for pairing informations #######
# Pipe read names for tagging
print ".....Begining piping.....\n";
open my $pipe,  "-|", "awk '\$1 !~ /@/ {print \$1}' $samfile";
while (my $line = <$pipe>){
    chomp $line;
    push(@names, $line);
    }
    close($pipe);
print ".....End piping.....\n";



print ".....Begining counting.....\n";

# How many for each category?
my $i;
while($i < scalar(@names)){
    if ($names[$i] eq $names[$i+4]){
    # Fives!
        push(@fives, $names[$i]);
        push(@append_name, ("1:5", "2:5", "3:5", "4:5", "5:5"));
        $i += 5;
    # Fours
    } elsif ($names[$i] eq $names[$i+3]){
        push(@fours, $names[$i]);
        push(@append_name, ("1:4", "2:4", "3:4", "4:4"));
        $i += 4;
    # Threes    
    } elsif ($names[$i] eq $names[$i+2]){
        push(@threes, $names[$i]);
        push(@append_name, ("1:3", "2:3", "3:3"));
        $i += 3;
    # Twos
    } elsif ($names[$i] eq $names[$i+1]){
        push(@twos, $names[$i]);
        push(@append_name, ("1:2", "2:2"));
        $i += 2;
    # Safety
    } else {
        print $i, $names[$i], "\n";
        die("Bad sam\n");
    }
}
my $total = @twos + @threes + @fours + @fives;
print "---------- Finished getting paired-end informations --------------\n",
"Total pairs = $total\n",
"Two alignments a pair    =  ", scalar(@twos), " (", sprintf("%.3f", @twos/$total*100), "%)\n",
"Three alignments a pair  = ", scalar(@threes), " (", sprintf("%.3f", @threes/$total*100), "%)\n",
"Four alignments a pair   = ", scalar(@fours), " (", sprintf("%.3f", @fours/$total*100), "%)\n",
"Five alignments a pair   = ", scalar(@fives), " (", sprintf("%.3f", @fives/$total*100), "%)\n",
"------------------------------------------------------------------\n\n",
"....Now extracting overhang reads from file...\n\n";

# Benchmark
$t1 = Benchmark->new;
####################################################

#### Line parsing ####
open(my $sam, "<", $samfile) || die("Samfile does not exist!\n");

while (my $line = <$sam>) {
# Now go on the main parsing
    if (substr($line, 0, 3) eq '@SQ'){
    # Is a header
        my @headerLine = split("\t", $line);
        my $sequence = substr($headerLine[1], 3);
        my $seqLength = substr($headerLine[2], 3);
        $LengthHeader{$sequence} = $seqLength;
    } elsif (substr($line, 0, 3) eq '@PG'){
        print $line, "\n"; # Is basic infos
    } else {
        # Is an alignment
        my @aln = split(/\t/, $line);
        my ($qname, $refname, $leftpos, $cigar, $fullseq) = @aln[0,2,3,5,9];
        my $qname2 = join("_", ($qname, shift(@append_name)));
       # print $qname, "\n";     
        # Is it a leftie?
        if($leftpos == 1 & $cigar =~ m/^(\d{1,3})(S)(\d{1,3}M[\dDMI]*)$/){
            # print "Leftie!\t", $line, "\n";
            push(@{$Lefties{"truncLength"}}, $1);
            push(@{$Lefties{"read"}}, $qname2);
            push(@{$Lefties{"refseq"}}, $refname);
            push(@{$Lefties{"sequence"}}, $fullseq);
        # Add 1 to leftie
        $leftcounter ++; # End of left ###################################################
        } elsif($cigar =~ m/^([\dDMI]*\d{1,3}M)(\d{1,3})(S)$/){
        # Is it a rightie?
        # In need of further indications!
            my $Sonly = $2;
            my $matched = length($fullseq) - $Sonly;
            #print $Sonly, "\n";
            my $length = $LengthHeader{$refname};
            if($length - $matched + 1 == $leftpos){
              #  print "Righto!\t", $line, "\n";
                push(@{$Rightoes{"truncLength"}}, $Sonly);
                push(@{$Rightoes{"read"}}, $qname2);
                push(@{$Rightoes{"refseq"}}, $refname);
                push(@{$Rightoes{"sequence"}}, $fullseq);
        # Add 1 to rightoes
                $rightcounter ++; # End of right #########################################
            } else {
        # Misc: hanging read in the middle
                $miscounter ++;
            } 
        } elsif ($cigar =~ m/^([\dDMI]*\d{1,3}M)/){
            $fullcounter ++;
        } elsif ($cigar =~ m/^([\dDMIS]*\d{1,3}M)(\d{1,3})(S)$/){
            $badcounter ++;
        } elsif ($cigar eq "*") {
            $nocounter ++;
        } else {
            $miscounter ++;
        }
}
}
close($sam);

#### Modifications ####
# Lefties need to be reverse complemented, then cropped, then ordered
revcompAll(\%Lefties);
cutter(\%Lefties);
thesorter(\%Lefties);
# Rightoes only needed to be cutted and sorted
cutter(\%Rightoes);
thesorter(\%Rightoes);

####### Stats ########
my $allcounter = $leftcounter + $rightcounter + $badcounter +
                 $fullcounter + $nocounter + $miscounter;

print "---------- Finished parsing your complete samfile --------------\n",
"Total reads = $allcounter\n",
"Not aligned reads    = $nocounter (", sprintf("%.3f", $nocounter/$allcounter*100), "%)\n",
"Fully aligned reads  = $fullcounter (", sprintf("%.3f", $fullcounter/$allcounter*100), "%)\n",
"Miss aligned reads   = ", $badcounter + $miscounter, " (", sprintf("%.3f", ($badcounter + $miscounter)/$allcounter*100), "%)\n",
"Left overhang reads  = $leftcounter (", sprintf("%.3f", $leftcounter/$allcounter*100), "%)\n",
"Right overhang reads = $rightcounter (", sprintf("%.3f", $rightcounter/$allcounter*100), "%)\n",
"-----------------------------------------------------------------\n\n",
"...Preparing databases for blastn...\n";

# Benchmark
$t2 = Benchmark->new;

###### Database and query preparation ########
my %leftFasta = map {$_ => { query => @{$Lefties{read}}[$_],
                             sequence => @{$Lefties{sequence}}[$_],
                             reference => @{$Lefties{refseq}}[$_] }  } 0..$#{$Lefties{read}};
my %rightFasta = map {$_ => { query => @{$Rightoes{read}}[$_],
                             sequence => @{$Rightoes{sequence}}[$_],
                             reference => @{$Rightoes{refseq}}[$_] }  } 0..$#{$Rightoes{read}};

# Database writing
print "...Creating .fna files for left and right...\n";

my $leftname = "left_abISsal.fna"; ##########
open(my $dbleft, ">", $leftname) || die("Could not create file ");
foreach my $j (0..keys(%leftFasta)-1){
    print $dbleft ">", $j, "\n",
                  $leftFasta{$j}{sequence}, "\n";
}
close($dbleft) || die ("Could not close dbleft\n");

my $rightname = "right_abISsal.fna"; ##########
open(my $dbright, ">", $rightname) || die("Could not create file ");
foreach my $k (0..keys(%rightFasta)-1){
    print $dbright ">", $k, "\n",
                  $rightFasta{$k}{sequence}, "\n";
}
close($dbright) || die ("Could not close dbleft\n");

# Database indexing
print "...Indexing databases...\n";
`makeblastdb -in $leftname -dbtype nucl > /dev/null`;  # So that messages are hidden
`makeblastdb -in $rightname -dbtype nucl > /dev/null`;

# Make expendable hashes for left and right
my %expL = map{ $_ => $leftFasta{$_}{sequence} } 0..keys(%leftFasta)-1;
my %expR = map{ $_ => $rightFasta{$_}{sequence} } 0..keys(%rightFasta)-1;


#########################
#### Parallelization ####
#########################

my $dirname = printableTime();
mkdir($dirname) || die("could not make dir.\n");

# parallelBlast($directory, "left", $leftname, $cores, \%hash);
parallelBlast($dirname, "left", $leftname, $cores, \%expL);
parallelBlast($dirname, "right", $rightname, $cores, \%expR);

print "Blast successful!\n";

# Benchmark
$t3 = Benchmark->new;

=head2
September 9th, 2015.
Now is maybe the time to merge samread_Grumpy.pl and bining_test.pl
Main points to check while merging:
	- references to hashes: are they really %leftFasta and %rightFasta?? #
	- proper ref attribution before subroutines #
	- proper variables calling in imported subroutines: global vs local #
	- get maxima array from scalar(keys(%dirFasta)) #
	- directories changing, data writing with $dirname, need to exist as argument of thebin #
	- create a new manager of only 2 forks mandatory for left and right analysis #
	- make new time
=cut

#########################
####   Bin making    ####
#########################

# Imported variables from bining_test.pl, with a twist.
my @sides = ("left", "right");
my @maxima = (scalar(keys(%leftFasta)), scalar(keys(%rightFasta)));
my @ref_to_hashes = (\%leftFasta, \%rightFasta);
my $newManager = new Parallel::ForkManager(2);

# Main bin making
for (my $n = 0; $n <=1; $n ++) {

    $newManager->start and next;
      my ($refP, $refO, $refH) = thebin($dirname, $sides[$n], $maxima[$n], $ref_to_hashes[$n]);
      my @printer = @$refP;
      my @orphans = @$refO;
      my %hash_print = %$refH;
      my $filename = join("", ("$dirname/", "abISsal_", $sides[$n], "report"));
      open(my $Report, ">", $filename) || die("Unable to create $filename!\n");
      foreach my $tag (@printer){
            print $Report join("\t", ($hash_print{$tag}{query},
                                  $hash_print{$tag}{reference},
                                  $hash_print{$tag}{bin},
                                  $hash_print{$tag}{sequence}), "\n");
      }
      print $Report "There are orphans too!\n", join("\t", @orphans), "\n", "END";
      close($Report);
    $newManager->finish;
}
$newManager->wait_all_children;


######################
#### Benchmarking ####
######################
$t4 = Benchmark->new;
my $pairing  = timediff($t1, $t0);
my $siding   = timediff($t2, $t1);
my $blasting = timediff($t3, $t2);
my $bining   = timediff($t4, $t3);
my $total    = timediff($t4, $t0);
print "---------- AbISsal -09/09/15- has finished analyzing your reads --------------\n",
"Total time = ", timestr($total), "\n",
"Retrieving pairing informations took", timestr($pairing), "\n",
"Assiging relevant reads to 5' or 3' end took", timestr($siding), "\n",
"Blasting 5' and 3' reads together using $cores cores took", timestr($blasting), "\n",
"Bining reads together based on blast results, in 2 cores, took", timestr($bining), "\n",
"Your results are in $dirname/abISsal_#side#report, and raw blast are still avaliable\n",
"-----------------------------------------------------------------\n\n",
"...We hope to see you again!...\n";
