#!/usr/bin/perl -w

=pod

The MIT License (MIT)

Copyright (c) <2015> <Haruka Ozaki>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

=cut

use Time::HiRes qw( time );
use File::Basename;
use POSIX; 
use FindBin;

print STDOUT "MOCCS version 1.3\n";

my $start = time();
my $current;

chomp(@ARGV);

if($#ARGV < 2){ 
    print STDERR "Error: No arguments are specified. Exiting.\n";
    exit;
}


my $k = $ARGV[1]; # length of k-mer 
my $regex = $ARGV[1]; # regex
my $label = $ARGV[2];
my $mask = 0;
my $ofile = "$label.mid.gz";
my $ofile2 = "$label.crf.gz";
my $ofile3 = "$label.auc_count.txt";
my $moccs_out = "$label.MOCCS_visualize.log.txt";
my $script_dir = $FindBin::Bin;
my $rscript = "${script_dir}/MOCCS_visualize.r";
my $fh;
my $is_regex_mode = 0;

if($#ARGV >= 3 && $ARGV[3] eq "-mask"){
    $mask = 1;
}

if($k !~ /^[0-9]+$/){
    $k = length($regex);
    $is_regex_mode = 1;
}

if(! -e $rscript){
    die("Error: $rscript is not found. Exiting.\n");
}

if(! -d dirname($label)){
    mkdir dirname($label)
}


print <<HERE;
[Arguments]
Sequence file: $ARGV[0]
is_regex_mode: $is_regex_mode
k or regex: $regex
Output files: $ofile, $ofile2, $ofile3

HERE



# Preparing hash(kmer) of hash(position) of count (Considering reverse compliment)
my %hash_of_kmer_position_hash = ();
my $rev;
my $mm;

if($is_regex_mode == 1){
    foreach $mm (&nucRegex($regex)){
        $rev = &revComp($mm);
        if(!defined $hash_of_kmer_position_hash{$rev}){
           my %a = ();
           $hash_of_kmer_position_hash{$mm} = \%a;
        }
    }
}else{
    foreach $mm (&nuc($k)){
        $rev = &revComp($mm);
        if(!defined $hash_of_kmer_position_hash{$rev}){
        	my %a = ();
        	$hash_of_kmer_position_hash{$mm} = \%a;
        }
    }
}

print STDOUT "# of ${k}mer: " . scalar(keys %hash_of_kmer_position_hash) . "\n";

$current = time();
printf("Elapsed time: %.2f\n", $current - $start);

# Count occurrence for each position for each kmer one sequence by one.
my $seq_count = 0;
# $/ = ">";
my $end;
my $seq_len;
my $seq;
my ($kmer, $rmer);
my $i;
my $s;

if($ARGV[0] eq "stdin"){
    $fh = STDIN;
}elsif($ARGV[0] =~ /\.gz$/){
    open($fh, "zcat $ARGV[0] |") or die("Error :$!");
}else{
    open($fh, $ARGV[0]) or die("Error :$!");
}


while($s = <$fh>){
    last if $s =~ /^>/;
}
my $header = "";
while(defined($s) && $s =~ /^>/){
    chomp($s);
    # $header = $s;
    $seq = "";
    while(1){
        if(!defined($s = <$fh>) || $s =~ /^>/){
            if($mask == 1){
                $seq =~ tr/acgt/N/;
            }
            $seq = uc($seq);

            # print $seq ."\n";

            for($i = 0; $i <= length($seq)-$k; $i++){
                $kmer = substr($seq, $i, $k);
                $rmer = &revComp($kmer);
                if(defined $hash_of_kmer_position_hash{$kmer}){
                    $hash_of_kmer_position_hash{$kmer}{$i}++;
                }elsif(defined $hash_of_kmer_position_hash{$rmer}){
                    $hash_of_kmer_position_hash{$rmer}{$i}++;
                }
            }

            $seq_len = length($seq);
            $end = $seq_len - $k;

            $seq_count++;
            # if($seq_count % 100000 ==0){
            #     print STDOUT "Finished $seq_count seqs\n";
            # }

            last;
        }
        chomp($s);
        $seq .= $s;
    }
}
close($fh);

print STDOUT "Number of sequences: $seq_count\n";
print STDOUT "Length of sequences: $seq_len\n";

my $cnk = &calcCnk($seq_count, $seq_len, $k);

print "C_nk (low count threshold): " . $cnk . "\n";

$current = time();
printf("Elapsed time: %.2f\n", $current - $start);

# Output Cumulative Relative Frequency
my $center_r = floor(($seq_len - $k + 1) / 2);
my $center_l = ($k % 2 != $seq_len % 2) ? ($center_r - 1) : $center_r;
$center_r -= 1;
$center_l += 1;

open(O, "| gzip > $ofile");
open(O1, "| gzip > $ofile2");
open(O2, "> $ofile3");

print O "kmer\t" . join("\t", (0..$end)) . "\n";
print O1 "kmer\t" . join("\t", (0..$center_l)) . "\n";
print O2 "kmer\tauc\tcount\n";

my $c_skip = 0;
foreach my $key (sort keys %hash_of_kmer_position_hash){
    my @res = ();
    my @crf = ();
    my $sum = 0;
    foreach my $pos (0..$end){
        if(defined $hash_of_kmer_position_hash{$key}{$pos}){
            push(@res, $hash_of_kmer_position_hash{$key}{$pos});
            $sum += $hash_of_kmer_position_hash{$key}{$pos};
        }else{
            push(@res, 0);
        }
    }
    if($sum < $cnk){
        $c_skip++;
        next;
    }

    push(@crf, 0);

    if($center_l - 1 == $center_r + 1){
        push(@crf, $res[$center_l-1]);
    }else{
        push(@crf, $res[$center_l - 1] + $res[$center_r + 1]);
    }

    for($i = 2; $i <= $center_l; $i++){
        push(@crf, $res[$center_l - $i] + $res[$center_r + $i] + $crf[$i - 1]);
    }

    if($sum != $crf[$#crf]){
        print STDERR "error: sum $sum != crf $crf[$#crf]\n";
    }

    my $auc = 0;
    for($i = 0; $i <= $#crf; $i++){
        $crf[$i] = $crf[$i] / $sum;
        $auc += $crf[$i] - (1/$center_l) * $i;
    }

    print O $key . "\t" . join("\t", @res) . "\n";
    print O1 $key . "\t" . join("\t", @crf) . "\n";
    print O2 $key . "\t" . $auc . "\t" . $sum . "\n";
}

close(O);
close(O1);
close(O2);

print "# of skipped kmer: $c_skip\n";

$current = time();
printf("Elapsed time: %.2f\n", $current - $start);


# Visualize cumulative relative frequency curves
my $rcmd = "Rscript --no-save --no-restore \"$rscript\" 'file=\"$ofile2\"' 'file2=\"$ofile3\"' 'label=\"$label\"' 'len=$seq_len' 'k=$k' > $moccs_out";

print STDOUT "\n" . "Command to visualize MOCCS result: " . $rcmd . "\n\n";
system($rcmd);

$current = time();
printf("Elapsed time: %.2f\n", $current - $start);


sub nuc(){
    my ($width) = @_;
    my @base = ("A", "T", "G", "C");
    if($width == 1){
    	return @base;
    }else{
	   my @new = ();
    	foreach my $mer (&nuc($width-1)){
    	    foreach my $b (@base){
    		push(@new, $mer . $b);
    	    }
    	}
	   return @new;
    }
}

sub revComp(){
    my ($sequence) = @_;
    chomp($sequence);
    my $reverse_complement = reverse $sequence;
    $reverse_complement =~ tr/ATGCatgc/TACGtacg/;
    return $reverse_complement;
}

sub nucRegex(){
    my ($regex) = @_;
    # my @base = ("A", "T", "G", "C");
    my @new = ();
    my $mer;
    my $b;
    if(length($regex) == 1){
        foreach $b (&regexBases($regex)){
            push(@new, $b);
        }
    }else{
        foreach $mer (&nucRegex( substr($regex, 0, -1) )){
            foreach $b (&regexBases($regex)){
                push(@new, $mer . $b);
            }
        }
    }
    return @new;
}

sub regexBases(){
    my ($regex) = @_;
    $regex = substr($regex, -1);
    if($regex eq "N"){
        return ("A", "C", "G", "T");
    }else{
        return ($regex);
    }
}

sub calcCnk(){
    my ($N, $w, $k) = @_;
    return $N * ($w - $k + 1) / (4 ** $k)
}


exit;
