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

use Getopt::Long qw(:config posix_default no_ignore_case gnu_compat);
use Time::HiRes qw( time );
use File::Basename;
use POSIX; 
use FindBin;

print STDOUT "MOCCS version 1.7\n";

my $start = time();
my $current;

my $input;
my $k;
my $regex; # regex
my $label = "MOCCS_result";
my $is_regex_mode = 0;
my $mask = 0;
my $use_threshold = 0;
my $threshold;
my $script_dir = $FindBin::Bin;
my $rscript = "${script_dir}/MOCCS_visualize.r";
my $fh;
my $stranded = 0;
my $cnk;

GetOptions(
    'i=s' => \$input,
    'k=i' => \$k,
    'regex=s' => \$regex,
    'label=s' => \$label,
    'mask' => \$mask,
    'threshold=f' => \$threshold,
    'stranded' => \$stranded,
    'low-count-threshold=f' => \$cnk
    );

if($input eq ''){
    print STDERR "Input FASTA name is not specified. Exiting.\n";
    exit;   
}

my $ofile1 = "$label.mid.gz";
my $ofile2 = "$label.crf.gz";
my $ofile3 = "$label.auc_count.txt";
my $moccs_out = "$label.MOCCS_visualize.log.txt";

print <<HERE;
[Arguments]
Sequence file: $input
Output files: $ofile1, $ofile2, $ofile3
HERE

if(defined $k){
    if(defined $regex){
        print STDERR "-k and --regex cannot be used at the same time. Exiting.\n";
        exit;
    }
    print "k: $k\n";
}elsif(defined $regex){
    $k = length($regex);
    $is_regex_mode = 1;
    print "regex: $regex\n";
}else{
    print STDERR "-k or --regex have to be specified. Exiting.\n";
    exit;
}
if(defined $threshold){
    $use_threshold = 1;
    print "threshold: $threshold\n";
}

if($mask == 1){
    print "mask: true\n";
}

if($stranded == 1){
    print "stranded: true\n";
}

if(! -e $rscript){
    die("Error: $rscript is not found. Exiting.\n");
}

if(! -d dirname($label)){
    mkdir dirname($label)
}

if(defined $threshold){
    $use_threshold = 1;
    print "threshold: $threshold\n";
}


# Preparing hash(kmer) of hash(position) of count (Considering reverse compliment)
my %hash_of_kmer_position_hash = ();
my $rev;
my $mm;

if($is_regex_mode == 1){
    foreach $mm (&nucRegex($regex)){
        $rev = &revComp($mm);
        if($stranded == 1){
           my %a = ();
           $hash_of_kmer_position_hash{$mm} = \%a;
           my %b = ();
           $hash_of_kmer_position_hash{$rev} = \%b;
        }else{
            if(!defined $hash_of_kmer_position_hash{$rev}){
               my %a = ();
               $hash_of_kmer_position_hash{$mm} = \%;
            }
        }
    }
}else{
    foreach $mm (&nuc($k)){
        $rev = &revComp($mm);
        if($stranded == 1){
           my %a = ();
           $hash_of_kmer_position_hash{$mm} = \%a;
           my %b = ();
           $hash_of_kmer_position_hash{$rev} = \%b;
        }else{
            if(!defined $hash_of_kmer_position_hash{$rev}){
               my %a = ();
               $hash_of_kmer_position_hash{$mm} = \%a;
            }
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

if($input eq "stdin"){
    $fh = STDIN;
}elsif($input =~ /\.gz$/){
    open($fh, "gunzip -c $input |") or die("Error :$!");
}else{
    open($fh, $input) or die("Error :$!");
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
            $seq_len = length($seq);

            for($i = 0; $i <= $seq_len-$k; $i++){
                $kmer = substr($seq, $i, $k);
                if(defined $hash_of_kmer_position_hash{$kmer}){
                    $hash_of_kmer_position_hash{$kmer}{$i}++;
                }else{
                    if($stranded == 0){
                        $rmer = &revComp($kmer);
                        if(defined $hash_of_kmer_position_hash{$rmer}){
                            $hash_of_kmer_position_hash{$rmer}{$i}++;
                        }
                    }
                }
            }

            
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

if(! defined $cnk){
    $cnk = &calcCnk($seq_count, $seq_len, $k);
}
print sprintf("C_nk (low count threshold): %.1f\n", $cnk) ;

$current = time();
printf("Elapsed time: %.2f\n", $current - $start);

# Output Cumulative Relative Frequency
my $center_r = floor(($seq_len - $k + 1) / 2);
my $center_l = ($k % 2 != $seq_len % 2) ? ($center_r - 1) : $center_r;
$center_r -= 1;
$center_l += 1;

open(O1, "| gzip > $ofile1");
open(O2, "| gzip > $ofile2");
open(O3, "> $ofile3");

print O1 "kmer\t" . join("\t", (0..$end)) . "\n";
print O2 "kmer\t" . join("\t", (0..$center_l)) . "\n";
print O3 "kmer\tauc\tcount\n";

my $c_skip = 0;
my @res = ();
my @crf = ();
my $sum = 0;
my $auc = 0;
foreach my $key (sort keys %hash_of_kmer_position_hash){
    @res = ();
    @crf = ();
    $sum = 0;
    foreach my $pos (0..$end){
        if(defined $hash_of_kmer_position_hash{$key}{$pos}){
            # push(@res, $hash_of_kmer_position_hash{$key}{$pos});
            $sum += $hash_of_kmer_position_hash{$key}{$pos};
        }
    }
    if($sum < $cnk){
        $c_skip++;
        next;
    }

    foreach my $pos (0..$end){
        if(defined $hash_of_kmer_position_hash{$key}{$pos}){
            push(@res, $hash_of_kmer_position_hash{$key}{$pos});
        }else{
            push(@res, 0);
        }
    }

    push(@crf, 0);
    $auc = 0;

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

    
    for($i = 0; $i <= $#crf; $i++){
        $auc += $crf[$i] / $sum - (1/$center_l) * $i;
    }

    if($use_threshold == 1 && $auc <= $threshold){
        $c_skip++;
        next;
    }
    print O1 $key . "\t" . join("\t", @res) . "\n";
    print O2 $key . "\t" . join("\t", @crf) . "\n";
    print O3 $key . "\t" . $auc . "\t" . $sum . "\n";
}

close(O1);
close(O2);
close(O3);

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
