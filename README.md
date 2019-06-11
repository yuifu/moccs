# MOCCS

Motif Centrality Analysis of ChIP-Seq (MOCCS) is a method for for clarifying DNA-binding motif ambiguity.

Given ChIP-Seq data of any DNA-binding proteins including transcription factors, MOCCS comprehensively analyzes and describes every $k$-mer that is bound by the DNA-binding proteins.

MOCCS (version 1.7) is written in Perl and R.  
MOCCS was tested on Perl version 5.18.2 and R version 3.2.1.


## Version history

- 2019/06/11: Version 2.0.
	- Added MOCCS2 score in the output (See below for MOCCS2 score).
- 2016/04/15: Version 1.7.
	- Added 'stranded' option. With this optin, MOCCS will count k-mers on the forward strand and will not reverse-complement k-mers. This is useful when analyzing RNA-binding proteins (RBPs) data (e.g. CLIP-Seq data).
	- Added '--low-count-threshold' option.

## Usage

	perl MOCCS.pl -i <string> -k <int>|--regex <string> \
		[--label <string>] [--mask] [--threshold <float>]

- `-i <string>`: FASTA file containing fixed-length sequences around transcription factor binding sites (TFBSs) identified by ChIP-Seq. If the file name ends with `.gz`, the file is treated as gzipped FASTA. If `-i stdin` is specified, MOCCS receive FASTA from the standard input.
	- We recommend to use genomic sequences ±350 bp around the TFBSs.
- `-k <int>`: Length of k-mers to calculate AUCs. Mutually exclusive with `--regex`.
	- We recommend to set k to 5, 6, 7, or 8.
- `--regex <string>`: (Optional) Regular expression (e.g. CANNTG). When this option is specified, MOCCS calculates AUCs for only k-mers corresponding to the regular expression. Currently, 'N' is interpreted as {A,C,G,T}. Mutually exclusive with `-k`.
- `--label <string>`: (Optional) Label used as the prefix of output files. Default is 'MOCCS_result'.
- `--mask`: (Optional) Mask lower-case characters in fasta file as 'N'.
- `--threshold <float>`: (Optional) Only print k-mers with AUCs higher than this value.
- `--stranded`: (Optional) Count k-mers on the forward strand. This is useful when you analyze sequences bound by RNA-binding proteins (RBPs) (e.g. CLIP-Seq data). By default, MOCCS sums counts of a k-mer and its reverse-complement before calculation of the AUC score.
- `--low-count-threshold <float>`: (Optional) Set the low-count threshold. Those k-mers whose counts are less than this value will be filtered. By default, MOCCS set the low-count threshold to `N * (w - k + 1) / 4^k`.

Note that `MOCCS_visualize.r` must be located on the same directory as `MOCCS.pl`.


## Example
### Calculate AUCs for all 6-mers

	perl MOCCS.pl -i test_data/test_701bp.fa -k 6 --label test_out_0/test_out_0

### Calculate AUCs for all 6-mers while thresholding k-mers with low AUCs

	perl MOCCS.pl -i test_data/test_701bp.fa -k 6 --label test_out_1/test_out_1 --threshold 10

### Calculate AUCs for all 6-mers while masking repeats (lower-case)

	perl MOCCS.pl -i test_data/test_701bp.fa -k 6 --label test_out_2/test_out_2 --mask

### Calculate AUCs for CANNTG

	perl MOCCS.pl -i test_data/test_701bp.fa --regex CANNTG --label test_out_3/test_out_3

### Calculate AUCs for CANNTG while masking fasta

	perl MOCCS.pl -i test_data/test_701bp.fa --regex CANNTG --label test_out_4/test_out_4 --mask

### Calculate AUCs for all 6-mers on the forward strand (Useful for RBP binding data)

	perl MOCCS.pl -i test_data/test_701bp.fa -k 6 --label test_out_5/test_out_5 --stranded

### Calculate AUCs for all 6-mers with a fixed low-count threshold

	perl MOCCS.pl -i test_data/test_701bp.fa -k 6 --label test_out_6/test_out_6 --low-count-threshold 100

## MOCCS2 score
If a k-mer sequence has a low appearance count, the AUC calculated by MOCCS for the k-mer is unstable. To overcome this issue, we (Hikari Yoshitane, Yoshimasa Asano, Wataru Iwasaki, and Haruka Ozaki) devised the MOCCS2 score, which is a relative value of AUC normalized by the standard deviation (SD) at its appearance count. From MOCCS version 2.0 (abbreviated as MOCCS2), MOCCS calculates MOCCS2 score for each k-mer sequence as well as an AUC.

Let $C$ be the appearance count of the k-mer sequence and Let $W$ be the size of the analyzed window where k-mer sequences are sought at around ChIP-peak positions, MOCCS2 score is defined as follows:

$$\textrm{[MOCCS2 score]} = \textrm{[AUC]} / \textrm{[SD of AUC]} = \textrm{[AUC]} \times \frac{\sqrt{12C}}{W}$$

The calculation of [SD of AUC] was mathematically derived as follows:

> If a k-mer sequence appears only once at a random position within the window, its coordinate follows the uniform distribution $U(0, W)$, whose variance is known to be $W^2 / 12$. Because i) AUC is calculated by subtracting $W / 2$ from the coordinate and ii) constant subtraction does not affect variance of probability distributions, variance of AUC is also $W^2 / 12$ if the appearance count is $1$.
>
> Next, assume that a k-mer sequence appears $C$ times at random positions within the window. The variance of the sum of their coordinates becomes $CW^2 / 12$, because variance of sum of random variables that follow the same probability distribution is proportional to the numbers of the variables.  Then, because AUC is calculated by dividing the sum of their coordinates by C and subtracting $W/2$, the variance of AUC is
> $$(CW^2 /12)/C^2 =W^2 /12C ,$$
> if the appearance count is $C$. Finally, we obtain [SD of AUC] by taking the square root of the variance:
>
> $$[𝑆𝐷 𝑜𝑓 𝐴𝑈𝐶] = \frac{𝑊}{\sqrt{12C}}$$

In the current implementation of MOCCS2, $W$ was set to `floor((sequence_length - 1)/2) + 1 – floor(k / 2)`. Note that since the value of $W$ is the same for every $k$-mer, the relative ranks among k-mer sequences are not affected by the choice of the definition of $W$.

## Citation

- Ozaki H, Iwasaki W. MOCCS: Clarifying DNA-binding motif ambiguity using ChIP-seq data. Computational biology and chemistry. 2016 Aug 1;63:62-72. [PubMed](https://www.ncbi.nlm.nih.gov/pubmed/26971251) [Preprint](http://yuifu.github.io/pdf/2016_moccs.pdf)
- A paper describing MOCCS2 score will be published soon.

## Contact

[Haruka Ozaki](https://yuifu.github.io/) <harukao.cb@gmail.com>

## License

Copyright (c) [2016] [Haruka Ozaki]
This software is released under the MIT License, see LICENSE.txt.
