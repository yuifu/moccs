# MOCCS

Motif Centrality Analysis of ChIP-Seq (MOCCS) is a method for for clarifying DNA-binding motif ambiguity.

Given ChIP-Seq data of any DNA-binding proteins including transcription factors, MOCCS comprehensively analyzes and describes every $k$-mer that is bound by the DNA-binding proteins. 

MOCCS (version 1.4) is written in Perl and R.  
MOCCS was tested on Perl version 5.18.2 and R version 3.2.1.

# Usage

	perl MOCCS.pl -i <string> -k <int>|--regex <string> [--label <string>] [--mask] [--threshold <float>]

- `-i <string>`: FASTA file containing fixed-length sequences around transcription factor binding sites (TFBSs) identified by ChIP-Seq. If the file name ends with `.gz`, the file is treated as gzipped FASTA. If `-i stdin` is specified, MOCCS receive FASTA from the standard input.
	- We recommend to use genomic sequences ±350 bp around the TFBSs.
- `-k <int>`: Length of k-mers to calculate AUCs. Mutually exclusive with `--regex`.
	- We recommend to set k to 5, 6, 7, or 8.
- `--regex <string>`: (Optional) Regular expression (e.g. CANNTG). When this option is specified, MOCCS calculates AUCs for only k-mers corresponding to the regular expression. Currently, 'N' is interpreted as {A,C,G,T}. Mutually exclusive with `-k`.
- `--label <string>`: (Optional) Label used as the prefix of output files. Default is 'MOCCS_result'.
- `--mask`: (Optional) Mask lower-case characters in fasta file as 'N'.
- `--threshold <float>`: (Optional) Only print k-mers with AUCs higher than this value.

Note that `MOCCS_visualize.r` must be located on the same directory as `MOCCS.pl`.


# Example
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


# Contact

Haruka Ozaki <harukao.cb@gmail.com>

# License 

Copyright (c) [2015] [Haruka Ozaki]
This software is released under the MIT License, see LICENSE.txt.
