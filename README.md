# MOCCS

Motif Centrality Analysis of ChIP-Seq (MOCCS) is a method for for clarifying DNA-binding motif ambiguity.

Given ChIP-Seq data of any DNA-binding proteins including transcription factors, MOCCS comprehensively analyzes and describes every $k$-mer that is bound by the DNA-binding proteins. 

# Usage

	perl MOCCS.pl <fasta> <k or regex> <out_prefix> [-mask]

- `<fasta>`: Fasta file containing fixed-length sequences around transcription factor binding sites (TFBSs) identified by ChIP-Seq.
- `<k>`: Length of k-mers to calculate AUCs.
- `<regex>`: Regular expression (e.g. CANNTG). Currently, 'N' is interpreted as {A,C,G,T}.
- `<out_prefix>`: Label used as the prefix of output files.
- `-mask`: (Optional) Mask lower-case characters in fasta file as 'N'.

Note that `MOCCS_visualize.r` must be located on the same directory as `MOCCS.pl`.


# Example
## Calculate AUCs for all 6-mers

	perl MOCCS.pl test_data/test_701bp.fa 6 test_out_1/test_out_1

## Calculate AUCs for all 6-mers while masking repeats (lower-case)

	perl MOCCS.pl test_data/test_701bp.fa 6 test_out_2/test_out_2 -mask

## Calculate AUCs for CANNTG

	perl MOCCS.pl test_data/test_701bp.fa CANNTG test_out_3/test_out_3

## Calculate AUCs for CANNTG while masking fasta
	
	perl MOCCS.pl test_data/test_701bp.fa CANNTG test_out_4/test_out_4 -mask

# Dependency

MOCCS is written in Perl and R.
MOCCS was tested on Perl version 5.18.2 and R version 3.2.1.

# Version 


# Contact

Haruka Ozaki <harukao.cb@gmail.com>

# License 

Copyright (c) [2015] [Haruka Ozaki]
This software is released under the MIT License, see LICENSE.txt.
