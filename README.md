# Simple Codon Optimizer README #

### Description ###
Program for converting nucleotide or amino acid sequences into an "optimized" nucleotide sequence using nucleotide triplet information.

Automatically downloads codon translation tables from NCBI.

Can parse several formats of codon usage tables (see included examples).

### Requirements ###

 * Python >= 3.6 ([source](https://www.python.org/downloads/), [binaries](https://www.python.org/downloads/), [documentation](https://docs.python.org/3/))
 * requests module ([source](https://github.com/kennethreitz/requests), [whls](https://pypi.org/project/requests/), [documentation](http://docs.python-requests.org/))

### Program usage ###
This will print out command line parameter descriptions and examples.
```sh
$ python3 simple-codom-optimizer.py --help
```

It returns the following output.
> ```
> usage: simple-codon-optimizer.py [-h] [--deterministic] [--samples N]
>                                  [--display N] [--suppress]
>                                  USAGE/EXPRESSION_TABLE TRANSLATION_TABLE
>                                  SEQUENCE
> 
> description:
>   Simple program for optimizing a protein-coding sequence.
> 
>   Several formats for codon usage table are supported (See included example files).
>   Additionally, a gene expression table can be provided to base codon optimality from.
> 
>   Output is in the format (FREQUENCY, SEQUENCE).
> 
>   The program will first check to see if the input SEQUENCE is composed
>   exclusively of 'nt' characters. If it is not, then it will check to
>   see if it is made of 'aa' characters. Space (' ') characters are
>   allowed in SEQUENCE.
> 
> positional arguments:
>   USAGE/EXPRESSION_TABLE
>                         File containing either the codon usage table (counts),
>                         or a gene expression table.
>   TRANSLATION_TABLE     The translation table id.
>   SEQUENCE              'nt' or 'aa' sequence to optimize.
> 
> optional arguments:
>   -h, --help            Show this help message and exit.
>   --deterministic       Instead of calculating a distribution of sequences,
>                         just find the single most-optimal sequence. (default:
>                         False)
>   --samples N           Number of sequences to generate. (default: 100000)
>   --display N           Number of output sequences to display. (default: 10)
>   --suppress            Suppress STDERR messages. (default: False)
> 
> examples:
>   python3 simple-codon-optimizer.py examples/5501_codons.txt 1 ASRWLAQC
>   python3 simple-codon-optimizer.py "examples/Codon usage table 5501.html" 1 "GCA TCA AGA TGG CTG GCG CAA TGT"
>   python3 simple-codon-optimizer.py examples/C_albicans_codon_usage.tab 12 EGRGSLLTCGDVEENPGP --deterministic
> ```

A basic program usage using DNA triplets as input would look like this.
```sh
$ python3 simple-codon-optimizer.py "examples/Codon usage table 5501.html" 1 "GCA TCA AGA TGG CTG GCG CAA TGT" --suppress
```

Which returns the following.
> ```
> 0.0009 GCCTCCCGCTGGCTCGCCCAGTGC
> 0.00073 GCCAGCCGCTGGCTCGCCCAGTGC
> 0.00072 GCCTCCCGCTGGCTCGCCCAGTGT
> 0.00066 GCCTCCCGCTGGCTCGCTCAGTGC
> 0.00065 GCCTCTCGCTGGCTCGCCCAGTGC
> 0.00065 GCCTCTCGTTGGCTCGCCCAGTGC
> 0.00062 GCCTCCCGTTGGCTCGCCCAGTGC
> 0.00061 GCCAGCCGCTGGCTCGCTCAGTGC
> 0.0006 GCCTCTCGCTGGCTCGCTCAGTGC
> 0.00059 GCTAGCCGCTGGCTCGCCCAGTGC
> ```

Here's a simple example using an amino acid sequence as input without the `--suppress` command line option.
```sh
$ python3 simple-codon-optimizer.py examples/5501_codons.txt 1 ASRWLAQC
```

This creates the following (verbose) output.
> ```
> Loading translation tabels.
> Parsing codon usage table.
> Usage table:
> TTT 291368 0.016864017795171243
> TTC 349064 0.02020339058391332
> TTA 140536 0.008134049054330558
> TTG 281524 0.016294259307019953
> CTT 337452 0.01953130245262392
> CTC 343944 0.019907051345866324
> CTA 163988 0.00949142167360363
> CTG 278148 0.016098860621932717
> ATT 329656 0.019080079659691426
> ATC 377884 0.021871456373076283
> ATA 168296 0.009740763360616607
> ATG 368400 0.02132253423760017
> GTT 313632 0.018152630450616224
> GTC 314600 0.018208657087809485
> GTA 121956 0.007058661741261581
> GTG 258528 0.014963279401135442
> TAT 237604 0.013752224280648075
> TAC 222076 0.012853482935258674
> TAA 13436 0.0007776589848436371
> TAG 9092 0.00052623366256314
> CAT 232416 0.013451949287095769
> CAC 190644 0.011034237831685795
> CAA 326356 0.01888907976016895
> CAG 370060 0.02141861297493572
> AAT 331276 0.019173843246729733
> AAC 319832 0.018511478746688757
> AAA 402628 0.023303608346950274
> AAG 481472 0.02786700109784427
> GAT 531948 0.030788489257934135
> GAC 419976 0.024307689030864194
> GAA 560712 0.032453313836680965
> GAG 524856 0.030378013110233103
> TCT 287384 0.016633428825565927
> TCC 298932 0.017301812716379733
> TCA 230612 0.013347536008690148
> TCG 229684 0.01329382452179413
> CCT 283644 0.016416962272773786
> CCC 248512 0.014383565766706009
> CCA 306240 0.017724790675685876
> CCG 240160 0.013900162384641849
> ACT 235828 0.013649431607450524
> ACC 289932 0.016780903899500252
> ACA 244556 0.01415459740230876
> ACG 206384 0.01194524947365058
> GCT 381676 0.02209093262125484
> GCC 381452 0.02207796777959028
> GCA 338976 0.01961950967894885
> GCG 282004 0.016322041110586858
> TGT 94896 0.005492462565177269
> TGC 131700 0.007622632353669768
> TGA 17112 0.0009904212971601903
> TGG 234756 0.013587385579484435
> CGT 168068 0.009727567003922327
> CGC 215888 0.012495329184275313
> CGA 210676 0.012193665100544662
> CGG 167304 0.009683347633245003
> AGT 172268 0.009970657785132753
> AGC 264248 0.01529434589364107
> AGA 208276 0.012054756082710134
> AGG 157392 0.009109653389588399
> GGT 272268 0.015758533528238118
> GGC 347092 0.020089253674259278
> GGA 315272 0.018247551612803153
> GGG 200944 0.011630389033225648
> Treating input sequence as 'aa'.
> Confusion matrix with potential codons for each aa in sequence.
>   aa   A   S   R   W   L   A   Q   C
> high  GCT TCC CGC TGG CTC GCT CAG TGC
>       GCC TCT CGA     CTT GCC CAA TGT
>       GCA AGC AGA     TTG GCA
>       GCG TCA CGT     CTG GCG
>           TCG CGG     CTA
> low       AGT AGG     TTA
> 0.00033 GCTTCTAGATGGCTTGCCCAGTGC
> 0.00031 GCATCTCGCTGGCTTGCTCAGTGC
> 0.0003 GCCTCTAGATGGCTTGCTCAATGC
> 0.00029 GCTTCCCGATGGCTCGCCCAGTGC
> 0.00028 GCCTCTCGATGGCTTGCCCAGTGC
> 0.00028 GCTTCTCGCTGGTTGGCCCAGTGC
> 0.00028 GCTTCTAGATGGCTCGCTCAGTGC
> 0.00027 GCATCCCGATGGCTCGCCCAGTGC
> 0.00026 GCTTCTCGATGGCTTGCTCAGTGC
> 0.00026 GCATCTCGCTGGCTCGCTCAATGC
> Simple Codon Optimizer finished.
> ```

An experimental feature allows you to codon-optimize using gene expression data instead of just the codon usage frequencies.
The example file `examples/expressions.txt` contains some simple test data:

> ```text
> # gene     expression  sequence
> gene_0001  10          GCA TAC GCA TAT AGG GCT CCT CAG ACT CAG CCA TAT GGC TCT TAA
> gene_0002  100         GCA AAA TAT AAT GCC TAG AAT GGC ATA GGC CAT GGT AAG AAT GAC TCT CGC TAT TAA
> gene_0003  1000        GTG ATG AAT GGT TTC GTT GAA TTA GCT ATA GGT AAG ATG AGA TAC GCT AAG CGC TAC GAT CAT AGC CAT TAG
> gene_0004  125         GCT CCG CTA GCT CAG ATA TAG AAT CGC CTC GCT CGC TGA
> ```

We specify an amino acid input with multiple alanine residues:

```sh
$ python3 simple-codon-optimizer.py examples/expressions.txt 1 ASRWLAQC
```

And we see the alanine codon associated with the highest expression changes depending on its location in the sequence.
> ```
> Loading translation tables.
> Treating input sequence as 'aa'.
> Parsing input file.
> No valid codon usage table discovered. Assuming input is an expression table.
> Confusion matrix with potential codons for each aa in sequence.
>   aa   A   S   R   W   L   A   Q   C
> high  GCA AGC AGA TGG TTA GCT CAA TGT
>       GCG TCG CGA     CTA GCG CAG TGC
>       GCT TCC CGG     CTT GCC
>       GCC TCA CGT     CTG GCA
>           AGT AGG     TTG
> low       TCT CGC     CTC
> Simple Codon Optimizer finished.
> ```

### Who do I talk to? ###

Thaddeus D. Seher ([@tdseher](https://twitter.com/tdseher))

