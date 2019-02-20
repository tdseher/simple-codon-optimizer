#!/usr/bin/env python3

'''
Simple Codon Optimizer
Copyright 2019 Thaddeus D. Seher (@tdseher)
'''

import sys
import random
from collections import defaultdict
import requests
import re
import itertools
import operator
import argparse

def parse_translation_tables(text):
    '''
    The format that is parsed:
        1. The Standard Code (transl_table=1)
          AAs  = FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
        Starts = ---M------**--*----M---------------M----------------------------
        Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
        Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
        Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
        ...
    Returns dict.
    '''
    
    matches = re.findall(
        r'\d+\.\s+(.*?)\s*\(transl_table=(\d+)\)' # Table name, id number
        r'(?:.|\n|\r\n?)*?' # Line break and anything else
        r'\s*AAs\s*\=\s*(\S+)' # Amino acids
        r'(?:\n|\r\n?)' # Line break
        r'\s*Starts\s*=\s*(\S+)' # Initiation and termination codons
        r'(?:\n|\r\n?)' # Line break
        r'\s*Base1\s*=\s*(\S+)' # First position
        r'(?:\n|\r\n?)' # Line break
        r'\s*Base2\s*=\s*(\S+)' # Second position
        r'(?:\n|\r\n?)' # Line break
        r'\s*Base3\s*=\s*(\S+)', # Third position
        text,
        re.MULTILINE
    )
    
    tables = {}
    
    for m in matches:
        name, table_id, aas, starts, base1, base2, base3 = m
        aa_to_codon_d = defaultdict(list)
        codon_to_aa = {}
        for i in range(len(aas)):
            aa = aas[i]
            codon = base1[i] + base2[i] + base3[i]
            aa_to_codon_d[aa].append(codon)
            codon_to_aa[codon] = aa
        
        # Convert to a regular dictionary with a tuple (instead of a list) as value
        aa_to_codon = {}
        for k, v in aa_to_codon_d.items():
            aa_to_codon[k] = tuple(v)
        
        tables[int(table_id)] = (aa_to_codon, codon_to_aa)
    
    return tables

def download(url):
    '''
    Use requests library to get the text from the input URL.
    '''
    r = requests.get(url, allow_redirects=True)
    return r.text

def nt_to_aa(sequence, table):
    '''
    Translate 'nt' string into 'aa' string given input codon to 'aa' table.
    '''
    aa_list = []
    for i in range(len(sequence)//3):
        codon = sequence[3*i:3*(i+1)]
        aa_list.append(table[codon])
    return ''.join(aa_list)

def is_url(text):
    '''
    Returns 'True' if the text is a URL.
    '''
    # Function unfinished
    return False

def is_file(text):
    '''
    Returns 'True' if the text is a file.
    '''
    # Function unfinished
    return True

def parse_format_a(text):
    '''
    Parses text resembling the following:
        fields: [triplet] [frequency: per thousand] ([number])
        UUU 15.3(   201)  UCU 15.7(   207)  UAU 15.0(   198)  UGU  6.2(    82)
        UUC 24.8(   327)  UCC 15.6(   205)  UAC 16.5(   217)  UGC  7.7(   101)
        UUA  4.7(    62)  UCA  7.8(   102)  UAA  0.9(    12)  UGA  1.0(    13)
        ...
    Returns dict.
    '''
    table = {}
    for m in re.findall(r'([ACGTU]{3})\s+(\d+(?:\.\d*))?\(\s*(\d+)\)', text):
        codon, freq_per_thousand, number = m[0], float(m[1]), int(m[2])
        codon = codon.replace('U', 'T')
        codon = codon.replace('u', 't')
        table[codon] = number
    return table

def parse_format_b(text):
    '''
    Parses this type of format:
        id,value
        taxid,5501
        collapse,4
        "#codon",17277496
        "#CDS",39640
        "GC%",50.78
        "GC1%",55.85
        "GC2%",44.36
        "GC3%",52.13
        TTT,291368
        TTC,349064
        TTA,140536
        TTG,281524
        CTT,337452
        ...
    Returns dict.
    '''
    table = {}
    for m in re.findall(r'([ACGTU]{3}),(\d+)', text):
        codon, number = m[0], int(m[1])
        table[codon] = number
    return table

def parse_format_c_old(text, taxid):
    '''
    Expects a header row, and one or more data rows
        Division	Assembly	Taxid	Species	Organelle	Translation Table	# CDS	# Codons	GC%	GC1%	GC2%	GC3%	TTT	TTC	TTA	TTG	CTT	CTC	CTA	CTG	ATT	ATC	ATA	ATG	GTT	GTC	GTA	GTG	TAT	TAC	TAA	TAG	CAT	CAC	CAA	CAG	AAT	AAC	AAA	AAG	GAT	GAC	GAA	GAG	TCT	TCC	TCA	TCG	CCT	CCC	CCA	CCG	ACT	ACC	ACA	ACG	GCT	GCC	GCA	GCG	TGT	TGC	TGA	TGG	CGT	CGC	CGA	CGG	AGT	AGC	AGA	AGG	GGT	GGC	GGA	GGG
        refseq	GCF_001560135.1	2285	Sulfolobus acidocaldarius	genomic	11	2222	630256	37.4	44.66	34.12	33.43	16426	10487	22724	8884	10561	4921	12006	5609	18559	6976	31791	13196	16278	5953	17766	9203	18246	10880	1057	404	5332	2753	8051	6151	18800	11435	25417	21808	19934	10723	22786	20113	9875	4893	11260	2699	9911	3599	8936	1975	11237	4662	11584	3020	13581	4705	13301	3109	2632	1070	759	6161	875	250	749	232	10720	4413	16339	10368	15101	4397	16981	5632
        refseq	GCF_000337915.1	1227484	Halorubrum saccharovorum DSM 1137	genomic	11	3122	944967	68.22	70.52	45.73	88.4	1744	28488	1174	3486	3128	47657	1322	23536	2635	32777	1396	14801	4504	55008	1983	23947	2050	22062	479	752	1364	16193	2625	17216	2022	18969	3398	12634	10035	74485	14521	66439	1839	11412	1816	20789	1231	15801	1849	26336	2542	25702	2433	27142	4205	43379	4495	57612	2423	3758	1886	10155	3425	27542	7166	23443	2417	12720	1237	1426	7294	44759	9016	24887
        ...
    Returns dict.
    '''
    table = {}
    header = None
    for line in text.splitlines():
        sline = line.split('\t')
        if (header == None):
            header = sline
        else:
            if (taxid == int(sline[2])):
                for p in itertools.product('ACGT', repeat=3):
                    codon = ''.join(p)
                    i = header.index(codon)
                    number = int(sline[i])
                    table[codon] = number
    return table

def parse_format_c(text):
    '''
    Expects a header row, and a single data row
        Division	Assembly	Taxid	Species	Organelle	Translation Table	# CDS	# Codons	GC%	GC1%	GC2%	GC3%	TTT	TTC	TTA	TTG	CTT	CTC	CTA	CTG	ATT	ATC	ATA	ATG	GTT	GTC	GTA	GTG	TAT	TAC	TAA	TAG	CAT	CAC	CAA	CAG	AAT	AAC	AAA	AAG	GAT	GAC	GAA	GAG	TCT	TCC	TCA	TCG	CCT	CCC	CCA	CCG	ACT	ACC	ACA	ACG	GCT	GCC	GCA	GCG	TGT	TGC	TGA	TGG	CGT	CGC	CGA	CGG	AGT	AGC	AGA	AGG	GGT	GGC	GGA	GGG
        refseq	GCF_001560135.1	2285	Sulfolobus acidocaldarius	genomic	11	2222	630256	37.4	44.66	34.12	33.43	16426	10487	22724	8884	10561	4921	12006	5609	18559	6976	31791	13196	16278	5953	17766	9203	18246	10880	1057	404	5332	2753	8051	6151	18800	11435	25417	21808	19934	10723	22786	20113	9875	4893	11260	2699	9911	3599	8936	1975	11237	4662	11584	3020	13581	4705	13301	3109	2632	1070	759	6161	875	250	749	232	10720	4413	16339	10368	15101	4397	16981	5632
    Returns dict.
    '''
    table = {}
    header = None
    for line in text.splitlines():
        sline = line.split('\t')
        if (header == None):
            header = sline
        else:
            for p in itertools.product('ACGT', repeat=3):
                codon = ''.join(p)
                i = header.index(codon)
                number = int(sline[i])
                table[codon] = number
    return table

def is_nt(sequence):
    '''
    Returns 'True' if the sequence uses exclusively nucleotide characters.
    '''
    # IUPAC Codes for Nucleotides
    # Symbol  Description           Bases Represented
    # A       adenosine/adenine     A---
    # C       cytidine/cytosine     -C--
    # G       guanosine/guanine     --G-
    # T       thymidine/thymine     ---T
    # U       uridine/uracil        ---U
    # W       weak                  A--T
    # S       strong                -CG-
    # M       amino                 AC--
    # K       keto                  --GT
    # R       purine                A-G-
    # Y       pyrimidine            -C-T
    # B       not A                 -CGT
    # D       not C                 A-GT
    # H       not G                 AC-T
    # V       not T                 ACGT
    # N or -  any base (not a gap)  ACGT
    # The most common, non-standard nucleotide codes are "I" (Inosine) then "X" (xanthine)
    for m in re.finditer(r'[^acgturymkwsbdhvnACGTURYMKWSBDHVN.-]', sequence):
        if m:
            return False
    return True

def is_aa(sequence):
    '''
    Returns 'True' if the sequence uses exclusively amino acid characters.
    '''
    # IUPAC Codes for Amino Acids
    # 1-Letter Code  3-Letter Code  Amino Acid
    # A              Ala            Alanine
    # C              Cys            Cysteine
    # D              Asp            Aspartic Acid
    # E              Glu            Glutamic Acid
    # F              Phe            Phenylalanine
    # G              Gly            Glycine
    # H              His            Histidine
    # I              Ile            Isoleucine
    # K              Lys            Lysine
    # L              Leu            Leucine
    # M              Met            Methionine
    # N              Asn            Asparagine
    # P              Pro            Proline
    # Q              Gln            Glutamine
    # R              Arg            Arginine
    # S              Ser            Serine
    # T              Thr            Threonine
    # V              Val            Valine
    # W              Trp            Tryptophan
    # X              Xaa            Unspecified or unknown
    # Y              Tyr            Tyrosine
    # *                             STOP
    
    # not included
    # B              Asx            Aspartic Acid or Asparagine
    # J              Xle            Leucine or Isoleucine
    # O              Pyl            Pyrrolysine
    # U              Sec            Selenocysteine
    # Z              Glx            Glutamic Acid or Glutamine
    
    for m in re.finditer(r'[^acdefghiklmnpqrstvwxyACDEFGHIKLMNPQRSTVWXY*.-]', sequence):
        if m:
            return False
    return True

def stochastic_aa(aa_sequence, table):
    '''
    Converts input amino acid sequence into a stochastic nucleotide string.
    '''
    #table = {
    #    'L': (('TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'), (14.3, 13.0, 11.9, 10.2, 4.2, 48.4))
    #}
    
    nt_list = []
    for aa in aa_sequence:
        codons, frequencies = table[aa]
        nt_list.append(random.choices(codons, weights=frequencies, k=1)[0])
        
    return ''.join(nt_list)

def deterministic_aa(aa_sequence, table):
    '''
    Returns highest frequency nucleotide string.
    '''
    nt_list = []
    for aa in aa_sequence:
        codons, frequencies = table[aa]
        nt_list.append(sorted(zip(frequencies, codons), reverse=True)[0][1])
    
    return ''.join(nt_list)

def potential_codons(aa_sequence, table):
    '''
    Returns multi-line string formatted as follows:
          aa   A   S   R   W   L   A   Q   C
        high  GCT TCC CGC TGG CTC GCT CAG TGC
              GCC TCT CGA     CTT GCC CAA TGT
              GCA AGC AGA     TTG GCA
              GCG TCA CGT     CTG GCG
                  TCG CGG     CTA
        low       AGT AGG     TTA
    '''
    columns = []
    
    for i, aa in enumerate(aa_sequence):
        codons, freqs = table[aa]
        cf = sorted(zip(freqs, codons), reverse=True)
        columns.append([x[1] for x in cf])
    
    outputs = []
    outputs.append('  aa ')
    for i, aa in enumerate(aa_sequence):
        outputs[-1] += '  ' + aa + ' '
    
    outputs.append('high ')
    for i, aa in enumerate(aa_sequence):
        if (len(columns[i]) > 0):
            outputs[-1] += ' ' + columns[i].pop(0)
        else:
            outputs[-1] += '    '
    
    while any([len(x)>1 for x in columns]):
        outputs.append('     ')
        for i, aa in enumerate(aa_sequence):
            if (len(columns[i]) > 0):
                outputs[-1] += ' ' + columns[i].pop(0)
            else:
                outputs[-1] += '    '
    
    outputs.append('low  ')
    for i, aa in enumerate(aa_sequence):
        if (len(columns[i]) > 0):
            outputs[-1] += ' ' + columns[i].pop(0)
        else:
            outputs[-1] += '    '
    
    return '\n'.join(outputs)

def process(args, usage_text, translation_table_number, sequence, ignore_mask=True):
    # Download the translation tables
    print("Loading translation tabels.", file=sys.stderr) if not args.suppress else None
    text = download('https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi')
    translation_tables = parse_translation_tables(text)
    aa_to_codon_table, codon_to_aa_table = translation_tables[translation_table_number]
    
    # Add 'X' 'any' character to table
    aa_to_codon_table['X'] = tuple(''.join(x) for x in itertools.product('ACGT', repeat=3))
    
    usage_table = {}
    
    # Try parsing the table each way to store codon usage
    print("Parsing codon usage table.", file=sys.stderr) if not args.suppress else None
    for parse_type in [parse_format_a, parse_format_b, parse_format_c]:
        if (len(usage_table) != 64):
            try:
                usage_table = parse_type(usage_text)
            except:
                pass
        else:
            break
    
    # Count the sum of all codons
    n = sum(usage_table.values())
    
    # Let user know if the table was parsed successfully
    print('Usage table:', file=sys.stderr) if not args.suppress else None
    for k, v in usage_table.items():
        print(k, v, v/n, file=sys.stderr) if not args.suppress else None
    
    # Remove any spaces
    sequence = sequence.replace(' ', '')
    
    # Deal with case-sensitivity
    if ignore_mask:
        sequence = sequence.upper()
    
    if is_nt(sequence):
        print("Treating input sequence as 'nt'.", file=sys.stderr) if not args.suppress else None
        # If it is a DNA sequence, then convert Ts to Us
        sequence = sequence.replace('U', 'T')
        sequence = sequence.replace('u', 't')
        
        # Convert it to an aa sequence
        aa_seq = nt_to_aa(sequence, codon_to_aa_table)
    
    elif is_aa(sequence):
        print("Treating input sequence as 'aa'.", file=sys.stderr) if not args.suppress else None
        aa_seq = sequence
    else:
        sys.exit("Sequence contains invalid characters.")
    
    # Link codon frequencies
    # aa_to_codon_freq_table = {
    #     'L': (('TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'), (14.3, 13.0, 11.9, 10.2, 4.2, 48.4)),
    # }
    aa_to_codon_freq_table = {}
    for aa, codons in aa_to_codon_table.items():
        freqs = []
        for codon in codons:
            freqs.append(usage_table[codon]/n)
        
        aa_to_codon_freq_table[aa] = (codons, tuple(freqs))
    
    # Link aa frequencies
    # codon_to_aa_freq_table = {
    #     'TTA': ('L', 14.3),
    #     'TTG': ('L', 13.0),
    #     'CTT': ('L', 11.9),
    #     'CTC': ('L', 10.2),
    #     'CTA': ('L', 4.2),
    #     'CTG': ('L', 48.4),
    # }
    codon_to_aa_freq_table = {}
    for codon, aa in codon_to_aa_table.items():
        codon_to_aa_freq_table[codon] = (aa, usage_table[codon]/n)
    
    # Let user know the aa sequence that will be translated
    print('Confusion matrix with potential codons for each aa in sequence.', file=sys.stderr) if not args.suppress else None
    print(potential_codons(aa_seq, aa_to_codon_freq_table), file=sys.stderr) if not args.suppress else None
    
    # Once there is an aa sequence, then for each character,
    # pick an optimal codon
    output = defaultdict(int)
    
    samples = args.samples
    
    if args.deterministic:
        samples = 1
        output[deterministic_aa(aa_seq, aa_to_codon_freq_table)] += 1
    else:
        for i in range(samples):
            output[stochastic_aa(aa_seq, aa_to_codon_freq_table)] += 1
        
    # Write results to STDOUT
    #for k, v in sorted(output, key=lambda x: output[x], reverse=True):
    for k, v in sorted(output.items(), key=operator.itemgetter(1), reverse=True)[:max(0, args.display)]:
        print(v/samples, k)

def main():
    # Create parser
    parser = argparse.ArgumentParser(
        description=(
            "Optimize 'aa' or 'nt' sequence." "\n"
            "3 formats for codon usage table are supported." "\n"
            "Output is in the format (FREQUENCY, SEQUENCE)." "\n"
        ),
        epilog=(
            "example:" "\n"
            "  (python3 simple-codon-optimizer.py 5501_codons.txt 1 ASRwLAQC)" "\n"
            "example:" "\n"
            '  (python3 simple-codon-optimizer.py "Codon usage table 5501.html" 1 "GCA TCA AGA TGG CTG GCG CAA TGT")' "\n"
        )
    )
    parser.add_argument('usage_table', metavar='USAGE_TABLE', type=str,
        help='File containing the codon usage table (counts).')
    parser.add_argument('translation_table', metavar='TRANSLATION_TABLE', type=int,
        help='The translation table id.')
    parser.add_argument('sequence', metavar='SEQUENCE', type=str,
        help="'nt' or 'aa' sequence to optimize.")
    parser.add_argument('--deterministic', action='store_true',
        help='Instead of calculating a distribution of sequences, \
             just find the single most-optimal sequence.')
    parser.add_argument('--samples', metavar='N', type=int, default=100000,
        help='Number of sequences to generate.')
    parser.add_argument('--display', metavar='N', type=int, default=10,
        help='Number of output sequences to display.')
    parser.add_argument('--suppress', action='store_true',
        help='Suppress STDERR messages.')
    
    args = parser.parse_args()
    
    # User specifies the table via text file or url
    # Format 'a' downloaded from r'https://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=5501'
    #usage_table = r'Codon usage table 5501.html'
    
    # Format 'b' downloaded from r'https://hive.biochemistry.gwu.edu/dna.cgi?cmd=ionTaxidCollapse&svcType=svc-refseq-processor&fileSource=refseq_species.tsv&taxid=5501&filterInColName=[%22Organelle%22]&filterIn=[%22genomic%22]&searchDeep=true&raw=1&raw=1'
    #usage_table = '5501_codons.txt'
    
    # Format 'c' downloaded from r'https://hive.biochemistry.gwu.edu/dna.cgi?cmd=objFile&ids=569942&filename=refseq_species.tsv&raw=1'
    #usage_table = r'o569942-refseq-GCF_000149335.2.txt'
    
    # User specifies the translation table
    #translation_table = 1
    
    # User specifies either the aa sequence, RNA exon sequence, or genomic DNA sequence
    #sequence = ' A   S   R   w   L   A   Q   C'
    #sequence = 'GCA TCA AGA TGG CTG GCG CAA TGT'
    
    if is_url(args.usage_table):
        # Download this url
        usage_text = download(url)
    elif is_file(args.usage_table):
        # Open this file
        with open(args.usage_table,'r') as flo:
            usage_text = flo.read()
    
    # Translate
    process(args, usage_text, args.translation_table, args.sequence)
    
if (__name__ == '__main__'):
    main()
