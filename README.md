# pairwise-alignment

This is a Python module to calculate a pairwise alignment between biological sequences (protein or nucleic acid). This module uses the [needle](https://www.ebi.ac.uk/Tools/psa/emboss_needle/) and [water](https://www.ebi.ac.uk/Tools/psa/emboss_water/) tools from the EMBOSS package to calculate an optimal, global/local pairwise alignment.

I wrote this module for two reasons. First, the needle and water tools are faster than any Python implementation. Second, Biopython has dropped support for tools from the EMBOSS package and recommends running them via the subprocess module directly.

## Table of contents

* [Introduction](#introduction)
* [Requirements](#requirements)
* [Quick Start](#quick-start)
* [Alignment object](#alignment-object)
   * [Attributes](#attributes)
   * [Methods](#methods)
* [Usage examples](#usage-examples)
   * [Alignment information](#alignment-information)
   * [Alignment in raw format](#alignment-in-raw-format)
   * [Alignment in FASTA format](#alignment-in-fasta-format)
   * [Alignment iteration](#alignment-iteration)
   * [Alignment iteration with index](#alignment-iteration-with-index)
   * [Query coverage](#query-coverage)
   * [Scoring scheme for alignment](#scoring-scheme-for-alignment)
   * [Statistical significance of alignment](#statistical-significance-of-alignment)
   * [All-against-all pairwise alignments](#all-against-all-pairwise-alignments)
* [Tests](#tests)
* [License](#license)


## Introduction
Pairwise sequence alignment is used to identify regions of similarity that may indicate functional, structural and/or evolutionary relationships between two sequences.

1. Global alignment (*needle*; [the Needleman-Wunsch algorithm](https://en.wikipedia.org/wiki/Needleman–Wunsch_algorithm)) aligns two sequences across their entire length, from beginning to end. It is most useful when sequences you are aligning are similar and roughly the same size.
2. Local alignment (*water*; [the Smith-Waterman algorithm](https://en.wikipedia.org/wiki/Smith–Waterman_algorithm)) finds the region with the highest level of similarity between the two sequences. It is suitable for sequences that are not assumed to be similar over the entire length.


## Requirements

1. Python >= 3.8
    > Check with `python3 --version`
2. EMBOSS >= 6.6.0
    > Check with `needle -version` or `water -version`.


## Quick Start

```python
import pairwise_alignment as pa

# Global alignment
aln = pa.needle(moltype='nucl', qseq='ATGCTAGTA', sseq='ATGCTAGTAGATGATGA')
aln = pa.needle(moltype='prot', qseq='MKSTVWSG', sseq='MKSSVLW')

# Local alignment
aln = pa.water(moltype='nucl', qseq='ATGCTAGTA', sseq='ATGCTAGTAGATGATGAT')
aln = pa.water(moltype='prot', qseq='MKSTVWSG', sseq='MKSSVLW')

print(aln.score)       # 20.0
print(aln.pidentity)   # 71.4
print(aln.psimilarity) # 85.7
print(aln.pgaps)       # 14.3
print(aln.qseq)        # MKSTV-W
print(aln.sseq)        # MKSSVLW
```

## Alignment object

### Attributes

| Attribute | Description |
| --- | --- |
| `qid` | Query sequence identifier |
| `sid` | Subject sequence identifier |
| `qseq` | Query unaligned sequence |
| `sseq` | Subject unaligned sequence |
| `qaln` | Query aligned sequence 
| `saln` | Subject aligned sequence |
| `qstart` | Start of alignment in query |
| `qend` | End of alignment in query |
| `sstart` | Start of alignment in subject |
| `send` | End of alignment in subject |
| `length` | Alignment length |
| `score` | Alignment score |
| `nidentity` | Number of identical matches in the alignment |
| `pidentity` | Percentage of identical matches in the alignment |
| `nsimilarity` | Number of positive-scoring matches in the alignment |
| `psimilarity` | Percentage of positive-scoring matches in the alignment |
| `ngaps` | Total number of gaps in the alignment |
| `pgaps` | Total percentage of gaps in the alignment |
| `moltype` | nucl/prot |
| `program` | needle/water |
| `gapopen` | Gap open penalty |
| `gapextend` | Gap extension penalty |
| `matrix` | Name of scoring matrix |  
| `raw` | Raw output obtained from EMBOSS' needle/water |


### Methods

| Method | Description |
| --- | --- |
| `query_coverage()` | Returns a query coverage [%] |
| `subject_coverage()` | Returns a subject coverage [%] |
| `pvalue()` | Returns a p-value of the alignment |
| `fasta()` | Returns pairwise alignment in FASTA/Pearson format |


## Usage examples

### Alignment information

```python
import pairwise_alignment as pa

aln = pa.needle(
    moltype='prot',
    qseq='MTSPSTKNSDDKGRPNLSSTEYFANTNVLTCRLKWVNPDTFIMDPRKPQLHSRT',
    sseq='MTTPSRENSDDKGRPIEEASNLSSTEYFANTNVLTCKLKYVNPDTFIMDPRKP',
    qid='seq1',
    sid='seq2'
)

print(aln.score)        # 225.0
print(aln.length)       # 59
print(aln.pidentity)    # 72.9
print(aln.psimilarity)  # 79.7
print(aln.pgaps)        # 18.6
print(aln.ngaps)        # 11
print(aln.qseq)         # MTSPSTKNSDDKGRP-----NLSSTEYFANTNVLTCRLKWVNPDTFIMDPRKPQLHSRT
print(aln.sseq)         # MTTPSRENSDDKGRPIEEASNLSSTEYFANTNVLTCKLKYVNPDTFIMDPRKP------
print(aln.qstart)       # 1
print(aln.qend)         # 54
print(aln.sstart)       # 1
print(aln.send)         # 53
print(aln.program)      # needle
print(aln.gapopen)      # 10
print(aln.gapextend)    # 0.5
print(aln.matrix)       # EBLOSUM62
```

### Alignment in raw format

```python
print(aln.raw)
```

Output:

```
#=======================================
#
# Aligned_sequences: 2
# 1: seq1
# 2: seq2
# Matrix: EBLOSUM62
# Gap_penalty: 10.0
# Extend_penalty: 0.5
#
# Length: 59
# Identity:      43/59 (72.9%)
# Similarity:    47/59 (79.7%)
# Gaps:          11/59 (18.6%)
# Score: 225.0
# 
#
#=======================================

seq1               1 MTSPSTKNSDDKGRP-----NLSSTEYFANTNVLTCRLKWVNPDTFIMDP     45
                     ||:||.:||||||||     ||||||||||||||||:||:||||||||||
seq2               1 MTTPSRENSDDKGRPIEEASNLSSTEYFANTNVLTCKLKYVNPDTFIMDP     50

seq1              46 RKPQLHSRT     54
                     |||      
seq2              51 RKP------     53


#---------------------------------------
#---------------------------------------
```

### Alignment in FASTA format

```python
print(aln.fasta(wrap=30))
```

Output:

```
>seq1 1-54
MTSPSTKNSDDKGRP-----NLSSTEYFAN
TNVLTCRLKWVNPDTFIMDPRKPQLHSRT
>seq2 1-53
MTTPSRENSDDKGRPIEEASNLSSTEYFAN
TNVLTCKLKYVNPDTFIMDPRKP------
```

### Alignment iteration

You can iterate over the alignment residue-by-residude:

```python
for aa1, aa2 in aln:
    if aa1 != aa2:
        print(aa1, aa2)
```

Output:

```
S T
T R
K E
- I
- E
- E
- A
- S
R K
W Y
Q -
L -
H -
S -
R -
T -
```

### Alignment iteration with index

You can also loop through residudes in an alignment with the information of their position.

```python
for i, (aa1, aa2) in enumerate(aln):
    if aa1 != aa2:
        print(i, aa1, aa2)
```

Output:

```
2  S T
5  T R
6  K E
15 - I
16 - E
17 - E
18 - A
19 - S
36 R K
39 W Y
53 Q -
54 L -
55 H -
56 S -
57 R -
58 T -
```

### Query coverage

Query coverage describes how much of the query sequence is covered in the alignment by the subject sequence. Specifically, query coverage is the percentage of the query sequence length that is included in the alignment. In global alignments, query coverage is always 100% because both the sequences, query and subject, are aligned from end to end. It is thus more useful to calculate query coverage from local alignments.

```python
import pairwise_alignment as pa

aln = pa.water(
    moltype='prot',
    qseq='MTSPSTKNSDDKGRPNLSSTEYFANTNVLTCRLKWVNPDTFIMDPRKPQLHSRT',
    sseq='NSDDKGRPIEEASNLSSTEYFANTNVLTCKLKYVNPDTFIMDPRKP',
    qid='seq1',
    sid='seq2'
)
# seq1               8 NSDDKGRP-----NLSSTEYFANTNVLTCRLKWVNPDTFIMDPRKP     48
#                      ||||||||     ||||||||||||||||:||:|||||||||||||
# seq2               1 NSDDKGRPIEEASNLSSTEYFANTNVLTCKLKYVNPDTFIMDPRKP     46

print(aln.query_coverage())
# 75.9
print(aln.subject_coverage())
# 100
```

### Scoring scheme for alignment

You can change a scoring matrix and penalties for the gap open and extension to calculate the alignment.

```python
import pairwise_alignment as pa

aln = pa.water(
    moltype='prot',
    qseq='MKSTWYERNST',
    sseq='MKSTGYWTRESA',
    matrix='EBLOSUM30',
    gapopen=5,
    gapextend=0.2
)
print(aln)
```

Output:

```
#=======================================
#
# Aligned_sequences: 2
# 1: query
# 2: subject
# Matrix: EBLOSUM30
# Gap_penalty: 5.0
# Extend_penalty: 0.2
#
# Length: 13
# Identity:       7/13 (53.8%)
# Similarity:     8/13 (61.5%)
# Gaps:           3/13 (23.1%)
# Score: 39.8
# 
#
#=======================================

query              1 MKST--WYERNST     11
                     ||||  |. |.|:
subject            1 MKSTGYWT-RESA     12


#---------------------------------------
#---------------------------------------
```

### Statistical significance of alignment

The Needleman-Wunsch and Smith-Waterman algorithms will always find an optimal alignment between two sequences, whether or not they are evolutionarily related. The strength of an alignment is determined by its score. However, often it is necessary to know if a score is high enough to indicate a biologically interesting alignment. The statistical significance of the score is assessed by the *P*-value, which describes how likely it is that two random sequences of similar length and composition will align with a score equal to or better than our target alignment. 

The `.pvalue()` method calculates the *P*-value of the alignment between query and subject sequences. The method shuffles a subject sequence many times (100 by default) and calculates the alignment score between the query and each shuffled subject sequence. It then counts how many times the alignment score was greater than or equal to the alignment score of the original query and subject sequences. For example, if 100 such shuffles all produce alignment scores that are lower than the observed alignment score, then one can say that the *P*-value is likely to be less than 0.01.

```python
import pairwise_alignment as pa

aln = pa.needle(moltype='prot', qseq='MKSTVILK', sseq='MKSRSLK')

print(aln.pvalue())   # 0.16
```


### All-against-all pairwise alignments

For more than two sequences, you can calculate alignments between every pair of your input sequences.

```python
import itertools
import pairwise_alignment as pa

# Input sequences
sequences = {
    'dna1': 'ATCGAGATCGAGATGGCGATAG',
    'dna2': 'ATGCTGATCGTAGGGGC',
    'dna3': 'GTCGGATCCTCGATGGAGA',
    'dna4': 'TTTGGGAATGCGTAGGAGCTA',
    'dna5': 'CCGTGATGCGATGCA'
}

# All-against-all pairwise alignments
for qid, sid in itertools.combinations(sequences, r=2):
    qseq = sequences[qid]
    sseq = sequences[sid]
    aln = pa.needle(moltype='nucl', qseq=qseq, sseq=sseq)
    print(f'{qid} {sid} {aln.pidentity:.1f}% {aln.score}')
```

Output:

```
dna1  dna2  38.5%  24.0
dna1  dna3  60.9%  34.0
dna1  dna4  35.7%  20.0
dna1  dna5  45.8%  27.0
dna2  dna3  43.5%  18.0
dna2  dna4  42.3%  39.0
dna2  dna5  39.1%  22.0
dna3  dna4  20.0%  14.0
dna3  dna5  52.4%  26.0
dna4  dna5  40.9%  14.0
```

### All-against-all pairwise alignments of sequences from a FASTA file

If you have multiple sequences in a FASTA file, you can use [Biopython](https://biopython.org) to read them and then calculate pairwise alignments.

```python
import itertools
import pairwise_alignment as pa

from Bio import SeqIO

# Input sequences
sequences = {}
for seq_record in SeqIO.parse('sequences.fasta', 'fasta'):
    sequences[seq_record.id] = str(seq_record.seq)

# All-against-all pairwise alignments
for qid, sid in itertools.combinations(sequences, r=2):
    qseq = sequences[qid]
    sseq = sequences[sid]
    aln = pa.needle(moltype='nucl', qseq=qseq, sseq=sseq)
    print(f'{qid} {sid} {aln.pidentity:.1f}% {aln.score}')
```

Output:

```
dna1  dna2  38.5%  24.0
dna1  dna3  60.9%  34.0
dna1  dna4  35.7%  20.0
dna1  dna5  45.8%  27.0
dna2  dna3  43.5%  18.0
dna2  dna4  42.3%  39.0
dna2  dna5  39.1%  22.0
dna3  dna4  20.0%  14.0
dna3  dna5  52.4%  26.0
dna4  dna5  40.9%  14.0
```


## Tests
This module contains automated tests. If you want to check that everything works as intended, just run:

```
python3 pairwise_alignment.py
```

## License

[GNU General Public License, version 3](https://www.gnu.org/licenses/gpl-3.0.html)