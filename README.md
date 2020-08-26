# CRISPRbuilder-TB
Reconstruct *Mycobacterium tuberculosis* CRISPR locus from WGS data.

## Purpose of the package

**CRISPRbuilder-TB** reconstructs the whole CRISPR locus from a 
*Mycobacterium tuberculosis* complex genome, by using Whole Genome Sequencing (WGS) 
data. This allows to deduce the real spoligotype, to detect mutants of spacers and
direct repeats, insertion of mobile elements, and duplications. The Cas locus 
is reconstructed too. This is a semi automatic approach that leads to a set of 
contigs to assemble manually. Depending on the number, length, and quality of
SRAs, the number of contigs can range from 1-2 patterns, in the best case scenario
where the good quality of sequences allows a quasi-automatic reconstruction of the 
CRISPR cut in mobile element positions, to several contigs difficult to process,
for too short or polluted reads.


## Requirements

**CRISPRbuilder-TB** needs Python 3, and the following dependencies to work:
*blastn* and *makeblastdb* from [blast+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download) NCBI package.

[fastq-dump](https://github.com/ncbi/sra-tools) from sra-tools (NCBI) can be used too, to collect new WGS fasta files.

A 'bin' directory has been added to the project with the 3 required executables, 
and for the GNU/Linux, Mac and Windows platforms. But the preferred method is
a clean installation from the source.


## Installation

You can either download a zip file by clicking on the 'Clone or download' green
button, or clone the repository by writing in a terminal:
<pre>
git clone https://github.com/cguyeux/TB-tools.git
</pre>

To install the required python libraries: 
<pre>
pip3 install -r requirements.txt
</pre>

You can try to replace pip3 by pip, if the command does not succeed. And you may
add the --user option in case where you don't have the root privilege.

## How to use TB-tools

To launch a first analysis: 
<pre>
python crisprbuilder.py -sra SRR1173284
</pre>

- By default, SRAs are looked for in the *sequences* directory.
- If this directory does not contain a subdirectory of this name, read files are
then downloaded by using fastq-dump.

Help about this script: 
<pre>
python crisprbuilder.py -h
</pre>


## Illustrative examples

### A first easy example

Here are the first contigs produced in the cas of a BCG strain:
<pre>
('TTGCACGTATCGACCGCTTCGCCATCGACTGCGACAACATCCGGATCTATAAGATAAGAGGTGTTGCGGCAGTTACGTTCTACGGAAGGGGACG*Rv2816c*starting_pattern1*DR0*esp1*DR0*esp2*DR0*esp3*DR0*esp5*DR0*esp6*DR0*esp7*DR0*esp8*DR0*esp11*DR0*esp12*DR0*esp13*DR0*esp14*DR0*esp15*DR0*esp16*DR0*esp17*DR0*esp18*DR0*esp20*DR0*esp21*DR0*esp22*DR0*esp23*DR0*esp24*DR0*esp25*DRb2*esp27*DR0*esp28*DR0*esp29*DR0*esp30*DR2*esp31*DR0*esp32*DR0*esp33*DR0*esp34*rDRa1*IS6110*GGTCATGTCAGGTGGTTCATCGAGGAGGTACCCGCCGGAGCTGCGTGAGCGGGCGGTGCGGATGGTCGCAGAGATCCGCGGTCAGCACGATTCGGAGT', 42964)
('CGCCGCCTCTACCAGTACTGCGGCGACGTCCCGCCGGTCGAACTCGAGGCTGCCTACTACGCTCAACGCCAGAGACCAGCCGCCGGCTGAGGTCTC*finIS6110*DRb1*esp35*DR0*esp36*DR0*esp37*DR0*esp38*DR0*esp39*DR0*esp40*DR0*esp41*DR0*esp35*DR0*esp42*DR0*esp43*DR0*', 14706)
('*DR0*esp48*DR0*esp50*DR0*esp51*DR0*esp52*DR0*ending_pattern1*ending_pattern2*ending_pattern3*Rv2813c*TGCGGGTGGTGGATTCGTCGACGATGGCCTTGTCGGCGGCGAAGGCGGCGACGAGGGCTTGCAGGGCG', 9009)
('DR0[9:]*esp43*DR0*esp44*DR0*esp46*DR15*esp47*DR0*esp48*DR0[:29]', 4391)
('DR0[9:]*esp43*DR0*esp43*DR0[:27]', 1145)
('DR0[8:]*esp48*DR0*esp48*DR0[:29]', 669)
('*esp6*DR29*esp7*DR0[:21]', 52)
('AAC*esp37*DR15*esp47*GTCGT', 48)
('DR0[30:]*esp5*DR0*esp48*DR0[:6]', 42)
('esp46[20:]*DR0*esp34*rDRa1*IS6110[:22]', 38)
...
</pre>
- We can see a gap in scores (the right component of the couples) after the 6 
first lines, meaning that what follows can be forgotten;
- the first contig is the beginning of the locus, until the ancestral IS6110;
- the second one corresponds to the locus part after this IS, starting with
spacer 35, until spacer 43;
- we can find the part from 43 to 48 in the fourth contig, and the part from 
spacers 48 to 52, and then the ending pattern, in the third contig; 
- the ruptures at position 43 and 48 are explained by the lines 5 and 6: they
are due to a tandem duplication of each of these spacers.

The reconstructed locus is then: 
<pre>
*Rv2816c*starting_pattern1*DR0*esp1*DR0*esp2*DR0*esp3*DR0*esp5*DR0*esp6*DR0*esp7*DR0*esp8*DR0*esp11*DR0*esp12*DR0*esp13*DR0*esp14*DR0*esp15*DR0*esp16*DR0*esp17*DR0*esp18*DR0*esp20*DR0*esp21*DR0*esp22*DR0*esp23*DR0*esp24*DR0*esp25*DRb2*esp27*DR0*esp28*DR0*esp29*DR0*esp30*DR2*esp31*DR0*esp32*DR0*esp33*DR0*esp34*rDRa1*IS6110*DRb1*esp35*DR0*esp36*DR0*esp37*DR0*esp38*DR0*esp39*DR0*esp40*DR0*esp41*DR0*esp35*DR0*esp42*DR0*esp43*DR0*esp43*DR0*esp44*DR0*esp46*DR15*esp47*DR0*esp48*DR0*esp48*DR0*esp50*DR0*esp51*DR0*esp52*DR0*ending_pattern*Rv2813c*
</pre>

and we can recover its spoligotype and SIT number as follows:
<pre>
In [1]: from tools import *
In [2]: str_to_spol('*Rv2816c*starting_pattern1*DR0*esp1*DR0*esp2*DR0*esp3*DR0*esp5*DR0*esp6*DR0*esp7*DR0*esp8*DR0*esp1
   ...: 1*DR0*esp12*DR0*esp13*DR0*esp14*DR0*esp15*DR0*esp16*DR0*esp17*DR0*esp18*DR0*esp20*DR0*esp21*DR0*esp22*DR0*esp23
   ...: *DR0*esp24*DR0*esp25*DRb2*esp27*DR0*esp28*DR0*esp29*DR0*esp30*DR2*esp31*DR0*esp32*DR0*esp33*DR0*esp34*rDRa1*IS6
   ...: 110*DRb1*esp35*DR0*esp36*DR0*esp37*DR0*esp38*DR0*esp39*DR0*esp40*DR0*esp41*DR0*esp35*DR0*esp42*DR0*esp43*DR0*es
   ...: p43*DR0*esp44*DR0*esp46*DR15*esp47*DR0*esp48*DR0*esp48*DR0*esp50*DR0*esp51*DR0*esp52*DR0*ending_pattern*Rv2813c
   ...: *')
[1, 2, 3, 5, 6, 7, 8, 11, 12, 13, 14, 15, 16, 17, 18, 20, 21, 22, 23, 24, 25, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 46, 47, 48, 50, 51, 52]
■■□■■■■■□■■■■■■□■■■■■■■■■■■■■■■■■■■■■■□□□□□
■■■□■■■■□□■■■■■■■■□■■■■■■□■■■■■■■■■■■■■■■■■■□■■■□■■■□□□□□□□□□□□□□□□□□□□□□□□□□□□□□□□□□□□□□□□□□□□□□□
SIT : 482

In [3]:
</pre>

### A more difficult one

The case of ERR028615 is more complicated, as the contigs are multiple:

<pre>
('CAGATCTGACGGAAAC*esp65*DR0*esp66*DR4*esp67*DR5*esp68*DR0*ending_pattern1*ending_pattern2*ending_pattern3[:29]', 5269)
('esp45[12:17]*DR0*esp53*DR0*esp62*DR0*esp63*DR0*esp64*DR6*esp65*DR0[:20]', 4347)
('esp28[31:]*DR0*esp29*DR0*esp38*DR0*esp39*DR0*esp40[:9]', 2576)
('esp24[30:]*DR0*esp25*DRb2*esp26*DR0*esp27*DR0*esp28[:19]', 2454)
('esp2[28:]*DR0*esp3*DR0*esp4*DR0*esp12*DR0[:24]', 2691)
('TGCGTGTTCAA*DR0*esp51*DR0*esp52*DR0*esp53[:20]', 1587)
('esp40[20:]*DR0*esp41*DR0*esp35*DR0*esp42[:21]', 1600)
('GGTTCTTTTGA*starting_pattern1*DR0*esp1*DR0*esp2[:9]', 1940)
('GGGAAAC*esp22*DR0*esp23*DR0*esp24[:18]', 1207)
('G*DR0*esp21*DR0*esp22*DR0[:20]', 1125)
('esp23[31:]*DR0*esp24*DR0*esp25[:16]', 715)
('TGATGACTCCG*DR0*esp28*DR0*esp29[:8]', 837)
('esp1[30:]*DR0*esp2*DR0*esp3[:13]', 882)
('esp14[30:]*DR0*esp15*DR0*esp18[:11]', 840)
('GACCCTGTCA*DR0*esp42*DR0*ATTT', 884)
('esp13[21:]*DR0*esp14*DR0*esp15[:9]', 893)
('esp12[25:]*DR0*esp13*DR0*esp14[:12]', 722)
('esp18[28:]*DR0*esp19*DR0*esp20[:14]', 712)
('TGC*DR0*esp40*DR0*esp41[:9]', 630)
('esp19[21:]*DR0*esp20*DR0*esp21[:9]', 577)
('DRb1*esp18*DR0*ACTCG', 240)
('ACAGTTACGTTCTACGGAAGGGGACG*Rv2816c*starting_pattern1*DR0[:8]', 965)
('CAGATCTGAAC*esp12*DR0*esp13[:11]', 305)
('ATGA*ending_pattern3*Rv2813c*TGCGGGTGGTGGATTCCAGATCT', 443)
('esp15[28:]*DR0*esp18*GTCG', 165)
('GCGGGTACCTCCTCGATGAACCACCTGACATGACCCCATCCTTTCC*finIS6110c*GA', 355)
('AACCACCTGACATGACCCCATCCTTTCC*finIS6110c*GATCGGTCATATCAAGTTTTGTCAGGAATGCGGGATTCGAAT', 452)
('CAGATCTGGGAAAC*esp22*rDRa1', 14)
('DR0[25:]*esp38*GTCGTCAGACCCAAAACAGATCT', 10)
('ending_pattern1[35:]*ending_pattern2*CCCCGCAGATCT', 10)
('esp15[13:]*DR0*GCGCAGATCT', 12)
('esp27[9:]*DR0*AGATCT', 10)
('CAGATCTGCT*DR0*esp63[:27]', 14)
('CAGATCTGTTCAGCCT*DR0*esp15[:22]', 24)
('CAGATCTGAACCCCGAGAGGGGACGGAAAC*esp13*DR0[:8]', 21)
('CAGATCTGGACATGACCCCATCCTTTCC*finIS6110c*GATCGGTCAT', 16)
('DR0[18:]*esp15*GTCGTCAGATCT', 9)
('DR0[15:]*esp38*GTCGTCACAGATCT', 12)
('CAGATCTGTCCTCCCA*DR5*esp68[:22]', 24)
('CAGATCTGCTTAG*DR0*esp14[:23]', 11)
('esp22[16:]*DR0*TTCGCAGATCT', 12)
('CAGATCTGGGACGGAAAC*esp62*DR0[:16]', 12)
('esp64[28:]*DR6*TGGACGCAGAATCGCACCGGCAGATCT', 12)
('CAGATCTGCCACCTGACATGACCCCATCCTTTCC*finIS6110c*GA', 24)
...
</pre>

However, a carefull study of the latter shows that all patterns overlap, except between spacers 42 and 51. And this problematic situation can be solved thanks to the other outputted files, leading to :

<pre>
CRISPR:
*Rv2816c*starting_pattern1*DR0*esp1*DR0*esp2*DR0*esp3*DR0*esp4*DR0*esp12*DR0*esp13*DR0*
esp14*DR0*esp15*DR0*esp18*DR0*esp19*DR0*esp20*DR0*esp21*DR0*esp22*DR0*esp23*DR0*
esp24*DR0*esp25*DRb2*esp26*DR0*esp27*DR0*esp28*DR0*esp29*DR0*esp38*DR0*esp39*DR0*
esp40*DR0*esp41*DR0*esp35*DR0*esp42*DR0*esp51*DR0*esp52*DR0*esp53*DR0*esp62*DR0*
esp63*DR0*esp64*DR6*esp65*DR0*esp66*DR4*esp67*DR5*esp68*DR0*ending_pattern*Rv2813c*

CAS:
Cas6-Csm1-Csm2-Csm3-Csm4-Csm5-Csm6-Cas1-Cas2-pattern_start*DR0
end_pattern-Rv2812c-Rv2811c-Rv2809c-Rv2808c-Rv2807c
</pre>


## Deeper explanations

In total, 6 files are created, which can be used to manually determine the locus constitution. They are detailed below, according to their suffix:
- .contig: as previously explained, contigs of the CRISPR locus;
- .not_consecutive and .reads_with_2_spacers: the number of reads that contain the end of spacer k, followed by a DR, followed by the beginning of spacer l (in .not_consecutive, l is not k+1); 
- .cas_locus: a tentative reconstruction of the Cas locus, based on the blast of beginning and end of Cas genes, plus part of such genes followed by a beginning (or end) of IS6110, based on a short list of known events;
- .cas_reads: number of reads for each event of the list used in .cas_locus;
- .is_around: investigation of possibly unknown IS6110 in the flanking regions of the CRISPR locus.

It is not rare that the produced results are too chaotic, to make it possible to proper reconstruct manually the locus. This may be due to a too short read length, to contamination, or to a bad quality of the reads. However, files above depend on the following parameters, whose modification may sometimes help to decipher the locus content:
- evalue: evalue when blasting spacers, DRs, etc.;
- kmer_size: kmer size in the reads decomposition (if 0: 4/5 of the reads length);
- overlap: spacer nucleotides required when looking reads with two pieces of spacers;
- limit: minimal number of reads that contains two pieces of spacers;
these parameters have been experimentally set at assumed values.


## Last comment

If you find new DR or spacer variant, or a new spacer related to *tuberculosis*, please send them with the SRA accession number for **CRISPRbuilder-TB** integration, to: *christophe* dot *guyeux* at *univ* minus *fcomte* dot *fr*.


## Citation

>Christophe Guyeux, Christophe Sola, Guislaine Refrégier. CRISPRbuilder-TB: “CRISPR-Builder for tuberculosis”. Exhaustive reconstruction of the CRISPR locus in Mycobacterium tuberculosis complex using SRA. Submitted article (2020).
