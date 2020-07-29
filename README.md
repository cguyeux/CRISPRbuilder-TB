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
and for the GNU/Linux, Mac and Windows platforms, but the preferred method is
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


## Produced results

- A file with ".contig" is written in the SRA directory, which contains all 
obtained contigs ordered by scores. In the best case scenario, only the 2-3 first
lines are needed to reconstruct the CRISPR locus.


## Deeper explanations

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
</pre>


## Citation

>Christophe Guyeux, Christophe Sola, Guislaine Refrégier. CRISPRbuilder-TB: “CRISPR-Builder for tuberculosis”. Exhaustive reconstruction of the CRISPR locus in Mycobacterium tuberculosis complex using SRA. Submitted article (2020).
