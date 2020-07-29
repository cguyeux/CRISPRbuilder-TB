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

You need first the reads of a given Sequence Read Archive (SRA) accession number,
e.g., by downloading them with fastq-dump:
<pre>
fastq-dump --fasta --split-files SRR1173284
</pre>

To launch a first analysis: 
<pre>
python crisprbuilder.py SRR1173284
</pre>
Help about this script: 
<pre>
python crisprbuilder.py -h
</pre>

## Citation

>Christophe Guyeux, Christophe Sola, Guislaine Refrégier. CRISPRbuilder-TB: “CRISPR-Builder for tuberculosis”. Exhaustive reconstruction of the CRISPR locus in Mycobacterium tuberculosis complex using SRA. Submitted article (2020).
