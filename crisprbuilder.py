import argparse
import Bio.SeqIO
import Bio.pairwise2
import collections
import copy
import logging
import networkx as nx
import os
import pathlib
import pydot
import sys
import subprocess as sp

from networkx.drawing.nx_pydot import write_dot

from tools import check_for_tools, prepare_sra, seq_info, rev_comp


class CRISPRbuilder:
    def __init__(self):
        logging.basicConfig(level=logging.DEBUG)
        self._logger = logging.getLogger()
        self._logger.info('Initialization')
        self._dir = pathlib.Path() / 'data'
        self._parse_args()
        self._fastq_dump, self._makeblastdb, self._blastn = check_for_tools()
        p = pathlib.Path(self._outdir)
        p.mkdir(exist_ok=True, parents=True)
        if self._sra.name not in [u.name for u in self._outdir.iterdir()]:
            self._logger.info(f'Downloading {self._sra.name} reads')
            p = pathlib.Path(self._sra)
            p.mkdir(exist_ok=True, parents=True)
            sp.run([self._fastq_dump,
                    '--split-files',
                    '--fasta',
                    '-O', self._sra,
                    self._sra.name
                    ])
        self._logger.info(f'Preparing {self._sra.name} sequences')
        self._sra_shuffled = prepare_sra(self._sra)
        self._coverage, self._len_reads, self._nb_reads = seq_info(self._sra)
        if not any([p.suffix == '.nin' for p in self._sra.iterdir()]):
            self._make_blast_db()

    def __str2bool(self, v):
        if isinstance(v, bool):
            return v
        if v.lower() in ('yes', 'true', 't', 'y', '1'):
            return True
        elif v.lower() in ('no', 'false', 'f', 'n', '0'):
            return False
        else:
            raise argparse.ArgumentTypeError('Boolean value expected.')

    def _parse_args(self):
        """
        Parses the arguments provided to the CRISPRbuilder-TB
        """
        parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        parser.add_argument("-sra",
                            help="accession number to deal with",
                            default='',
                            type=str)
        parser.add_argument("-out",
                            "--output_directory",
                            default='sequences',
                            help="directory of outputs",
                            type=str)
        parser.add_argument("-num_threads",
                            default='1',
                            help="number of threads",
                            type=str)
        parser.add_argument("-evalue",
                            default='1e-7',
                            help="evalue when blasting spacers, DRs, etc.",
                            type=str)
        parser.add_argument("-overlap",
                            default=12,
                            help="spacer nucleotides required when looking reads with two pieces of spacers",
                            type=int)
        parser.add_argument("-limit",
                            default=3,
                            help="minimal number of reads that contains two pieces of spacers",
                            type=int)
        args = parser.parse_args()
        self._outdir = pathlib.Path(args.output_directory)
        self._sra = self._outdir / args.sra
        self._evalue = args.evalue
        self._num_threads = args.num_threads
        self._overlap = args.overlap
        self._limit = args.limit

    def _make_blast_db(self):
        completed = sp.run([self._makeblastdb,
                            '-in', self._sra_shuffled,
                            '-dbtype', 'nucl',
                            '-title', self._sra.name,
                            '-out', self._sra / self._sra.name])
        assert completed.returncode == 0

    def parse(self):
        print("Investigations started...")
        self._sequences_of_interest()
        self._get_contigs()
        self._find_duplicates()
        self._cas_investigation()
        self._find_IS_around()

    def _get_contigs(self):
        # TODO: Code cleaning
        self._dicofind = {}
        with open(pathlib.Path('data') / 'fastas' / 'crispr_patterns.fasta') as f:
            txt = f.read()
        for k in txt.split('>')[1:]:
            self._dicofind[k.split('\n')[0]] = k.split('\n')[1]
        kmers = int(4 * self._len_reads / 5)
        print('Read length:', self._len_reads)
        print('k-mers length:', kmers)
        sequences = []
        for dd in self._SEQS:
            sequences.extend([dd[u:u + kmers] for u in range(len(dd) - kmers)])
        total = []
        SEQUENCES = []
        while len(sequences) > 0:
            target = sequences.pop(0)
            nb_reads = 1 + sequences.count(target)
            sequences = list(filter(lambda a: a != target, sequences))
            avance = True
            while avance:
                prochain = [k for k in sequences if k[:-1] == target[-len(k) + 1:]]
                U = [[u[-1] for u in prochain].count(v) for v in 'ACGT']
                avance = ((sorted(U)[-1] >= 2 * sorted(U)[-2]) and (sorted(U)[-1] > 1))
                if avance:
                    nb_reads += sorted(U)[-1]
                    target += 'ACGT'[U.index(max(U))]
                    for k in prochain:
                        sequences.remove(k)
            recule = True
            while recule:
                prochain = [k for k in sequences if k[1:] == target[:len(k) - 1]]
                U = [[u[0] for u in prochain].count(v) for v in 'ACGT']
                recule = ((sorted(U)[-1] >= 2 * sorted(U)[-2]) and (sorted(U)[-1] > 1))
                if recule:
                    nb_reads += sorted(U)[-1]
                    target = 'ACGT'[U.index(max(U))] + target
                    for k in prochain:
                        sequences.remove(k)
            SEQUENCES.append(target)
            u = copy.deepcopy(target)
            for k in self._dicofind:
                u = u.replace(self._dicofind[k], '*' + k + '*')
            u = u.replace('**', '*')
            v = []
            for k in u.split('*'):
                if set(list(k)).issubset(set('ACGT')):
                    if len(k) >= 6:
                        trouve = False
                        for l in self._dicofind:
                            if k in self._dicofind[l]:
                                debut = self._dicofind[l].index(k)
                                if debut == 0:
                                    chaine = l + '[:'
                                else:
                                    chaine = l + '[' + str(self._dicofind[l].index(k) + 1) + ':'
                                fin = self._dicofind[l].index(k) + len(k)
                                if fin == len(self._dicofind[l]):
                                    chaine += ']'
                                else:
                                    chaine += str(self._dicofind[l].index(k) + len(k)) + ']'
                                v.append(chaine)
                                trouve = True
                                break
                        if not trouve:
                            v.append(k)
                    else:
                        v.append(k)
                else:
                    v.append(k)
            u = '*'.join(v)
            u = u.replace('DR0[13:]', 'DRb2').replace('DR0[17:]', 'DRb1').replace('DR0[:19]', 'rDRa1').replace(
                'AACCGAGAGGGGACGGAAAC', 'DRb(1)')
            total.append((u, nb_reads))
        filename = self._sra / (self._sra.name + '.contig')
        print(f"Writing deduced contigs in {filename}")
        with open(filename, 'w') as f:
            for k in sorted([u for u in total if u[1] > 1], key=lambda x: x[0].count('*'), reverse=True):
                f.write(str(k)+os.linesep)

    def _sequences_of_interest(self):
        completed = sp.run([self._blastn,
                            '-num_threads', self._num_threads,
                            '-query', pathlib.Path('data') / 'fastas' / 'crispr_patterns.fasta',
                            '-evalue', self._evalue,
                            '-task', "blastn",
                            '-db', self._sra / self._sra.name,
                            '-outfmt', '10 sseqid sstart send',
                            '-out', self._sra / (self._sra.name + '_crispr_patterns.blast')])
        assert completed.returncode == 0
        seqs = {}
        with open(self._sra / (self._sra.name + '_crispr_patterns.blast'), 'r') as f:
            for u in f.read().split('\n')[:-1]:
                nom, deb,fin = u.split(',')
                deb,fin = eval(deb),eval(fin)
                seqs[nom] = deb<fin
        self._SEQS = []
        fasta_sequences = Bio.SeqIO.parse(open(self._sra_shuffled),'fasta')
        for fasta in fasta_sequences:
            if fasta.id in seqs:
                if seqs[fasta.id]:
                    self._SEQS.append(str(fasta.seq))
                else:
                    self._SEQS.append(str(fasta.seq.reverse_complement()))

    def _find_duplicates(self):
        couples = []
        for k in self._SEQS:
            for l in [m for m in self._dicofind if m.startswith('DR')]:
                if self._dicofind[l] in k:
                    for i in range(k.count(self._dicofind[l])):
                        u, v = k.split(self._dicofind[l])[i], k.split(self._dicofind[l])[i + 1]
                        U, V = '', ''
                        for m in self._dicofind:
                            if u[-self._overlap:] == self._dicofind[m][-self._overlap:]:
                                U = m
                            if v[:self._overlap] == self._dicofind[m][:self._overlap]:
                                V = m
                        couples.append((U, V))
            kk = rev_comp(k)
            for l in [m for m in self._dicofind if m.startswith('DR')]:
                if self._dicofind[l] in kk:
                    for i in range(kk.count(self._dicofind[l])):
                        u, v = kk.split(self._dicofind[l])[i], kk.split(self._dicofind[l])[i + 1]
                        U, V = '', ''
                        for m in self._dicofind:
                            if u[-self._overlap:] == self._dicofind[m][-self._overlap:]:
                                U = m
                            if v[:self._overlap] == self._dicofind[m][:self._overlap]:
                                V = m
                        couples.append((U, V))
        txt = ''
        for k in sorted([(u, couples.count(u)) for u in list(set(couples)) if
                         couples.count(u) >= self._limit and '' not in u and 'pattern' not in ''.join(u)
                         and 'IS' not in ''.join(u)],
                        key=lambda x: eval(x[0][0].replace('esp', '').split('(')[0])):
            if eval(k[0][1].replace('esp', '').replace('rev(', '').split('(')[0].split(')')[0]) != eval(
                    k[0][0].replace('esp', '').replace('rev(', '').split('(')[0].split(')')[0]) + 1:
                txt += f"{k[0][0].split('(')[0]}-{k[0][1].split('(')[0]}: {k[1]}" + os.linesep
        filename = self._sra / (self._sra.name + '.not_consecutive')
        print(f"Writing reads number with non-consecutive spacers in {filename}")
        with open(filename, 'w') as f:
            f.write(txt)
        txt = ''
        for k in sorted([(u, couples.count(u)) for u in list(set(couples)) if
                couples.count(u) >= self._limit and '' not in u and 'pattern' not in ''.join(u)
                         and 'IS' not in ''.join(u)],
               key=lambda x: eval(x[0][0].replace('esp', '').split('(')[0])):
            txt += f"{k[0][0].split('(')[0]}-{k[0][1].split('(')[0]}: {k[1]}" + os.linesep
        filename = self._sra / (self._sra.name + '.reads_with_2_spacers')
        print(f"Writing reads number with two spacers in {filename}")
        with open(filename, 'w') as f:
            f.write(txt)

    def _cas_investigation(self):
        self._dico_cas = {}
        with open(pathlib.Path('data') / 'fastas' / 'cas_patterns.fasta') as f:
            txt = f.read()
        for k in txt.split('>')[1:]:
            self._dico_cas[k.split('\n')[0]] = k.split('\n')[1]
        completed = sp.run([self._blastn,
                            '-num_threads', self._num_threads,
                            '-query', pathlib.Path('data') / 'fastas' / 'cas_patterns.fasta',
                            '-evalue', self._evalue,
                            '-task', "blastn",
                            '-db', self._sra / self._sra.name,
                            '-outfmt', '10 qseqid sseqid sstart send sseq',
                            '-out', self._sra / (self._sra.name + 'cas_patterns.blast')])
        assert completed.returncode == 0
        s = ''
        with open(self._sra / (self._sra.name + 'cas_patterns.blast'), 'r') as f:
            txt = f.read()
            for u in self._dico_cas:
                if txt.count(self._dico_cas[u]) > 0:
                    s += u + '-'
        s = s.replace('start_pattern-start_pattern', 'start_pattern')
        s = s.replace('DR0-DR0', 'DR0')
        for u in ['Cas6', 'Csm1', 'Csm2', 'Csm3', 'Csm4', 'Csm5', 'Csm6', 'Cas1', 'Cas2',
                  'pattern', 'Rv2812c', 'Rv2811c', 'Rv2809c', 'Rv2808c', 'Rv2807c']:
            for i in range(2, 8):
                s = s.replace(u + '_start' + str(i), u + '_start')
                s = s.replace(u + '_end' + str(i), u + '_end')
            s = s.replace(u + '_start-' + u + '_start', u + '_start')
            s = s.replace(u + '_end-' + u + '_end', u + '_end')
            s = s.replace(u + '_start-' + u + '_end', u)

        s = s.replace('end_pattern_start-end_pattern_middle-end_pattern_end', 'end_pattern')
        s = s.rstrip('-')
        s = s.replace('*DR0*esp1', '-CRISPR')
        filename = self._sra / (self._sra.name + '.cas_locus')
        print(f"Writing CAS locus in {filename}")
        with open(filename, 'w') as f:
            f.write(s)
        s = ''
        for u in self._dico_cas:
            if 'IS' in u:
                s += f"{txt.count(self._dico_cas[u])} {u}"
            else:
                s += f"  {txt.count(self._dico_cas[u])} {u}"
            s += os.linesep
        filename = self._sra / (self._sra.name + '.cas_reads')
        print(f"Writing nb reads per CAS pattern in {filename}")
        with open(filename, 'w') as f:
            f.write(s)

    def _find_IS_around(self):
        genes_cas = {}
        with open(pathlib.Path('data') / 'fastas' / 'cas_genes.fasta') as f:
            txt = f.read()
        for k in txt.split('>')[1:]:
            genes_cas[k.split('\n')[0]] = k.split('\n')[1]
        dico_IS = {}
        with open(pathlib.Path('data') / 'fastas' / 'IS6110.fasta') as f:
            txt = f.read()
        for k in txt.split('>')[1:]:
            dico_IS[k.split('\n')[0]] = k.split('\n')[1]
        with open(pathlib.Path('data') / 'fastas' / 'CASs.fasta', 'w') as f:
            for k in genes_cas:
                cpt, mot = 0, 50
                for l in range(0, len(genes_cas[k]), mot):
                    if len(genes_cas[k][l:]) >= mot / 2:
                        f.write(f">{k + str(cpt)}\n{genes_cas[k][l:l + mot]}\n")
                        cpt += 1
                for l in range(mot // 2, len(genes_cas[k]), mot):
                    if len(genes_cas[k][l:]) >= mot / 2:
                        f.write(f">{k + str(cpt)}\n{genes_cas[k][l:l + mot]}\n")
                        cpt += 1

        completed = sp.run([self._blastn,
                            '-num_threads', self._num_threads,
                            '-query', pathlib.Path('data') / 'fastas' / 'CASs.fasta',
                            '-evalue', self._evalue,
                            '-task', "blastn",
                            '-db', self._sra / self._sra.name,
                            '-outfmt', '10 sseqid sstart send',
                            '-out', self._sra / (self._sra.name + '_cas')])
        assert completed.returncode == 0

        seqs = {}
        with open( self._sra / (self._sra.name + '_cas'), 'r') as f:
            for u in f.read().split('\n')[:-1]:
                nom, deb, fin = u.split(',')
                deb, fin = eval(deb), eval(fin)
                seqs[nom] = deb < fin

        fasta_sequences = {}
        SEQS = []
        fasta_sequences = Bio.SeqIO.parse(self._sra_shuffled, 'fasta')
        for fasta in fasta_sequences:
            if fasta.id in seqs:
                if seqs[fasta.id]:
                    SEQS.append(str(fasta.seq))
                else:
                    SEQS.append(rev_comp(str(fasta.seq)))

        is_autour = {}
        txt = ''
        for k in SEQS:
            for l in dico_IS:
                debut = dico_IS[l][:15]
                if debut in k:
                    prefixe = k.split(debut)[0]
                    for m in genes_cas:
                        if prefixe in genes_cas[m]:
                            desc = f"{m}[:{genes_cas[m].find(prefixe) + len(prefixe)}]+{l}"
                            if desc not in is_autour:
                                is_autour[desc] = 1
                            else:
                                is_autour[desc] += 1
                fin = dico_IS[l][-15:]
                if fin in k:
                    suffixe = k.split(fin)[1]
                    for m in genes_cas:
                        if suffixe in genes_cas[m]:
                            desc = f"{l}+{m}[{genes_cas[m].find(suffixe) + len(suffixe)}:]"
                            if desc not in is_autour:
                                is_autour[desc] = 1
                            else:
                                is_autour[desc] += 1
        filename = self._sra / (self._sra.name + '.is_around')
        print(f"Writing ISs around the locus in {filename}")
        with open(filename, 'w') as f:
            f.write(os.linesep.join([f"{k:<40} : {is_autour[k]:>10} times" for k in is_autour]))


if __name__ == '__main__':
    cb = CRISPRbuilder()
    cb.parse()