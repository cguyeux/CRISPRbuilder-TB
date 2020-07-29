import argparse
import Bio.SeqIO
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
        self._dir_data = pathlib.Path() / 'data'
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
        parser.add_argument("-graph_method",
                            type=self.__str2bool,
                            nargs = '?',
                            default=False,
                            const=True,
                            help="use the experimental graph method")
        parser.add_argument("-size_tuple",
                            default=3,
                            help="tuple size in graph method",
                            type=int)
        args = parser.parse_args()
        self._outdir = pathlib.Path(args.output_directory)
        self._sra = self._outdir / args.sra
        self._evalue = args.evalue
        self._num_threads = args.num_threads
        self._graph_method = args.graph_method
        self._taille_tuple = args.size_tuple

    def _make_blast_db(self):
        completed = sp.run([self._makeblastdb,
                            '-in', self._sra_shuffled,
                            '-dbtype', 'nucl',
                            '-title', self._sra.name,
                            '-out', self._sra / self._sra.name])
        assert completed.returncode == 0

    def parse(self):
        self._sequences_of_interest()
        self._get_matches()
        if self._graph_method:
            print(self._get_contigs2())
        else:
            self._get_contigs1()
        self._find_duplicates()

    def _get_contigs1(self):
        # TODO: Code cleaning
        kmers = int(4 * self._len_reads / 5)
        print('Read length:', self._len_reads)
        print('k-mers length:', kmers)
        sequences = []
        for dd in self._SEQS:
            sequences.extend([dd[u:u + kmers] for u in range(len(dd) - kmers)])
        total = []
        SEQUENCES = []
        while len(sequences) > 0:
            cible = sequences.pop(0)
            nb_reads = 1 + sequences.count(cible)
            sequences = list(filter(lambda a: a != cible, sequences))
            avance = True
            while avance:
                prochain = [k for k in sequences if k[:-1] == cible[-len(k) + 1:]]
                U = [[u[-1] for u in prochain].count(v) for v in 'ACGT']
                avance = ((sorted(U)[-1] >= 2 * sorted(U)[-2]) and (sorted(U)[-1] > 1))
                if avance:
                    nb_reads += sorted(U)[-1]
                    cible += 'ACGT'[U.index(max(U))]
                    for k in prochain:
                        sequences.remove(k)
            recule = True
            while recule:
                prochain = [k for k in sequences if k[1:] == cible[:len(k) - 1]]
                U = [[u[0] for u in prochain].count(v) for v in 'ACGT']
                recule = ((sorted(U)[-1] >= 2 * sorted(U)[-2]) and (sorted(U)[-1] > 1))
                if recule:
                    nb_reads += sorted(U)[-1]
                    cible = 'ACGT'[U.index(max(U))] + cible
                    for k in prochain:
                        sequences.remove(k)
            SEQUENCES.append(cible)
            u = copy.deepcopy(cible)
            for k in self._dicofind:
                u = u.replace(self._dicofind[k], '*' + k + '*')
            # u = u.replace('GTCGTCAGACCCAAAACCC', 'rDRa1').replace('AAAACCCCGAGAGGGGACGGAAAC', 'DRb2').replace('CCCCGAGAGGGGACGGAAAC','DRb1')
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
        with open(self._sra / (self._sra.name + '.contig'), 'w') as f:
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

    def _get_chaine(self, G, node):
        chaine = node.split('*')[0]
        successor = list(G.neighbors(node))
        suite = ''
        score = 0
        while len(successor) == 1:
            score += G[node][successor[0]]['weight']
            chaine += '*' + successor[0].split('*')[0]
            suite = '*'.join(successor[0].split('*')[1:])
            node = successor[0]
            successor = list(G.neighbors(successor[0]))
        if suite != '':
            chaine += '*' + suite
        return (chaine, successor, score)

    def _formate_contigs(self, chaine):
        chaine = chaine.replace('AAAACCCCGAGAGGGGACGGAAAC', '*DRb2')
        chaine = chaine.replace('motif_fin1*motif_fin2*motif_fin3', 'motif_fin')
        chaine = chaine.replace('**', '*')
        cpt, txt = 0, ''
        for k in chaine.split('*'):
            txt += k + '*'
            cpt += len(k) + 1
            if cpt > 75:
                cpt = 0
                txt += '\n'
        return txt

    def _find_duplicates(self):
        couples = []
        chevauche = 12
        for k in self._SEQS:
            for l in ['DR' + str(m) for m in range(0, 16)]:
                if self._dicofind[l] in k:
                    for i in range(k.count(self._dicofind[l])):
                        u, v = k.split(self._dicofind[l])[i], k.split(self._dicofind[l])[i + 1]
                        U, V = '', ''
                        for m in self._dicofind:
                            if u[-chevauche:] == self._dicofind[m][-chevauche:]:
                                U = m
                            if v[:chevauche] == self._dicofind[m][:chevauche]:
                                V = m
                        couples.append((U, V))
            kk = rev_comp(k)
            for l in ['DR' + str(m) for m in range(0, 16)]:
                if self._dicofind[l] in kk:
                    for i in range(kk.count(self._dicofind[l])):
                        u, v = kk.split(self._dicofind[l])[i], kk.split(self._dicofind[l])[i + 1]
                        U, V = '', ''
                        for m in self._dicofind:
                            if u[-chevauche:] == self._dicofind[m][-chevauche:]:
                                U = m
                            if v[:chevauche] == self._dicofind[m][:chevauche]:
                                V = m
                        couples.append((U, V))
        limite = 3
        for k in sorted([(u, couples.count(u)) for u in list(set(couples)) if
                         couples.count(u) >= limite and '' not in u and 'pattern' not in ''.join(u)],
                        key=lambda x: eval(x[0][0].replace('esp', '').split('(')[0])):
            if eval(k[0][1].replace('esp', '').replace('rev(', '').split('(')[0].split(')')[0]) != eval(
                    k[0][0].replace('esp', '').replace('rev(', '').split('(')[0].split(')')[0]) + 1:
                print(k)

    def _get_contigs2(self, sensibility=20):
        liste = [u for u in [[k[l:l + self._taille_tuple]
                              for l in range(len(k) - (self._taille_tuple - 1))]
                             for k in S
                             if len(k) >= self._taille_tuple - 1] if u != []]
        triplets = [item for sublist in liste for item in sublist]

        occurrences_minimales = 3
        sommets = [item for sublist in triplets for item in sublist]
        sommets = sorted(set([k for k in sommets if sommets.count(k) >= occurrences_minimales]))

        self._triplets = [k for k in triplets if all([u in sommets for u in k])]
        arcs = {}
        for k in self._triplets:
            for l in self._triplets:
                if k[1:] == l[:self._taille_tuple - 1]:
                    cle = '*'.join(k) + '+' + '*'.join(l)
                    if cle not in arcs:
                        arcs[cle] = 1
                    else:
                        arcs[cle] += 1
        G = nx.DiGraph()
        for arc in arcs:
            k, l = arc.split('+')
            if arcs[arc] >= sensibility:
                G.add_edge(k, l, weight=arcs[arc])
        write_dot(G, self._sra / (self._sra.name + '.dot'))
        investigated = []
        roots = collections.deque([k for k in G.nodes() if list(G.predecessors(k)) == []])
        contigs = []
        while len(roots) != 0:
            node = roots.popleft()
            investigated.append(node)
            chaine, suite, score = self._get_chaine(G, node)
            if suite != []:
                roots.extend([u for u in suite if u not in investigated])
                for k in suite:
                    if chaine + k not in contigs:
                        contigs.append((chaine + '*' + k, score))
            else:
                if chaine not in [u[0] for u in contigs]:
                    contigs.append((chaine, score))
        if len(contigs) == 2:
            if contigs[0][0].endswith('IS6110') and contigs[1][0].startswith('finIS6110'):
                txt = contigs[0][0] + contigs[1][0]
                contigs = [(txt.replace('IS6110finIS6110', 'IS6110'), contigs[0][1] + contigs[1][1])]
            elif contigs[1][0].endswith('IS6110') and contigs[0][0].startswith('finIS6110'):
                txt = contigs[1][0] + contigs[0][0]
                contigs = [(txt.replace('IS6110finIS6110', 'IS6110'), contigs[0][1] + contigs[1][1])]
        return contigs

    def _get_matches(self):
        self._dicofind, self._dico_cas = {}, {}
        with open(pathlib.Path('data') / 'fastas' / 'crispr_patterns.fasta') as f:
            txt = f.read()
        for k in txt.split('>')[1:]:
            self._dicofind[k.split('\n')[0]] = k.split('\n')[1]
        with open(pathlib.Path('data') / 'fastas' / 'cas_patterns.fasta') as f:
            txt = f.read()
        for k in txt.split('>')[1:]:
            self._dico_cas[k.split('\n')[0]] = k.split('\n')[1]
        S = []
        for k in self._SEQS:
            s = copy.deepcopy(k)
            for l in self._dicofind:
                s = s.replace(self._dicofind[l], '*' + l + '*')
            for l in self._dico_cas:
                s = s.replace(self._dico_cas[l], '*' + l + '*')
            s = s.replace('**', '*')
            s = s.replace('*GTCGTCAGACCCAAAACCC*', '*rDRa1*')
            s = s.replace('*CCCCGAGAGGGGACGGAAAC*', '*DRb1*')
            s = s.replace('*AAAACCCCGAGAGGGGACGGAAAC*', '*DRb2*')
            if '*' in s:
                S.append(s.split('*')[1:-1])

if __name__ == '__main__':
    cb = CRISPRbuilder()
    cb.parse()