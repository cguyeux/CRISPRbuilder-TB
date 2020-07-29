import Bio.SeqIO
import distutils.spawn
import os
import pathlib
import subprocess as sp
import sys
from xlrd import open_workbook


def check_for_tools():
    for name in ['blastn', 'fastq-dump', 'makeblastdb']:
        if distutils.spawn.find_executable(f"{name}.exe") is not None:
            if name == 'blastn':
                blastn = 'blastn.exe'
            elif name == 'fastq-dump':
                fastqdump = 'fastq-dump.exe'
            else:
                makeblastdb = 'makeblastdb.exe'
        elif distutils.spawn.find_executable(name) is not None:
            if name == 'blastn':
                blastn = 'blastn'
            elif name == 'fastq-dump':
                fastqdump = 'fastq-dump'
            else:
                makeblastdb = 'makeblastdb'
        else:
            if sys.platform == 'win32':
                if name == 'blastn':
                    blastn = pathlib.Path('bin') / 'windows' / 'blastn.exe'
                elif name == 'fastq-dump':
                    fastqdump = pathlib.Path('bin') / 'windows' / 'fastq-dump.exe'
                else:
                    makeblastdb = pathlib.Path('bin') / 'windows' / 'makeblastdb.exe'
            elif sys.platform == 'linux':
                if name == 'blastn':
                    blastn = 'bin/linux/./blastn'
                elif name == 'fastq-dump':
                    fastqdump = 'bin/linux/./fastq-dump'
                else:
                    makeblastdb = 'bin/linux/./makeblastdb'
            elif sys.platform in ['darwin']:
                if name == 'blastn':
                    blastn = 'bin/mac/./blastn'
                elif name == 'fastq-dump':
                    fastqdump = 'bin/mac/./fastq-dump'
                else:
                    makeblastdb = 'bin/mac/./makeblastdb'
    return fastqdump, makeblastdb, blastn


def prepare_sra(sra):
    sra_shuffled = sra.parent / sra.name / (sra.name + '_shuffled')
    sra_shuffled = sra_shuffled.with_suffix('.fasta')
    if not pathlib.Path.is_file(sra_shuffled):
        if sys.platform in ['linux', 'darwin']:
            files = [k for k in sra.iterdir()
                     if k.name.startswith(sra.name)
                     and k.suffix == '.fasta'
                     and 'shuffled' not in k.stem]
            for fic in files:
                sp.run(["sed",
                        "-i", f"s/{sra.name}./{fic.stem}./g",
                        fic])
            with open(sra_shuffled, "w+") as f:
                sp.call(["cat", *files], stdout=f)
        else:
            # TODO: In windows, only the first sra is considered
            # TODO: log the process
            sra_shuffled = sra
    return sra_shuffled

def h37Rv():
    dir_data = pathlib.Path('data')
    h37Rv = Bio.SeqIO.read(dir_data / 'fastas' / 'NC_000962.3.fasta', 'fasta')
    if not any([p.suffix == '.nin' for p in dir_data.iterdir()]):
        _, makeblastdb, _ = check_for_tools()
        completed = sp.run([makeblastdb,
                            '-in', dir_data / 'fastas' / 'NC_000962.3.fasta',
                            '-dbtype', 'nucl',
                            '-title', 'h37Rv',
                            '-out', dir_data / 'h37Rv'])
        assert completed.returncode == 0
    return h37Rv


def seq_info(sra):
    sra_shuffled = prepare_sra(sra)
    if sys.platform in ['linux', 'darwin']:
        proc1 = sp.Popen(["cat", sra_shuffled],
                         stdout=sp.PIPE)
        proc2 = sp.Popen(["grep", ">"],
                         stdin=proc1.stdout,
                         stdout=sp.PIPE)
        proc3 = sp.run(["wc", "-l"],
                       stdin=proc2.stdout,
                       stdout=sp.PIPE)
        nb_reads = int(proc3.stdout.decode('utf8'))
    else:
        with open(sra_shuffled) as f:
            nb_reads = f.read().count('>')
    # Getting read length
    head_seq = open(sra_shuffled).read(10000).split('>')[1]
    len_reads = len(''.join(head_seq.splitlines()[1:]))
    # Getting reads coverture
    exact_cov = nb_reads * len_reads / len(h37Rv())
    coverage = round(exact_cov, 2)
    return coverage, len_reads, nb_reads


def str_to_spol(ch):
    liste = [k for k in range(1,99) if 'esp'+str(k)+'*' in ch or 'esp'+str(k)+'(' in ch 
             or 'esp'+str(k)+')' in ch or 'esp'+str(k)+'[' in ch]
    print(liste)
    spol,spol_old = '',''
    for k in range(1,99):
        if k in liste:
            spol += '■'
        else:
            spol += '□'
        if k in list(range(2,5))+list(range(12,16))+list(range(18,45))+list(range(46,48))+list(range(51,54))+list(range(62,66)):
            if k in liste:
                spol_old += '■'
            else:
                spol_old += '□'
    print(spol_old)
    print(spol)
    if spol_old in spol_sit:
        print('SIT :',spol_sit[spol_old])
    else:
        print('SIT : ND')


wb = open_workbook('data/SIT.xls')
ws = wb.sheet_by_index(0)
spol_sit = {}
for row in range(1, ws.nrows):
    spol, sit = ws.cell_value(row, 2).replace('n','\u25A0').replace('o','\u25A1'), ws.cell_value(row, 8)
    if spol not in spol_sit:
        spol_sit[spol] = sit
