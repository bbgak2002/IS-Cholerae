import Bio.SeqIO
import distutils.spawn
import os.path
import pathlib
import pickle
import subprocess as sp
import sys


def spol_octal_to_bin(dd):
    return ''.join([(('0'*3)+bin(int(k))[2:])[-3:] for k in dd[:-1]])+dd[-1]


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
    print(f"Locate sequences in {sra} directory")
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
    else:
        files = [k for k in sra.iterdir()
                 if k.name.startswith(sra.name)
                 and k.suffix == '.fasta'
                 and 'shuffled' not in k.stem]
        for fic in files:
            os.remove(fic)
    return sra_shuffled


def h37Rv():
    dir_data = pathlib.Path('data')
    h37Rv = Bio.SeqIO.read(dir_data / 'fastas' / 'NC_002505.1.fasta','fasta') #'NC_000962.3.fasta', 'fasta')
    if not any([p.suffix == '.nin' for p in dir_data.iterdir()]):
        _, makeblastdb, _ = check_for_tools()
        completed = sp.run([makeblastdb,
                            '-in', dir_data / 'fastas' / 'NC_002505.1.fasta',
                            '-dbtype', 'nucl',
                            '-title', 'h37Rv',                    #'h37Rv'
                            '-out', dir_data / 'h37Rv'])          #'h37Rv'
        assert completed.returncode == 0
    return h37Rv

def h37RvCDS():
    dir_data = pathlib.Path('data')
    h37Rv = dir_data / 'fastas' / 'NC_002505.1_CDS.fasta'        #'NC_000962.3_CDS.fasta
    data = open(h37Rv).read().split('locus_tag=')[1:]
    return {u.split(']')[0]: ''.join(u.split('\n')[1:]).split('>')[0] for u in data}

def mt18b():
    dir_data = pathlib.Path('data')
    mt18b = Bio.SeqIO.read(dir_data / 'fastas' / 'mt18b.fasta', 'fasta')
    # Todo, changer la ligne ci-dessous
    if not any([p.suffix == '.nin' for p in dir_data.iterdir()]):
        _, makeblastdb, _ = check_for_tools()
        completed = sp.run([makeblastdb,
                            '-in', dir_data / 'fastas' / 'mt18b.fasta',
                            '-dbtype', 'nucl',
                            '-title', 'mt18b',
                            '-out', dir_data / 'mt18b'])
        assert completed.returncode == 0
    return mt18b

def rev_comp(s):
    u=[change(x) for x in s]
    u.reverse()
    return ''.join(u)

def change(x):
    if x == 'A':
        return 'T'
    elif x == 'T':
        return 'A'
    elif x == 'C':
        return 'G'
    elif x == 'G':
        return 'C'
    return x

h37Rv_CDS = open(pathlib.Path('data') / 'fastas' / 'NC_002505.1_CDS.fasta').read().split('locus_tag=')[1:] #NC_000962.3_CDS.fasta
dico_genes = {}
for k in h37Rv_CDS:
    name = k.split(']')[0]
    location = k.split('location=')[1].split(']')[0]
    reverse = 'complement' in location
    if reverse:
    	debut = int(location.split('..')[0].replace('complement(', ''))
    	fin = int(location.split('..')[1].split(']')[0].replace(')', ''))
    else:
    	debut = int(location.split('..')[0])
    	fin = int(location.split('..')[1].split(']')[0])
    protein = k.split('protein=')[1].split(']')[0]
    seq = ''.join(k.split('\n')[1:])
    """if reverse:
        seq = rev_comp(''.join(k.split('\n')[1:]))
    else:
        seq = ''.join(k.split('\n')[1:])"""
    dico_genes[name] = {
        'debut': debut,
        'fin': fin,
        'reverse': reverse,
        'seq': seq,
        'protein': protein
    }



def get_gene_name(pos):
    return [k for k in dico_genes if dico_genes[k]['debut']<=pos<= dico_genes[k]['fin']]

def get_gene_func(pos):
    return [dico_genes[k]['protein'] for k in dico_genes if dico_genes[k]['debut']<=pos<= dico_genes[k]['fin']]




def seq_info(sra):
    print(f"Sra path is {sra}")
    sra_shuffled = prepare_sra(sra)
    for k in range(1,10):
        sra_shuffled = pathlib.Path(str(sra_shuffled).replace('sequences'+str(k), 'sequences'))
    print(sra_shuffled)
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


def write_dict(dico):
    rep = dico['directory']
    for k in range(1,10):
        rep = pathlib.Path(str(rep).replace('sequences'+str(k), 'sequences'))
    #rep = pathlib.Path('sequences') / rep.name
    with open( rep / (rep.name + '.pkl'), 'wb') as f:
        pickle.dump(dico, f)


def read_dict(dico, verbose=True):
    p = dico / (dico.name + '.pkl')
    if os.path.isfile(p):
        with open(p, 'rb') as f:
            dico_info = pickle.load(f)
        if verbose:
            if 'directory' in dico_info:
                print(f"Directory: {dico_info['directory']}")
            if 'name' in dico_info:
                print(f"Name: {dico_info['name']}")
            if 'strain' in dico_info:
                print(f"Strain: {dico_info['strain']}")
            if 'taxid' in dico_info:
                print(f"TaxID: {dico_info['taxid']}")
            if 'scientific_name' in dico_info:
                print(f"Scientific name: {dico_info['scientific_name']}")
            if 'biosample' in dico_info:
                print(f"Biosample: {dico_info['biosample']}")
            if 'biosample_title' in dico_info:
                print(f"Biosample title: {dico_info['biosample_title']}")
            if 'study' in dico_info:
                print(f"study: {dico_info['study']}")
            if 'bioproject' in dico_info:
                print(f"BioProject: {dico_info['bioproject']}")
            if 'date' in dico_info:
                print(f"Date of isolation: {dico_info['date']}")
            if 'location' in dico_info:
                print(f"Location of isolation: {dico_info['location']}")
            if 'coverage' in dico_info:
                print(f"Coverage: {dico_info['coverage']}")
            if 'nb_reads' in dico_info:
                print(f"Number of reads: {dico_info['nb_reads']}")
            if 'len_reads' in dico_info:
                print(f"Length of reads: {dico_info['len_reads']}")
            if 'lineage_Coll' in dico_info:
                if dico_info['lineage_Coll'] is not None:
                    print(f"Lineage Coll: {','.join(dico_info['lineage_Coll'])}")
            if 'lineage_Coll_full' in dico_info:
                print(f"Lineage Coll (full): {dico_info['lineage_Coll_full']}")
            if 'lineage_L6_animal' in dico_info:
                print(f"Lineage L6 animal: {dico_info['lineage_L6_animal']}")
            if 'lineage_PGG' in dico_info:
                print(f"Lineage PGG: {dico_info['lineage_PGG']}")
            if 'lineage_Palittapongarnpim' in dico_info:
                print(f"Lineage Palittapongarnpim: {dico_info['lineage_Palittapongarnpim']}")
            if 'lineage_Palittapongarnpim_full' in dico_info:
                print(f"Lineage Palittapongarnpim (full): {dico_info['lineage_Palittapongarnpim_full']}")
            if 'lineage_Shitikov' in dico_info:
                print(f"Lineage Shitikov: {dico_info['lineage_Shitikov']}")
            if 'lineage_Stucki' in dico_info:
                print(f"Lineage Stucki: {dico_info['lineage_Stucki']}")
            if 'spoligotype43_in_silico' in dico_info:
                print(f"Spoligotype 43 in silico: ")
                print(dico_info['spoligotype43_in_silico'])
            if 'spoligotype43_in_silico_count' in dico_info:
                print(f"Matches per spacer: {dico_info['spoligotype43_in_silico_count']}")
            if 'SIT' in dico_info:
                print(f"SIT: {dico_info['SIT']}")
            if 'spoligotype98_in_silico' in dico_info:
                print(f"Spoligotype 98 in silico: ")
                print(dico_info['spoligotype98_in_silico'])
            if 'spoligotype98_in_silico_count' in dico_info:
                print(f"Matches per spacer: {dico_info['spoligotype98_in_silico_count']}")
        return dico_info
    else:
        return {'directory': dico}
        
        
        
        
        
        
        
        
        
        
