import Bio.pairwise2
import Bio.SeqIO
import collections
import logging
import os.path
import pathlib
import pickle
import subprocess as sp

from .tools import check_for_tools, prepare_sra, seq_info, h37Rv, read_dict, write_dict


class ISprofiler:
    # TODO: Rename filename without capital letter
    def __init__(self, sra="", outdir="sequences", loglevel=logging.INFO,
                 IS_evalue='1e-7', IS_prefix_size=20,
                 percentage_of_coverage=0.05, prefix_similarity=0.85,
                 nb_threads=12):
        logging.basicConfig(level=loglevel)
        self._logger = logging.getLogger()
        self._logger.info('Initialization')
        self._dir_data = pathlib.Path() / 'data'
        self._outdir = outdir
        self._sra = pathlib.Path(self._outdir) / sra
        self._dico_info = read_dict(self._sra)
        self._IS_evalue = IS_evalue
        self._IS_prefix_size = IS_prefix_size
        self._percentage_of_coverage = percentage_of_coverage
        self._prefix_similarity = prefix_similarity
        self._num_threads = str(nb_threads)
        self._results = {}
        _, self._makeblastdb, self._blastn = check_for_tools()
        self._h37Rv = h37Rv()
        self._logger.info(f'Preparing {self._sra.name} sequences')
        self._sra_shuffled = prepare_sra(self._sra)
        completed = sp.run(['makeblastdb',
                            '-in', self._sra_shuffled,
                            '-dbtype', 'nucl',
                            '-title', self._sra.name,
                            '-out', self._sra / self._sra.name])
        self._coverage, self._len_reads, self._nb_reads = seq_info(self._sra)

    def __similarity(self, x, y):
        alignments = Bio.pairwise2.align.globalxx(x, y)
        return alignments[0][2]/alignments[0][4]

    def __get_position_in_h37Rv(self, name, seq):
        pos = 0
        with open(self._sra / 'snp.fasta','w') as f:
            f.write(f'>{name}\n{seq}')
        result = sp.run([self._blastn,
                         "-num_threads", self._num_threads,
                         "-query", self._sra / 'snp.fasta',
                         "-evalue", '1e-2',
                         "-task", "blastn",
                         "-db", self._dir_data / "h37Rv",
                         "-outfmt", "10 sstart send sseq"],
                        stdout=sp.PIPE)
        formated_results = result.stdout.decode('utf8').splitlines()
        if formated_results:
            pos = formated_results[0]
        return pos

    def _find_IS(self):
        p = pathlib.Path(self._dir_data) / "dict_IS.pkl"
        if os.path.isfile(p):
            with open(p, "rb") as f:
                dico = pickle.load(f)
        else:
            dico = {}
        p = pathlib.Path(self._dir_data) / 'IS'
        for IS in sorted(p.iterdir()):
            IS_name = IS.stem
            if IS_name not in dico:
                dico[IS_name] = {}
            self._logger.info(f"Looking for {IS_name}")
            if IS_name not in self._dico_info:
                self._dico_info[IS_name] = []
                completed = sp.run([self._blastn,
                                    '-num_threads', self._num_threads,
                                    '-query', IS,
                                    '-evalue', self._IS_evalue,
                                    '-task', "blastn",
                                    '-db', self._sra / self._sra.name,
                                    "-max_target_seqs", "2000000",
                                    '-outfmt',
                                    '10 qstart sstart send sseqid'],
                                   stdout=sp.PIPE)
                assert completed.returncode == 0
                matches = completed.stdout.decode('utf8').splitlines()
                good_start_matches = []
                for k in range(self._IS_prefix_size, self._len_reads):
                    good_start_matches.extend([u for u in matches
                                               if u.startswith('1,'+str(k)+',')])
                good_start_matches = [(u.split(',')[1], u.split(',')[3])
                                      for u in good_start_matches
                                      if eval(u.split(',')[1]) < eval(u.split(',')[2])]
                noisy_prefixes = []
                fasta_sequences = Bio.SeqIO.parse(self._sra_shuffled, 'fasta')
                if good_start_matches != []:
                    for fasta in fasta_sequences:
                        name, sequence = fasta.id, str(fasta.seq)
                        for match in good_start_matches:
                            if name == match[1]:
                                seq = sequence[eval(match[0])-self._IS_prefix_size:eval(match[0])-1]
                                noisy_prefixes.append(seq)
                counter = collections.Counter(noisy_prefixes)
                self._dico_info[f"{IS_name}_blasted"] = counter
                if 'coverage' not in self._dico_info:
                    self._dico_info['coverage']=0
                cleaned_prefixes = [k[0] for k in counter.items()
                                    if k[1] >= self._dico_info['coverage']*self._percentage_of_coverage]
                self._logger.info(f"{len(cleaned_prefixes)} copie(s) found")
                for seq in cleaned_prefixes:
                    formatted_sra_name = f"{self._sra.name}({counter[seq]})"
                    if seq in dico[IS_name]:
                        dico[IS_name][seq]['sra'].append(formatted_sra_name)
                        self._dico_info[IS_name].append(f"{seq}({dico[IS_name][seq]['position']})")
                    else:
                        for prefix in dico[IS_name]:
                            print(f"Prefix:{prefix}")
                            if seq in dico[IS_name][prefix]['variants']:
                                dico[IS_name][prefix]['sra'].append(formatted_sra_name)
                                self._dico_info[IS_name].append(f"{prefix}({dico[IS_name][prefix]['position']})")
                                break
                            elif self.__similarity(prefix, seq) > self._prefix_similarity:
                                dico[IS_name][prefix]['variants'].append(seq)
                                dico[IS_name][prefix]['sra'].append(formatted_sra_name)
                                self._dico_info[IS_name].append(f"{prefix}({dico[IS_name][prefix]['position']})")
                                break
                            elif [u for u in dico[IS_name][prefix]['variants']
                                  if self.__similarity(u,seq) > self._prefix_similarity]:
                                dico[IS_name][prefix]['variants'].append(seq)
                                dico[IS_name][prefix]['sra'].append(formatted_sra_name)
                                self._dico_info[IS_name].append(f"{prefix}({dico[IS_name][prefix]['position']})")
                                break
                            else:
                                ret = self.__get_position_in_h37Rv(IS_name, seq)
                            if ret == 0:
                                dico[IS_name][seq] = {
                                    'variants': [],
                                    'sens': '',
                                    'sra': [formatted_sra_name],
                                    'position': 0
                                }
                                self._dico_info[IS_name].append(f"{prefix}(0)")
                            else:
                                pos_deb, pos_fin, seq0 = ret.split(',')
                                if eval(pos_deb) < eval(pos_fin):
                                    sens = '-'
                                    pos = eval(pos_fin)+1
                                else:
                                    sens = '+'
                                    pos = eval(pos_deb)-1
                                found = [k for k in dico[IS_name]
                                         if pos == dico[IS_name][k]['position']]
                                if found:
                                    assert len(found) == 1
                                    dico[IS_name][found[0]]['variants'].append(seq)
                                    dico[IS_name][found[0]]['sra'].append(formatted_sra_name)
                                    self._dico_info[IS_name].append(f"{found[0]}({dico[IS_name][found[0]]['position']})")
                                else:
                                    dico[IS_name][seq0] = {
                                        'variants': [],
                                        'sens': sens,
                                        'sra': [formatted_sra_name],
                                        'position': pos
                                    }
                                    self._dico_info[IS_name].append(f"{seq0}({pos})")
                with open(pathlib.Path(self._dir_data) / "dict_IS.pkl", "wb") as f:
                    pickle.dump(dico, f)
                #self._dico_info.update(dico)
                write_dict(self._dico_info)
            else:
                self._logger.info(f"{self._dico_info[IS_name]}")

    def parse(self):
        self._logger.info('Looking for the insertion sequences')
        self._find_IS()
        
        
        
