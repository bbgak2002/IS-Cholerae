 
import logging
from os import system, listdir
from tbtools.ISprofiler import ISprofiler

SRAS = open('data/SRR_Acc_List.txt').read().split('\n')

for sra in SRAS:
  system(f'fastq-dump --split-files --fasta 0 -A {sra}')
  system(f'mkdir {sra}')
  system(f'mv {sra}*.fasta {sra}')
 
  profiler = ISprofiler(sra=sra,
                        outdir='',
                        nb_threads=4,
                        loglevel=logging.DEBUG)
  profiler.parse()
