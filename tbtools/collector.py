#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 18 07:27:17 2020

@author: Christophe Guyeux
"""

import Bio.Entrez
import logging
import logging.config
import os.path
import pathlib
import subprocess as sp
import time
import xmltodict
from os import system

from .tools import check_for_tools, read_dict, write_dict


class Collector:

    def __init__(self, sras=[], bioprojects=[], biosamples=[], outdir="sequences",
                 loglevel=logging.INFO):
        logging.basicConfig(level=loglevel)
        self._logger = logging.getLogger()
        self._logger.info('Initialization')
        self._dico_info = {}
        self.sras = sras
        self._bioprojects = bioprojects
        self._biosamples = biosamples
        self._outdir = pathlib.Path(outdir)
        self._fastq_dump, _, _ = check_for_tools()

    def _set_sras(self):
        self._logger.debug('Collecting all SRA accession numbers')
        self._add_list_to_sras()
        Bio.Entrez.email = ""
        self._add_bioprojects_to_sras()
        self._add_bio_samples_to_sras()
        self._add_taxid_to_sras()

    def _add_list_to_sras(self):
        # TODO: _add_list_to_sras
        self._logger.debug('Adding SRA accessions from list')
        pass

    def _add_bioprojects_to_sras(self):
        self._logger.debug('Adding SRA accessions from bioprojects')
        for bioproject in self._bioprojects:
            self._logger.debug(f"Getting SRA IDs from bioproject {bioproject}")
            ret = Bio.Entrez.esearch(db="sra",
                                     term=bioproject,
                                     retmode="xml")
            dico = xmltodict.parse(ret.read())
            id_sras = dico['eSearchResult']['IdList']['Id']
            if isinstance(id_sras, str):
                id_sras = [id_sras]
            for item in id_sras:
                self._logger.debug(f"Getting SRA accession number from {item} ID")
                ret = Bio.Entrez.efetch(db="sra", id=item,
                                        retmode="xml", retmax=10000)
                dico = xmltodict.parse(ret.read())
                sra = dico['EXPERIMENT_PACKAGE_SET']['EXPERIMENT_PACKAGE']['RUN_SET']['RUN']['@accession']
                self.sras.append(sra)
                time.sleep(1)

    def _add_bio_samples_to_sras(self):
        self._logger.debug('Adding SRA accessions from biosamples')
        for biosample in self._biosamples:
            self._logger.debug(f"Getting SRA IDs from biosample {biosample}")
            ret = Bio.Entrez.esearch(db="sra",
                                     term=biosample,
                                     retmode="xml")
            dico = xmltodict.parse(ret.read())
            id_sras = dico['eSearchResult']['IdList']['Id']
            if isinstance(id_sras, str):
                id_sras = [id_sras]
            for item in id_sras:
                self._logger.debug(f"Getting SRA accession number from {item} ID")
                ret = Bio.Entrez.efetch(db="sra", id=item,
                                        retmode="xml", retmax=10000)
                dico = xmltodict.parse(ret.read())
                sra = dico['EXPERIMENT_PACKAGE_SET']['EXPERIMENT_PACKAGE']['RUN_SET']['RUN']['@accession']
                self.sras.append(sra)
                time.sleep(1)

    def _add_taxid_to_sras(self):
        self._logger.debug('Adding SRA accessions from taxids')
        # TODO: _add_taxid_to_sras
        pass

    def _prepare_directory(self):
        """
        Prepare the structure of sequence directory
        """
        self._logger.info(f'Preparing directory structure for {self._sra}')
        self._dir = self._outdir / self._sra
        self._dir.mkdir(exist_ok=True, parents=True)
        #self._dico_info[self._sra]['directory'] = self._dir

    def _collect_and_prepare_sra(self):
        files = [k for k in os.listdir(self._dir)
                 if k.startswith(self._sra) and k.endswith('.fasta')]
        if not files:
            system(
                f'fastq-dump --split-files --fasta 0 -O {self._dir} {self._sra}') #prefetch --fasta -O {self._dir / self._sra} {self._sra} &&
            '''
            self._logger.info('Downloading SRA fasta files')
            completed = sp.run([self._fastq_dump,
                                '--split-files',
                                '--fasta',
                                '-O', self._dir,
                                self._sra
                                ])
            while completed.returncode != 0:
                completed = sp.run([self._fastq_dump,
                                    '--split-files',
                                    '--fasta',
                                    '-O', self._dir,
                                    self._sra
                                    ])'''
        else:
            self._logger.info('SRA fasta files already downloaded')

    def get(self):
        self._set_sras()
        for sra in self.sras:
            self._sra = sra
            self._dico_info[sra] = read_dict(self._outdir / self._sra)
            self._prepare_directory()
            self._collect_and_prepare_sra()
            write_dict(self._dico_info[sra])
            
            
            
            
            
            
            
