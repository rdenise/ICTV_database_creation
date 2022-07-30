#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import logging
from tqdm import tqdm
from pathlib import Path
import Bio
from Bio import Entrez, SeqIO
import re

##########################################################################################

Entrez.email = "rdenise@ucc.ie"
Entrez.tool = "create_ictv.py"
Entrez.api_key = "c7271a66422f38febdbfbc6169ae65764108"


##########################################################################################
##########################################################################################
##
##                                Functions
##
##########################################################################################
##########################################################################################


def efetch_accession2gbk(accGenBank_nameFile):
    """
    Function that will get the assembly file from the ncbi server

    :params url: the url of the file (e.g. "https://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank.txt")
    :type: str
    :return: A dataframe with the information of the assembly file
    :rtype: pandas.dataframe
    """

    accGenBank = accGenBank_nameFile["Virus GENBANK accession"]
    nameFile = accGenBank_nameFile["File_identifier"]

    gbk_file = GenBank / f"{nameFile}.gbk"

    if gbk_file.is_file():
        logging.debug(f"-> Reading: {nameFile}.gbk")

        gbk = SeqIO.read(gbk_file, "genbank")

        with open(gbk_file, "rt") as r_file:
            line = r_file.readline()
            bp_file = int(re.search(r"([0-9]+) bp", line).group(1))

        if len(gbk) == bp_file:
            return

    logging.debug(f"-> Downloading the identifier {accGenBank}")

    handle = Entrez.efetch(
        db="nucleotide", rettype="gbwithparts", retmode="text", id=accGenBank
    )
    gbk = SeqIO.read(handle, "genbank")

    logging.debug(f"-> Creating: {nameFile}.gbk")

    SeqIO.write(gbk, gbk_file, format="genbank")

    # global counter_gbk
    # counter_gbk.increment()

    return


##########################################################################################
