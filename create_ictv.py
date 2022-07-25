#!/usr/bin/env python3
# -*- coding: utf-8 -*-

##########################################################################################
##########################################################################################
##
##                                Library
##
##########################################################################################
##########################################################################################

import argparse
from textwrap import dedent
import sys, os
import pandas as pd
import requests
from io import StringIO, BytesIO
import re
import hashlib
import gzip
import Bio
from Bio import SeqIO  
from Bio import Seq 
import logging
from tqdm import tqdm

##########################################################################################
##########################################################################################
##
##                                Functions
##
##########################################################################################
##########################################################################################


def create_folder(mypath):

    """
    Created the folder that I need to store my result if it doesn't exist
    :param mypath: path where I want the folder (write at the end of the path)
    :type: string
    :return: Nothing
    """

    try:
        os.makedirs(mypath)
    except OSError:
        pass

    return


##########################################################################################
##########################################################################################


def get_assembly_file(url) :
    """
    Function that will get the assembly file from the ncbi server

    :params url: the url of the file (e.g. "https://ftp.ncbi.nlm.nih.gov/genomes/genbank/viral/assembly_summary.txt")
    :type: str
    :return: A dataframe with the information of the assembly file
    :rtype: pandas.dataframe
    """

    response = requests.get(url)

    assembly_summary = pd.read_table(StringIO(response.text), skiprows=1)

    logging.debug("Finish to read the assembly file : {assembly}".format(assembly=url))

    assembly_summary.rename(columns={"# assembly_accession":"assembly_accession"}, inplace=True)

    return assembly_summary


##########################################################################################
##########################################################################################


def select_assembly_rows(assembly_df, kingdom) :
    """
    Function that will select the interresting columns of the assembly file

    :params assembly_df: Dataframe that contains the informations of the assembly ncbi file
    :type: pandas.dataframe
    :params kingdom: the name of the kingdom of the species of the assemnly file
    :type: str
    :return: A sub dataframe with the selected columns and the kingdom added at the end
    :rtype: pandas.dataframe
    """

    assembly_df = assembly_df[(assembly_df.genome_rep == 'Full') & (assembly_df.assembly_level == 'Complete Genome')].reset_index(drop=True)

    assembly_df.loc[:,"ftp_file"] = assembly_df.ftp_path.apply(lambda x : "{}_genomic.gbff.gz".format(x.split('/')[-1]))
    assembly_df.loc[:,"kingdom"] = kingdom

    logging.debug('New columns "ftp_file" and "kingdom" created for {kingdom}'.format(kingdom=kingdom))

    return assembly_df[['assembly_accession', 'species_taxid', 'ftp_path', 'ftp_file', 'organism_name', 'kingdom']]

##########################################################################################
##########################################################################################

def fetch_genbank_file(species, Genomes) :

    '''
    Function that will fetch the concatenated and gzipped genbank file and create the id 
    for the strain number curr_No_taxidsp (the id of the species)

    :params species: row of the assembly dataframe
    :type: pandas.Series
    :params Genomes: Path to the Genomes folder
    :type: str
    '''

    gbff_url = "{}/{}".format(species.ftp_path, species.ftp_file).replace('ftp:', 'https:')
    
    # print("-> Using or downloading/using {}".format(gbff_url))
    loging.debug("-> Using or downloading/using {}".format(gbff_url))

    gbff_response = requests.get(gbff_url)
    md5_gbff = hashlib.md5(gbff_response.content)

    # check integrity of the file after download (note: if file was dezipped and rezipped it will not be recognized anymore)
    md5_gbk = pd.read_csv(StringIO(requests.get(
                        os.path.join(species.ftp_path.replace('ftp:', 'https:'),'md5checksums.txt')).text), 
                        names=['md5','assembly_files'], sep='\s+')

    if md5_gbk[md5_gbk.assembly_files.str.contains(species.ftp_file)].md5.values == md5_gbff.hexdigest() :
        
        # print("\n-> md5 CHECKED OK")
        logging.debug('MD5 OK and checked for -> {genome}'.formar(genome=gbff_url))

    else :
        logging.debug('Did not check, erasing the file -> {genome}'.formar(genome=gbff_url))

    return 


##########################################################################################
##########################################################################################
##
##                                Main
##
##########################################################################################
##########################################################################################

parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
     description=dedent("""Gembases Microbial Creation""") )

general_option = parser.add_argument_group(title = "General input dataset options")
general_option.add_argument("-o",'--output',
                            default=None,
                            dest="output",
                            metavar='<FOLDER>',
                            help="Folder where to put the database (default: $PWD)")
general_option.add_argument("-d",'--date_stamp',
                            metavar="<DATE>",
                            dest="date_stamp",
                            help="The time stamp: MMYY (e.g. 1016 for october 2016)",
                            required=True)                            
parser.add_argument("-v", "--verbosity", 
                            default=1,
                            action="count",
                            dest="verbosity",
                            help="Increase log verbosity could be : -v or -vv (default: -v)")

args = parser.parse_args()

##########################################################################################

if args.output :
    OUTPUT = args.output
else :
    OUTPUT = os.getcwd()

##########################################################################################

if args.verbosity == 1 :
    level = logging.INFO
else : 
    level = logging.DEBUG

LOG_FORMAT = '%(levelname)s::%(asctime)s - %(message)s'
logging.basicConfig(filename = os.path.join(INPUT_DIRNAME, 'ictv_downloading.log'),
                    level = level,
                    format = LOG_FORMAT,
                    filemode = 'w')

logger = logging.getLogger()

##########################################################################################

taxa = os.path.join(OUTPUT, f"ICTV_database_{args.date_stamp}")
Genomes = os.path.join(taxa, "Genomes")
GenBank = os.path.join(taxa, "GenBank")
Genes = os.path.join(taxa, "Genes")
Proteins = os.path.join(taxa, "Proteins")
Gff = os.path.join(taxa, "Gff")
Assembly_report = os.path.join(taxa, "Assembly_Reports")


create_folder(taxa)
create_folder(Genomes)
create_folder(GenBank)
create_folder(Genes)
create_folder(Proteins)
create_folder(Gff)
create_folder(Assembly_report)

##########################################################################################

if args.verbosity == 1 :
    level = logging.INFO
elif args.verbosity == 2 :
    level = logging.DEBUG
else : 
    level = logging.NOTSET

LOG_FORMAT = '%(levelname)s::%(asctime)s - %(message)s'
logging.basicConfig(filename = os.path.join(taxa, 'do_gembase.log'),
                    level = level,
                    format = LOG_FORMAT,
                    filemode = 'w')

logger = logging.getLogger()

logging.info(f'Gembase creation logging for version : ICTV_database_{args.date_stamp}')

##########################################################################################

list_viral = "viral_assembly_summary.txt"
viral_table_assembly = "https://ftp.ncbi.nlm.nih.gov/genomes/genbank/viral/assembly_summary.txt"
assembly_summary_viral_file = os.path.join(taxa, list_viral)

##########################################################################################

logging.info(f"-> Fetching {viral_table_assembly}")
print(f"-> Fetching {viral_table_assembly}")

# Test if the file exist before to not have problem if the assembly file of NCBI change between two run for the same MICROBIAL
if os.path.isfile(assembly_summary_viral_file) :
    assembly_summary_viral.read_table(os.path.join(taxa, list_bacteria))
else :
    assembly_summary_viral = get_assembly_file(viral_table_assembly)
    # Writing the assembly summary
    assembly_summary_viral.to_csv(os.path.join(taxa, list_viral), index=False, sep='\t')

    logging.debug(f'{assembly_summary_viral_file} had been written')
 
logging.info('Done!')
print('\nDone!\n')

##########################################################################################

logging.info('-> Creating all the files for each genomes in assembly summary')
print('-> Creating all the files for each genomes in assembly summary')


stats_df = assembly_summary_modify.apply(fetch_genbank_file, axis=1)

stats_df.sort_values('Mnemo',inplace=True)
stats_df.to_csv(STATS_file, index=False, header=header_stats_write, sep='\t')

logging.info('Done!')
print('\nDone!\n')

if args.concat_prot :
    logging.info('-> Concatenating all the files for each genomes')
    print('-> Concatenating all the files for each genomes')

    all_prt = glob.glob(os.path.join(Proteins, '*prt'))

    num_file = len(all_prt)

    # Check number if files is less than 26
    FINAL_NUM_FILE = args.concat_prot[1] if  1 < args.concat_prot[1] < 26 else 1

    num_per_file = num_file // FINAL_NUM_FILE

    if FINAL_NUM_FILE != 1 :
        for index in range(FINAL_NUM_FILE) :
            begin = index * num_per_file
            end = (index + 1) * num_per_file if (index + 1) < 26 else -1

            file_prot = os.path.join(Genomes, '{}.{}.{}.prot'.format(args.concat_prot[0], chr(ord('@')+number),VERSION_ID))

            with open(file_prot, 'w') as w_file:
                for prt_file in all_prt[begin:end] :
                    seqs = SeqIO.parse(prt_file, 'fasta')
                    SeqIO.write(seqs, w_file, 'fasta')
    else :
        file_prot = os.path.join(Genomes, '{}.{}.prot'.format(args.concat_prot[0], VERSION_ID))

        with open(file_prot, 'w') as w_file:
            for prt_file in all_prt :
                seqs = SeqIO.parse(prt_file, 'fasta')
                SeqIO.write(seqs, w_file, 'fasta') 

    logging.info('Concatenation -> Done!')
    print('\nDone!\n')
##########################################################################################
