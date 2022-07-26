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
import hashlib
import gzip
import logging
from tqdm import tqdm
import multiprocessing

##########################################################################################

session = None
data_to_be_processed = [...]

def init_process():
    global session
    session = requests.Session()

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


def write_file(gbff_response, local_filename) :

    '''
    Function that will read the (multi) genbank file and extract the infornations
    about each replicon.

    :params species: row of the assembly dataframe
    :type: pandas.Series
    :params Genomes: Path to the Genomes folder
    :type: str
    :params gbff_response: The gz file from the ftp of assembly database
    :type: requests.models.Response
    '''

    with open(local_filename, 'wt') as w_file:      
        w_file.write(gbff_response.text)

    return

##########################################################################################


def write_file_uncompress(gbff_response, local_filename) :

    '''
    Function that will read the (multi) genbank file and extract the infornations
    about each replicon.

    :params species: row of the assembly dataframe
    :type: pandas.Series
    :params Genomes: Path to the Genomes folder
    :type: str
    :params gbff_response: The gz file from the ftp of assembly database
    :type: requests.models.Response
    '''

    with gzip.open(BytesIO(gbff_response.content), mode='rt') as r_file:  
        with open(local_filename, 'wt') as w_file:      
            for line in r_file:        
                w_file.write(line)

    return


##########################################################################################


def fetch_genbank_file(species) :

    '''
    Function that will fetch the concatenated and gzipped genbank file and create the id 
    for the strain number curr_No_taxidsp (the id of the species)

    :params species: row of the assembly dataframe
    :type: pandas.Series
    :params Genomes: Path to the Genomes folder
    :type: str
    '''

    global session

    # check integrity of the file after download (note: if file was dezipped and rezipped it will not be recognized anymore)
    md5_gbk = pd.read_csv(StringIO(session.get(
                        os.path.join(species['ftp_path'].replace('ftp:', 'https:'),'md5checksums.txt')).text), 
                        names=['md5','assembly_files'], sep='\s+')

    ftp_location = {'ftp_gbff':GenBank,
                     'ftp_fna':Genomes,
                     'ftp_faa':Proteins, 
                     'ftp_gff':Gff, 
                     'ftp_gene':Genes, 
                     'ftp_report':Assembly_report,
                     }

    for ftp_file, local_folder in ftp_location.items():

        gbff_url = "{}/{}".format(species['ftp_path'], species[ftp_file]).replace('ftp:', 'https:')
        
        logging.debug(f"-> Using or downloading/using {gbff_url}")

        gbff_response = session.get(gbff_url)
        md5_gbff = hashlib.md5(gbff_response.content)

        sub_md5_gbff = md5_gbk[md5_gbk.assembly_files.str.contains(species[ftp_file])]

        if not sub_md5_gbff.empty and sub_md5_gbff.md5.values[0] == md5_gbff.hexdigest() :
            
            # print("\n-> md5 CHECKED OK")
            logging.debug(f'MD5 OK and checked for -> {gbff_url}')


            file_name = os.path.join(local_folder, species[ftp_file].replace('.gz',''))

            if species[ftp_file].endswith(".gz"):
                write_file_uncompress(gbff_response, file_name)
            else:
                write_file(gbff_response, file_name)

        else :
            logging.debug(f'Did not check, erasing the file -> {gbff_url}')

    return 


##########################################################################################
##########################################################################################
##
##                                Main
##
##########################################################################################
##########################################################################################

parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
     description=dedent("""ICTV Download and Creation""") )

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
general_option.add_argument(
    "-t",
    "--threads",
    metavar="<num_threads>",
    dest="threads",
    help="Number of threads to use (default:1)",
    required=True,
    default=1,
    type=int,
    choices=range(1, multiprocessing.cpu_count()),
)

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
logging.basicConfig(filename = os.path.join(OUTPUT, 'ictv_downloading.log'),
                    level = level,
                    format = LOG_FORMAT,
                    filemode = 'w')

logger = logging.getLogger()

##########################################################################################

taxa = os.path.join(OUTPUT, "ICTV_database", f"{args.date_stamp}")
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
    assembly_summary_viral = pd.read_table(assembly_summary_viral_file)
else :
    assembly_summary_viral = get_assembly_file(viral_table_assembly)
    # Writing the assembly summary

    assembly_summary_viral['ftp_gbff'] = assembly_summary_viral.ftp_path.apply(lambda x : "{}_genomic.gbff.gz".format(x.split('/')[-1]))
    assembly_summary_viral['ftp_fna'] = assembly_summary_viral.ftp_path.apply(lambda x : "{}_genomic.fna.gz".format(x.split('/')[-1]))
    assembly_summary_viral['ftp_faa'] = assembly_summary_viral.ftp_path.apply(lambda x : "{}_protein.faa.gz".format(x.split('/')[-1]))
    assembly_summary_viral['ftp_gff'] = assembly_summary_viral.ftp_path.apply(lambda x : "{}_genomic.gff.gz".format(x.split('/')[-1]))
    assembly_summary_viral['ftp_gene'] = assembly_summary_viral.ftp_path.apply(lambda x : "{}_cds_from_genomic.fna.gz".format(x.split('/')[-1]))
    assembly_summary_viral['ftp_report'] = assembly_summary_viral.ftp_path.apply(lambda x : "{}_assembly_report.txt".format(x.split('/')[-1]))

    assembly_summary_viral.to_csv(assembly_summary_viral_file, index=False, sep='\t')

    logging.debug(f'{assembly_summary_viral_file} had been written')
 
logging.info('Done!')
print('\nDone!\n')

##########################################################################################

logging.info('-> Creating all the files for each genomes in assembly summary')
print('-> Creating all the files for each genomes in assembly summary')

##### MULTIPROCESS ACTION
args_func = assembly_summary_viral.to_dict('records')

num_rows = assembly_summary_viral.shape[0]



pool = multiprocessing.Pool(args.threads, initializer=init_process)
results = list(
    tqdm(
        pool.imap(fetch_genbank_file, args_func), total=num_rows
    )
)
pool.close()


logging.info('Done!')
print('\nDone!\n')

##### TODO RENAME OF FILES + MERGE_REPORT

##########################################################################################
