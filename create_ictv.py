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
from unicodedata import name
from textwrap import dedent
import sys, os
import pandas as pd
import requests
from io import StringIO, BytesIO
import hashlib
import gzip
import logging
from logging.handlers import QueueHandler, QueueListener
from tqdm import tqdm
import multiprocessing
from pathlib import Path

##########################################################################################

session = None

##########################################################################################
##########################################################################################
##
##                                Functions
##
##########################################################################################
##########################################################################################

def logger_init():
    mpQueue = multiprocessing.Queue()

    LOG_FORMAT = '%(levelname)s::%(asctime)s - %(message)s'
    logging.basicConfig(filename = os.path.join(args.output, 'ictv_downloading.log'),
                        level = level,
                        format = LOG_FORMAT,
                        filemode = 'w')

    # this is the handler for all log records
    handler = logging.StreamHandler()
    LOG_FORMAT_HANDLER = "%(levelname)s: %(asctime)s - %(process)s - %(message)s"
    handler.setFormatter(logging.Formatter(LOG_FORMAT_HANDLER))

    # queueListerner gets records from the queue and sends them to the handler
    queueListerner = QueueListener(mpQueue, handler)
    queueListerner.start()

    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    # add the handler to the logger so records from this process are handled
    logger.addHandler(handler)

    return queueListerner, mpQueue

##########################################################################################

def init_process(mpQueue):
    global session
    session = requests.Session()

    # all records from worker processes go to queueHandler and then into mpQueue
    queueHandler = QueueHandler(mpQueue)
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    logger.addHandler(queueHandler)

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
    ftp_md5_path = f"{species['ftp_path'].replace('ftp:', 'https:')}/md5checksums.txt"

    try: 
        md5_gbk = pd.read_csv(StringIO(session.get(
                            ftp_md5_path).text), 
                            names=['md5','assembly_files'], sep='\s+')
    except:
        # Because some url seems to not be updated in the assembly file
        old_ftp = ftp_md5_path
        new_ftp = old_ftp.replace(old_ftp.split('_')[-1], species['asm_name'])

        logging.debug(f"Exception:: Change last url to correct {old_ftp} -> {new_ftp}")

        md5_gbk = pd.read_csv(StringIO(session.get(
                            new_ftp).text), 
                            names=['md5','assembly_files'], sep='\s+')        

    ftp_location = {
                'ftp_gbff':GenBank,
                'ftp_fna':Genomes,
                'ftp_faa':Proteins, 
                'ftp_gff':Gff, 
                'ftp_gene':Genes, 
                'ftp_report':Assembly_report,
                }

    for ftp_file, local_folder in ftp_location.items():

        file_name = local_folder / species[ftp_file].replace('.gz','')

        if not file_name.is_file():
            gbff_url = "{}/{}".format(species['ftp_path'], species[ftp_file]).replace('ftp:', 'https:')
            
            logging.debug(f"-> Using or downloading/using {gbff_url}")

            gbff_response = session.get(gbff_url)
            md5_gbff = hashlib.md5(gbff_response.content)

            sub_md5_gbff = md5_gbk[md5_gbk.assembly_files.str.contains(species[ftp_file])]

            if not sub_md5_gbff.empty and sub_md5_gbff.md5.values[0] == md5_gbff.hexdigest() :
                
                # print("\n-> md5 CHECKED OK")
                logging.debug(f'MD5 OK and checked for -> {gbff_url}')

                if species[ftp_file].endswith(".gz"):
                    write_file_uncompress(gbff_response, file_name)
                else:
                    write_file(gbff_response, file_name)

            else :
                logging.debug(f'Did not check, erasing the file -> {gbff_url}')
        else:
            logging.debug(f'This file already exists -> {file_name}')

    return 

##########################################################################################

def get_info_report(report_files):

    header_report = [
                "Sequence-Name", 
                "Sequence-Role",
                "Assigned-Molecule",
                "Assigned-Molecule-Location/Type",
                "GenBank-Accn",
                "Relationship",
                "RefSeq-Accn",
                "Assembly-Unit",
                "Sequence-Length",
                "UCSC-style-name",
                ]

    dict_genbank = {}
    dict_ncid = {}
    for report in report_files:
        name_GCA = '_'.join(report.stem.split('_')[:2])
        ids = pd.read_table(report, comment='#', names=header_report)[['GenBank-Accn', 'RefSeq-Accn']]
        for index, genbank, ncid in ids.itertuples():
            dict_genbank[name_GCA] = genbank.split('.')[0]
            dict_ncid[name_GCA] = ncid

    return dict_genbank, dict_ncid

##########################################################################################

def create_metadata(ictv_xlsx, assembly_table, dict_ncid, dict_genbank, output_metadata):

    ictv_file = pd.read_excel(ictv_xlsx)

    # Change the list of GENBANK accession to list
    ictv_file['Virus GENBANK accession'] = ictv_file['Virus GENBANK accession'].apply(lambda x: x.split('; ') if x == x else '')
    # Create one line per GENBANK accession ids
    ictv_file = ictv_file.explode('Virus GENBANK accession')
    # Take only the important part of the name
    ictv_file['Virus GENBANK accession'] = ictv_file['Virus GENBANK accession'].apply(lambda x: x.split(': ')[-1] if x == x else '')

    # Get the assembly table and add the new columns
    assembly_table['NC_Id'] = assembly_table.assembly_accession.map(dict_ncid)
    assembly_table['Virus GENBANK accession'] = assembly_table.assembly_accession.map(dict_genbank)

    name2keep_from_assembly = [
                    'assembly_accession',
                    'NC_Id',
                    'Virus GENBANK accession',
                    'bioproject',
                    'taxid',
                    'species_taxid',
                    'asm_name',
                    'gbrs_paired_asm'
                    ]

    ictv_file = ictv_file.merge(assembly_table[name2keep_from_assembly], on='Virus GENBANK accession', how='left')

    ictv_file.to_csv(output_metadata, index=False, sep='\t')

    return ictv_file

##########################################################################################

def rename_file(ictv_df, all_folders):

    ictv_df['new_name'] = ictv_df.apply(lambda x: f"{x['Species']}_{x['Virus GENBANK accession']}_{x['assembly_accession']}")
    GCA2name = ictv_df.set_index('assembly_accession').new_name.to_dict()

    for folder in all_folders:
        for myfile in folder.glob('*'):
            myfile_GCA = "_".join(myfile.stem.split("_")[:2])
            myfile_ext = myfile.suffix

            myfile_GCA = GCA2name(myfile_GCA)

            logging.info(f"Renaming {myfile.name} to {myfile_GCA + myfile_ext}")

            myfile.rename(Path(myfile.parent, myfile_GCA + myfile_ext))

    return 

##########################################################################################

def rename_gene_files(gene_files):

    for gene_file in gene_files:
        parser = SeqIO.parse(gene_file, 'fasta')

        tmp_file = str(gene_file) + ".tmp"
        gene_file_name = str(gene_file)

        with open(tmp_file, 'wt') as w_file:
            for gene in parser:
                gene.id = " gene_index:".join(gene.id.split("_")[-2:])
                gene.name = ""

                SeqIO.write(gene, w_file, "fasta")

        tmp_file = Path(gene_file.parent, tmp_file)
        gene_file.unlink()
        tmp_file.rename(gene_file_name)

    return
##########################################################################################

def concat_files():

    folder2concat = {
                Genomes:"ICTV_all_genomes.fna",
                Proteins:"ICTV_all_proteins.faa", 
                Genes:"ICTV_all_genes.fna", 
    }

    for folder, concat_file in folder2concat.items():
        with open(concat_file, "wt") as w_file:
            for myfile in folder.glob("*"):
                with open(myfile) as r_file:
                    for line in r_file:
                        w_file.write(line)

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
                            default=os.getcwd(),
                            dest="output",
                            metavar='<FOLDER>',
                            type=str,
                            help=f"Folder where to put the database (default: {os.getcwd()})")
general_option.add_argument("-d",'--date_stamp',
                            metavar="<DATE>",
                            dest="date_stamp",
                            help="The time stamp: YYYYMMDD (e.g. 20221016 for 16 october 2022)",
                            type=int,
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
general_option.add_argument(
    "-ictv",
    "--ictv_metadata",
    metavar="<ictv_metadata_table>",
    dest="ictv_metadata",
    help="Path to the ICTV taxonomy table that contain the Genbank Id",
    required=True,
    default="",
    type=str,
)

args = parser.parse_args()

##########################################################################################

taxa = Path(args.output) / "ICTV_database" / f"{args.date_stamp}"
Genomes = taxa / "Genomes"
GenBank = taxa / "GenBank"
Genes = taxa / "Genes"
Proteins = taxa / "Proteins"
Gff = taxa / "Gff"
Assembly_report = taxa / "Assembly_Reports"

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

queueListerner, mpQueue = logger_init()

logging.info(f'Gembase creation logging for version : ICTV_database_{args.date_stamp}')


##########################################################################################

list_viral = "viral_assembly_summary.txt"
viral_table_assembly = "https://ftp.ncbi.nlm.nih.gov/genomes/genbank/viral/assembly_summary.txt"
assembly_summary_viral_file = taxa / list_viral

##########################################################################################

logging.info(f"-> Fetching {viral_table_assembly}")
print(f"-> Fetching {viral_table_assembly}")

# Test if the file exist before to not have problem if the assembly file of NCBI change between two run for the same MICROBIAL
if assembly_summary_viral_file.is_file() :
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

pool = multiprocessing.Pool(processes=args.threads, initializer=init_process, initargs=[mpQueue])
results = list(
    tqdm(
        pool.imap(fetch_genbank_file, args_func), total=num_rows
    )
)
pool.close()
queueListerner.stop()

logging.info('Done!')
print('\nDone!\n')

##### TODO RENAME OF FILES + MERGE_REPORT
logging.info('-> Creating metadata for the database')
print('-> Creating metadata for the database')

Path()

dict_genbank, dict_ncid = get_info_report(Assembly_report.glob("*txt"))
ictv_df = create_metadata(
                    ictv_xlsx=args.ictv_metadata,
                    assembly_table=assembly_summary_viral,
                    dict_ncid=dict_ncid,
                    dict_genbank=dict_genbank,
                    output_metadata=taxa / "ICTV_metadata.tsv")

logging.info('Done!')
print('\nDone!\n')

logging.info('-> Rename files based on virus species name, genbank accession and assembly accession')
print('-> Rename files based on virus species name, genbank accession and assembly accession')

rename_file(ictv_df=ictv_df, all_folders=[Genomes, Genes, GenBank, Proteins, Gff, Assembly_report])

rename_gene_files(Genes)

logging.info('Done!')
print('\nDone!\n')

logging.info('-> Concatenation of Genomes, Genes and Proteins into one file each')
print('-> Concatenation of Genomes, Genes and Proteins into one file each')

concat_files()

logging.info('Done!')
print('\nDone!\n')

##########################################################################################
