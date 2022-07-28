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
import logging
from logging.handlers import QueueHandler, QueueListener
from tqdm import tqdm
import multiprocessing
from pathlib import Path
import datetime
from Bio import Entrez, SeqIO
from BCBio import GFF

##########################################################################################
##########################################################################################
##
##                                Classes
##
##########################################################################################
##########################################################################################


class Counter(object):
    def __init__(self, initval=0):
        self.val = multiprocessing.RawValue('i', initval)
        self.lock = multiprocessing.Lock()

    def increment(self):
        with self.lock:
            self.val.value += 1

    @property
    def value(self):
        return self.val.value

##########################################################################################

currentDate = datetime.date.today()
Entrez.email = "rdenise@ucc.ie"
Entrez.tool = "create_ictv.py"

counter_gbk = Counter()
counter_faa = Counter()
counter_fna = Counter()
counter_gene_fna = Counter()
counter_gff = Counter()
counter_lst = Counter()

##########################################################################################
##########################################################################################
##
##                                Functions
##
##########################################################################################
##########################################################################################


def logger_init(level):
    mpQueue = multiprocessing.Queue()

    LOG_FORMAT = "%(levelname)s::%(asctime)s - %(message)s"
    logging.basicConfig(
        filename=os.path.join(args.output, "ictv_downloading.log"),
        level=level,
        format=LOG_FORMAT,
        filemode="w",
    )

    # this is the handler for all log records
    handler = logging.FileHandler(filename=os.path.join(args.output, "ictv_downloading.log"), mode="w")
    LOG_FORMAT_HANDLER = "%(levelname)s: %(asctime)s - %(process)s - %(message)s"
    handler.setFormatter(logging.Formatter(LOG_FORMAT_HANDLER))

    # queueListerner gets records from the queue and sends them to the handler
    queueListerner = QueueListener(mpQueue, handler)
    queueListerner.start()

    logger = logging.getLogger()
    logger.setLevel(level)
    # add the handler to the logger so records from this process are handled
    logger.addHandler(handler)
    logger.addHandler(handler)

    return queueListerner, mpQueue


##########################################################################################


def init_process(mpQueue, level):
    
    # all records from worker processes go to queueHandler and then into mpQueue
    queueHandler = QueueHandler(mpQueue)
    logger = logging.getLogger()
    logger.setLevel(level)
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


def efetch_accession2gbk(accGenBank, nameFile):
    """
    Function that will get the assembly file from the ncbi server

    :params url: the url of the file (e.g. "https://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank.txt")
    :type: str
    :return: A dataframe with the information of the assembly file
    :rtype: pandas.dataframe
    """

    global counter_gbk
    global pbar_gbk

    gbk_file = GenBank / f"{nameFile}.gbk"
    
    if gbk_file.is_file():
        logging.debug(f"-> Reading: {nameFile}.gbk")
    
        gbk = SeqIO.read(gbk_file, "genbank")
    else:
        logging.debug(f"-> Downloading the identifier {accGenBank}")

        handle = Entrez.efetch(db="nucleotide", rettype="gbwithparts", retmode="text",
                            id=accGenBank)

        gbk = SeqIO.read(handle, 'genbank')

        logging.debug(f"-> Creating: {nameFile}.gbk")


        SeqIO.write(gbk, gbk_file, format='genbank')

    gff_file = Gff / f"{nameFile}.gff"
    
    if not gff_file.is_file():
        logging.debug(f"-> Creating: {nameFile}.gff")

        gbk2gff3(gbk_file=gbk_file, outfile=gff_file)

    fasta_file = Genomes / f"{nameFile}.fna"

    if not fasta_file.is_file():
        logging.debug(f"-> Creating: {nameFile}.fna")

        gbk2fasta(replicon=gbk, fasta_file=fasta_file)

    lst_file = Lst / f"{nameFile}.lst"

    if lst_file.is_file():
        logging.debug(f"-> Reading: {nameFile}.gbk")

        lst_df = pd.read_table(lst_file)
    else:
        # create a file to quickly have the informations 
        logging.debug(f"-> Creating: {nameFile}.lst")

        lst_df = gbk2lst(gbk=gbk, lst_file=lst_file)

    gen_file = Genes / f"{nameFile}.genes.fna"
    
    if not prt_file.is_file(): 
        logging.debug(f"-> Creating: {nameFile}.gene.fna")

        valid_df = gbk2gen(df_lst=lst_df, gen_file=gen_file)
        
        logging.debug(f"-> Creating: {nameFile}.faa")

        prt_file = Proteins / f"{nameFile}.faa"
        gbk2prt(prt_file=prt_file, df_lst_Valid_CDS=valid_df)

    counter_gbk.increment()
    pbar_gbk.update(counter_gbk.value())

    return 


##########################################################################################


def concat_files():

    folder2concat = {
        Genomes: taxa / "ICTV_all_genomes.fna",
        Proteins: taxa / "ICTV_all_proteins.faa",
        Genes: taxa / "ICTV_all_genes.fna",
    }

    for folder, concat_file in folder2concat.items():
        if not concat_file.is_file():
            with open(concat_file, "wt") as w_file:
                for myfile in tqdm(folder.glob("*"), colour="GREEN"):
                    with open(myfile) as r_file:
                        for line in tqdm(r_file, colour="BLUE", leave=False):
                            w_file.write(line)

    return

##########################################################################################

def gbk2fasta(replicon, fasta_file) :

    '''
    Function that will create the .fst file for a specific replicon
    Converting the genbank to fasta for the full sequence of the chromosome/plasmid

    :params replicon: The actual replicon with all the informations of the genbank 
    :type: Bio.SeqRecord.SeqRecord
    :params replicon_id: The gembase id of the replicon
    :type: str
    :params fasta_file: Path of the new created file
    :type: str    
    '''

    global counter_fna
    global pbar_fna

    date = replicon.annotations['date'] 

    len_bp = '{bp_size}bp'.format(bp_size=len(replicon))

    replicon.id = replicon.name

    replicon.name = ''

    replicon.description = f'{len_bp} {replicon.description} [{date}]'

    SeqIO.write(replicon, fasta_file, 'fasta')

    counter_fna.increment()
    pbar_fna.update(counter_fna.value())

    return
    
##########################################################################################

def gbk2lst(replicon, lst_file) :

    '''
    Function that will create the .lst file for a specific replicon

    :params replicon: The actual replicon with all the informations of the genbank 
    :type: Bio.SeqRecord.SeqRecord
    :params replicon_id: The gembase id of the replicon
    :type: str
    :params lst_file: Path of the new created file
    :type: str
    :return: The dataframe that was write in .lst file
    :rtype: pandas.DataFrame
    '''

    global counter_lst
    global pbar_lst

    tmp_dict = {
        'start':[],
        'end':[],
        'strand':[],
        'type':[],
        'gembase_Id':[],
        'gene':[],
        'status':[],
        'synonyms':[],
        'locus_tag':[],
        'nexons':[],
        'position':[],
        'sequence_nt':[],
        'sequence_aa':[],
        'product_note':[]
        } 

    # Allow to only have the interesting features of the gbk without redondancy
    new_list = [x for x in replicon.features if x.type not in ('source', 'gene')]    

    for sequence in new_list : 

        tmp_dict['status'].append('')

        # As the gene could be at multiple position because of programmed frameshift or at the extremity of the chromosome
        list_position = []

        for position in sequence.location.parts :
            if Bio.SeqFeature.BeforePosition == type(position.start) or Bio.SeqFeature.AfterPosition == type(position.end) :
                tmp_dict['status'][-1] = 'Partial'
            
            # Add position to the list of position and start +1 because python begin at 0
            list_position.extend([position.nofuzzy_start + 1, position.nofuzzy_end])
            
        # I do not know why but even if biopython know where is the start and stop, he order weirdly when it is in the complementary strand
        list_position = sorted(list_position)

        # Because bug of biopython, when gene on the virtual cut of the chromosome it does not know how to identify the start and the stop
        if sequence.location_operator == 'join' and sequence.location.start == 0 :
            # put the first 2 element at the end
            for j in range(2) :
                list_position.append(list_position.pop(0))
        
        list_position = [str(position) for position in list_position]

        tmp_dict['position'] = ' '.join(list_position)

        tmp_dict['start'] = int(list_position[0])
        tmp_dict['end'] = int(list_position[-1])

        # The only states that I found in the .lst of MICROBIAL_D Partial, Pseudo, Invalid_Size, Valid : 
        if tmp_dict['status'] == '' :
            if 'pseudo' in sequence.qualifiers :
                tmp_dict['status'] = 'Pseudo'
            elif 'frameshift' in ' '.join(sequence.qualifiers['note']).lower() : 
                tmp_dict['status'] = 'Invalid_Size'
            else :
                tmp_dict['status'] = 'Valid'

        # D if in direct strand (=> strand = 1), C in complementary strand (=> strand = -1)
        tmp_dict['strand'] = '+' if sequence.strand == 1 else '-' 
        tmp_dict['type'] = sequence.type
        tmp_dict['nexons'] = len(sequence.location.parts)

        # Because in the gembase file the gene name and the locus_tag are the same so the information is twice in the file
        # But I have no clue why. Maybe because the file from assembly have a GenBank file without the gene part
        if 'gene' in sequence.qualifiers :
            tmp_dict['gene'] = ' '.join(sequence.qualifiers['gene'])
            tmp_dict['locus_tag'] = ' '.join(sequence.qualifiers['locus_tag'])
        else :
            tmp_dict['gene'] = tmp_dict['locus_tag'] = ' '.join(sequence.qualifiers['locus_tag'])

        if 'protein_id' in sequence.qualifiers :
            tmp_dict['synonyms'] = ' '.join(sequence.qualifiers['protein_id'])
        else :
            tmp_dict['synonyms'] = ''

        product = ' '.join(sequence.qualifiers['product'])
        note = ' '.join(sequence.qualifiers['note'])
        tmp_dict['product_note'] = f'|{product}|{note}'

        tmp_dict['sequence_aa'] = sequence.qualifiers['translation']
        tmp_dict['sequence_nt'] = sequence.extract(replicon)


    df = pd.DataFrame(tmp_dict)

    df.to_csv(lst_file, index=False, sep='\t')

    counter_lst.increment()
    pbar_lst.update(counter_lst.value())

    return df


##########################################################################################

def gbk2gen(df_lst, gen_file) :

    '''
    Function that will create the .gen file for a specific replicon
    Converting the genbank to multifasta for each genes in the replicon

    :params replicon: The actual replicon with all the informations of the genbank 
    :type: Bio.SeqRecord.SeqRecord
    :params df_lst: The dataframe that was write in .lst file
    :type: pandas.DataFrame  
    :params gen_file: Path of the new created file
    :type: str   
    :return: The sub dataframe of the genes CDS and Valid to have the same description without
             calculating it again
    :rtype: pandas.DataFrame
    '''

    global counter_gene_fna
    global pbar_gene_fna

    df_lst_Valid_CDS = df_lst[(df_lst.type == 'CDS') & (df_lst.status == 'Valid')].reset_index(drop=True)
    dict_lst_Valid_CDS = df_lst_Valid_CDS.to_dict['records']
    # List of all the sequences to write
    sequences = []
    
    for index, sequence in dict_lst_Valid_CDS.items():

        start_codon = str(sequence['sequence_nt'][:3])
        stop_codon = str(sequence['sequence_nt'][-3:])      

        product_note = sequence['product_note']

        gene_seq = sequence['sequence_nt']

        size_gene = len(gene_seq)

        gene_seq.id = sequence['synonyms']
        gene_seq.name = ''

        list_description = [sequence.strand, start_codon, stop_codon, str(sequence.start), str(sequence.end),
                       sequence.gene, size_gene, sequence.synonyms, sequence.locus_tag,
                       sequence.nexons, sequence.positions, product_note]

        gene_seq.description = ' '.join(list_description)

        # We add the gene modify in description to the list
        sequences.append(gene_seq)

        df_lst_Valid_CDS['description_gembase'] = ' '.join(list_description)
        df_lst_Valid_CDS['start_codon'] = start_codon
        df_lst_Valid_CDS['stop_codon'] = stop_codon

    SeqIO.write(sequences, gen_file, 'fasta')

    counter_gene_fna.increment()
    pbar_gene_fna.update(counter_gene_fna.value())

    return df_lst_Valid_CDS

##########################################################################################

def gbk2gff3(gbk_file, outfile):

    global counter_gff
    global pbar_gff

    with open(outfile, "wt") as w_file:
        GFF.write(SeqIO.parse(gbk_file, "genbank"), w_file)

    counter_gff.increment()
    pbar_gff.update(counter_gff.value())

    return

##########################################################################################

def gbk2prt(prt_file, df_lst_Valid_CDS) :

    '''
    Function that will create the .prt file for a specific replicon
    Converting the genbank to multi fasta for each proteins in the replicon

    :params replicon: The actual replicon with all the informations of the genbank 
    :type: Bio.SeqRecord.SeqRecord
    :params prt_file: Path of the new created file
    :type: str      
    :params df_lst_Valid_CDS: The sub dataframe of the genes CDS and Valid to have the same description 
             of .gen file without calculating it again
    :type: pandas.DataFrame   
    '''

    global counter_faa
    global pbar_faa

    # Transform the features of the gbk into a dict with the name of the locus_tag as key and feature as value
    # Allow to keep the order of the .lst file but using the translation of the .gbk without calculate it
    dict_lst_Valid_CDS = df_lst_Valid_CDS.to_dict['records']

    # List of all the proteins to write
    proteins = []

    for index, sequence in dict_lst_Valid_CDS.items():
        seq = Seq(sequence['sequence_aa'])

        seq.id = sequence['synonyms']

        # There is two spaces before the parenthesis in the .prt format so I keep it
        seq.description = '{} (translation)'.format(sequence['description_gembase'])

        proteins.append(seq)

    SeqIO.write(proteins, prt_file, 'fasta')

    counter_faa.increment()
    pbar_faa.update(counter_faa.value())

    return

##########################################################################################
##########################################################################################
##
##                                Main
##
##########################################################################################
##########################################################################################

parser = argparse.ArgumentParser(
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description=dedent("""ICTV Download and Creation"""),
)

general_option = parser.add_argument_group(title="General input dataset options")
general_option.add_argument(
    "-o",
    "--output",
    default=os.getcwd(),
    dest="output",
    metavar="<FOLDER>",
    type=str,
    help=f"Folder where to put the database (default: {os.getcwd()})",
)
general_option.add_argument(
    "-d",
    "--date_stamp",
    metavar="<DATE>",
    dest="date_stamp",
    help=f"The time stamp: YYYYMMDD (e.g. {currentDate.strftime('%Y%m%d')} for {currentDate.strftime('%d %B %Y')})",
    type=int,
    default=currentDate.strftime("%Y%m%d"),
)
parser.add_argument(
    "-v",
    "--verbosity",
    default=1,
    action="count",
    dest="verbosity",
    help="Increase log verbosity could be : -v or -vv (default: -v)",
)
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
    help="Path to the ICTV taxonomy xlsx table that contain the Genbank Id, can be found: https://talk.ictvonline.org/taxonomy/vmr/m/vmr-file-repository/",
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
Lst = taxa / "Lst"

create_folder(taxa)
create_folder(Genomes)
create_folder(GenBank)
create_folder(Genes)
create_folder(Proteins)
create_folder(Gff)
create_folder(Lst)

##########################################################################################

if args.verbosity == 1:
    level = logging.INFO
elif args.verbosity == 2:
    level = logging.DEBUG
else:
    level = logging.NOTSET

queueListerner, mpQueue = logger_init(level)

logging.info(f"ICTV creation logging for version : ICTV_database_{args.date_stamp}\n")


##########################################################################################

# File for the metadata of ICTV could be found here: https://talk.ictvonline.org/taxonomy/vmr/m/vmr-file-repository/

logging.info("\n-> Reading all the information on ICTV and exploding the identifiers\n")
print("-> Reading all the information on ICTV and exploding the identifiers")


ictv_df = pd.read_excel(args.ictv_metadata)
ictv_df = ictv_df[~ictv_df["Virus GENBANK accession"].isna()].reset_index(drop=True)

tqdm.pandas(desc="Step1", colour="GREEN")

# Change the list of GENBANK accession to list
ictv_df["Virus GENBANK accession"] = ictv_df["Virus GENBANK accession"].progress_apply(
    lambda x: x.split("; ") if x == x else ""
)
# Create one line per GENBANK accession ids
ictv_df = ictv_df.explode("Virus GENBANK accession")

tqdm.pandas(desc="Step2", colour="GREEN")

# Take only the important part of the name
ictv_df["Virus GENBANK accession"] = ictv_df["Virus GENBANK accession"].progress_apply(
    lambda x: x.split(": ")[-1] if x == x else ""
)

tqdm.pandas(desc="Creating names", colour="GREEN")

# Changing the name to have a good one Species.Notes.GenBankAcc
ictv_df["File_identifier"] = ictv_df.progress_apply(
    lambda x: f"{x.Species}.{x.Sort}.{x['Virus GENBANK accession']}",
    axis = 1
)

ictv_df.to_csv(taxa / "ICTV_metadata.tsv", index=False, sep="\t")

##########################################################################################

num_rows = ictv_df.shape[0]

# RED, GREEN, YELLOW, BLUE, MAGENTA, CYAN, WHITE
pbar_gbk =tqdm(desc="Completely done", colour="GREEN", position=0, total=num_rows) 
pbar_gff =tqdm(desc="Gff processed", colour="CYAN", position=1, total=num_rows) 
pbar_fna =tqdm(desc="Lst processed", colour="BLUE", position=2, total=num_rows) 
pbar_lst =tqdm(desc="Fna processed", colour="MAGENTA", position=3, total=num_rows) 
pbar_gene_fna =tqdm(desc="Gen processed", colour="RED", position=4, total=num_rows) 
pbar_faa =tqdm(desc="Prt processed", colour="YELLOW", position=5, total=num_rows) 

##########################################################################################

logging.info("\n-> Creating all the files for each genomes in ICTV\n")
print("-> Creating all the files for each genomes in ICTV")

##### MULTIPROCESS ACTION
args_func = ictv_df[["Virus GENBANK accession", "File_identifier"]].to_dict("records")


pool = multiprocessing.Pool(
    processes=args.threads, initializer=init_process, initargs=[mpQueue, level]
)
results = list(pool.starmap(efetch_accession2gbk, args_func))
pool.close()
queueListerner.stop()

logging.info("Done!")
print("\nDone!\n")

logging.info("\n-> Concatenation of Genomes, Genes and Proteins into one file each\n")
print("-> Concatenation of Genomes, Genes and Proteins into one file each")

concat_files()

logging.info("Done!")
print("\nDone!\n")

##########################################################################################
