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

from common.download_genomes import *
from common.utils import *
from common.format_genomes import *

##########################################################################################

currentDate = datetime.date.today()

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
    choices=range(1, multiprocessing.cpu_count() + 1),
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

queueListerner, mpQueue = logger_init(level=level, args=args)

logging.info(f"ICTV creation logging for version : ICTV_database_{args.date_stamp}\n")


##########################################################################################

# File for the metadata of ICTV could be found here: https://talk.ictvonline.org/taxonomy/vmr/m/vmr-file-repository/

logging.info("\n-> Reading all the information on ICTV and exploding the identifiers\n")
print("-> Reading all the information on ICTV and exploding the identifiers")


ictv_df = pd.read_excel(args.ictv_metadata)
ictv_df = ictv_df[~ictv_df["Virus GENBANK accession"].isna()].reset_index(drop=True)

# Change the list of GENBANK accession to list remove space
ictv_df["Virus GENBANK accession"] = ictv_df["Virus GENBANK accession"].apply(
    lambda x: x.replace("; ", ";") if x == x else ""
)

# Change the list of GENBANK accession to list remove space
ictv_df["Virus GENBANK accession"] = ictv_df["Virus GENBANK accession"].apply(
    lambda x: x.replace(": ", ":") if x == x else ""
)

# Change the list of GENBANK accession to list
ictv_df["Virus GENBANK accession"] = ictv_df["Virus GENBANK accession"].apply(
    lambda x: x.split(";") if x == x else ""
)

# Create one line per GENBANK accession ids
ictv_df = ictv_df.explode("Virus GENBANK accession")

# Some have a " and "
ictv_df["Virus GENBANK accession"] = ictv_df["Virus GENBANK accession"].apply(
    lambda x: x.split(" and ") if x == x else ""
)

# Create one line per GENBANK accession ids
ictv_df = ictv_df.explode("Virus GENBANK accession")

# Some have a ", "
ictv_df["Virus GENBANK accession"] = ictv_df["Virus GENBANK accession"].apply(
    lambda x: x.split(",") if x == x else ""
)

# Create one line per GENBANK accession ids
ictv_df = ictv_df.explode("Virus GENBANK accession")

# Take only the important part of the name
ictv_df["Virus GENBANK accession"] = ictv_df["Virus GENBANK accession"].apply(
    lambda x: x.split(":")[-1] if x == x else ""
)

# Make sure identifier don't have space
ictv_df["Virus GENBANK accession"] = ictv_df["Virus GENBANK accession"].apply(
    lambda x: x.strip()
)

# Some have a " "
ictv_df["Virus GENBANK accession"] = ictv_df["Virus GENBANK accession"].apply(
    lambda x: x.split(" ") if x == x else ""
)

# Create one line per GENBANK accession ids
ictv_df = ictv_df.explode("Virus GENBANK accession")

# Remove empty
ictv_df = ictv_df[~(ictv_df["Virus GENBANK accession"] == "")].reset_index(drop=True)

# Changing the name to have a good one Species.Notes.GenBankAcc
ictv_df["File_identifier"] = ictv_df.apply(
    lambda x: f"{x.Species.replace(' ', '_')}.{x.Sort}.{x['Virus GENBANK accession']}",
    axis=1,
)

ictv_df.to_csv(taxa / "ICTV_metadata.tsv", index=False, sep="\t")

logging.info("-> Creating all the files for each genomes in ICTV\n")
print("-> Creating all the files for each genomes in ICTV")

num_rows = ictv_df.shape[0]

# RED, GREEN, YELLOW, BLUE, MAGENTA, CYAN, WHITE
# counter_gbk = Counter(desc="Gbk processed", colour="GREEN", total=num_rows)
# counter_gff = Counter(desc="Gff processed", colour="CYAN", total=num_rows)
# counter_fna = Counter(desc="Lst processed", colour="BLUE", total=num_rows)
# counter_lst = Counter(desc="Fna processed", colour="MAGENTA", total=num_rows)
# counter_gene_fna = Counter(desc="Gen processed", colour="RED", total=num_rows)
# counter_faa = Counter(desc="Prt processed", colour="YELLOW", total=num_rows)

##### MULTIPROCESS ACTION
args_func = ictv_df[["Virus GENBANK accession", "File_identifier"]].to_dict("records")

args_func = [[i, taxa] for i in args_func]

if args.threads > 1:
    subthreads = 9 if args.threads > 9 else args.threads

    pool = multiprocessing.Pool(
        processes=subthreads, initializer=init_process, initargs=[mpQueue, level]
    )
    results = list(
        tqdm(
            pool.imap(efetch_accession2gbk, args_func),
            desc="Genomes download",
            total=num_rows,
            colour="GREEN",
        )
    )
    pool.close()

    pool = multiprocessing.Pool(
        processes=args.threads, initializer=init_process, initargs=[mpQueue, level]
    )
    results = list(
        tqdm(
            pool.imap(gbk2file, args_func),
            desc="Genomes processed",
            total=num_rows,
            colour="BLUE",
        )
    )
    pool.close()

    queueListerner.stop()
else:
    for pair_acc in tqdm(
        args_func, desc="Genomes download", total=num_rows, colour="GREEN"
    ):
        efetch_accession2gbk(pair_acc)

    for pair_acc in tqdm(
        args_func, desc="Genomes processed", total=num_rows, colour="BLUE"
    ):
        gbk2file(pair_acc)


# counter_gbk.close()
# counter_gff.close()
# counter_fna.close()
# counter_lst.close()
# counter_gene_fna.close()
# counter_faa.close()

logging.info("Done!")
print("\nDone!\n")

logging.info("-> Concatenation of Genomes, Genes and Proteins into one file each\n")
print("-> Concatenation of Genomes, Genes and Proteins into one file each")

concat_files(Genomes=Genomes, Proteins=Proteins, Genes=Genes, taxa=taxa)

logging.info("Done!")
print("\nDone!\n")

##########################################################################################
