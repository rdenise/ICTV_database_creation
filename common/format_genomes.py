#!/usr/bin/env python3
# -*- coding: utf-8 -*-

##########################################################################################
##########################################################################################
##
##                                Library
##
##########################################################################################
##########################################################################################

import sys, os
import pandas as pd
import logging
from tqdm import tqdm
import Bio
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from BCBio import GFF


##########################################################################################
##########################################################################################
##
##                                Functions
##
##########################################################################################
##########################################################################################


def gbk2file(accGenBank_nameFile_taxa):
    """
    Function that will get the assembly file from the ncbi server

    :params url: the url of the file (e.g. "https://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank.txt")
    :type: str
    :return: A dataframe with the information of the assembly file
    :rtype: pandas.dataframe
    """

    accGenBank_nameFile, taxa = accGenBank_nameFile_taxa

    Genomes = taxa / "Genomes"
    GenBank = taxa / "GenBank"
    Genes = taxa / "Genes"
    Proteins = taxa / "Proteins"
    Gff = taxa / "Gff"
    Lst = taxa / "Lst"

    nameFile = accGenBank_nameFile["File_identifier"]

    gbk_file = GenBank / f"{nameFile}.gbk"

    gbk = SeqIO.read(gbk_file, "genbank")

    # GFF
    gff_file = Gff / f"{nameFile}.gff"

    if not gff_file.is_file():
        logging.debug(f"-> Creating: {nameFile}.gff")

        gbk2gff3(gbk_file=gbk_file, outfile=gff_file)

    # global counter_gff
    # counter_gff.increment()

    # Genome fasta
    fasta_file = Genomes / f"{nameFile}.fna"

    if not fasta_file.is_file():
        logging.debug(f"-> Creating: {nameFile}.fna")

        gbk2fasta(replicon=gbk, fasta_file=fasta_file)

    # global counter_fna
    # counter_fna.increment()

    # Lst file
    lst_file = Lst / f"{nameFile}.lst"

    if not lst_file.is_file():
        # create a file to quickly have the informations
        logging.debug(f"-> Creating: {nameFile}.lst")

        lst_df = gbk2lst(replicon=gbk, lst_file=lst_file)

    logging.debug(f"-> Reading: {nameFile}.lst")

    dtype_lst = {
        "start": "int32",
        "end": "int32",
        "strand": "string",
        "type": "string",
        "gene": "string",
        "status": "string",
        "synonyms": "string",
        "locus_tag": "string",
        "nexons": "int32",
        "positions": "string",
        "sequence_nt": "string",
        "sequence_aa": "string",
        "product_note": "string",
    }

    lst_df = pd.read_table(lst_file, dtype=dtype_lst)
    lst_df.fillna("", inplace=True)

    # global counter_lst
    # counter_lst.increment()

    # Gene and protein fasta
    gen_file = Genes / f"{nameFile}.genes.fna"
    prt_file = Proteins / f"{nameFile}.faa"

    if not prt_file.is_file():
        # Gene fasta
        logging.debug(f"-> Creating: {nameFile}.gene.fna")

        valid_df = gbk2gen(df_lst=lst_df, gen_file=gen_file)

        # global counter_gene_fna
        # counter_gene_fna.increment()

        # Protein fasta
        logging.debug(f"-> Creating: {nameFile}.faa")

        gbk2prt(prt_file=prt_file, df_lst_Valid_CDS=valid_df)

        # global counter_faa
        # counter_faa.increment()

    return


##########################################################################################


def gbk2fasta(replicon, fasta_file):

    """
    Function that will create the .fst file for a specific replicon
    Converting the genbank to fasta for the full sequence of the chromosome/plasmid

    :params replicon: The actual replicon with all the informations of the genbank
    :type: Bio.SeqRecord.SeqRecord
    :params replicon_id: The gembase id of the replicon
    :type: str
    :params fasta_file: Path of the new created file
    :type: str
    """

    date = replicon.annotations["date"]

    len_bp = "{bp_size}bp".format(bp_size=len(replicon))

    replicon.id = replicon.name

    replicon.name = ""

    replicon.description = f"{len_bp} {replicon.description} [{date}]"

    SeqIO.write(replicon, fasta_file, "fasta")

    return


##########################################################################################


def gbk2lst(replicon, lst_file):

    """
    Function that will create the .lst file for a specific replicon

    :params replicon: The actual replicon with all the informations of the genbank
    :type: Bio.SeqRecord.SeqRecord
    :params replicon_id: The gembase id of the replicon
    :type: str
    :params lst_file: Path of the new created file
    :type: str
    :return: The dataframe that was write in .lst file
    :rtype: pandas.DataFrame
    """

    tmp_dict = {
        "start": [],
        "end": [],
        "strand": [],
        "type": [],
        "gene": [],
        "status": [],
        "synonyms": [],
        "locus_tag": [],
        "nexons": [],
        "positions": [],
        "sequence_nt": [],
        "sequence_aa": [],
        "product_note": [],
    }

    # Allow to only have the interesting features of the gbk without redondancy
    new_list = [x for x in replicon.features if x.type not in ("source", "gene")]

    for sequence in new_list:

        tmp_dict["status"].append("")

        # As the gene could be at multiple position because of programmed frameshift or at the extremity of the chromosome
        list_position = []

        for position in sequence.location.parts:
            if Bio.SeqFeature.BeforePosition == type(
                position.start
            ) or Bio.SeqFeature.AfterPosition == type(position.end):
                tmp_dict["status"][-1] = "Partial"

            # Add position to the list of position and start +1 because python begin at 0
            list_position.extend([position.nofuzzy_start + 1, position.nofuzzy_end])

        # I do not know why but even if biopython know where is the start and stop, he order weirdly when it is in the complementary strand
        list_position = sorted(list_position)

        # Because bug of biopython, when gene on the virtual cut of the chromosome it does not know how to identify the start and the stop
        if sequence.location_operator == "join" and sequence.location.start == 0:
            # put the first 2 element at the end
            for j in range(2):
                list_position.append(list_position.pop(0))

        list_position = [str(position) for position in list_position]

        tmp_dict["positions"].append(" ".join(list_position))

        tmp_dict["start"].append(int(list_position[0]))
        tmp_dict["end"].append(int(list_position[-1]))

        # The only states that I found in the .lst of MICROBIAL_D Partial, Pseudo, Invalid_Size, Valid :
        if tmp_dict["status"][-1] == "":
            if "pseudo" in sequence.qualifiers or "pseudogene" in sequence.qualifiers:
                tmp_dict["status"][-1] = "Pseudo"
            elif (
                "note" in sequence.qualifiers
                and "frameshift" in " ".join(sequence.qualifiers["note"]).lower()
            ):
                tmp_dict["status"][-1] = "Invalid_Size"
            else:
                tmp_dict["status"][-1] = "Valid"

        # D if in direct strand (=> strand = 1), C in complementary strand (=> strand = -1)
        tmp_dict["strand"].append("+" if sequence.strand == 1 else "-")
        tmp_dict["type"].append(sequence.type)
        tmp_dict["nexons"].append(len(sequence.location.parts))

        # Because in the gembase file the gene name and the locus_tag are the same so the information is twice in the file
        # But I have no clue why. Maybe because the file from assembly have a GenBank file without the gene part
        if "gene" in sequence.qualifiers:
            tmp_dict["gene"].append(" ".join(sequence.qualifiers["gene"]))
            if "locus_tag" in sequence.qualifiers:
                tmp_dict["locus_tag"].append(" ".join(sequence.qualifiers["locus_tag"]))
            else:
                tmp_dict["locus_tag"].append(" ".join(sequence.qualifiers["gene"]))
        elif "locus_tag" in sequence.qualifiers:
            tmp_dict["gene"].append(" ".join(sequence.qualifiers["locus_tag"]))
            tmp_dict["locus_tag"].append(" ".join(sequence.qualifiers["locus_tag"]))
        else:
            tmp_dict["gene"].append("NoLocusTag")
            tmp_dict["locus_tag"].append("NoLocusTag")

        if "protein_id" in sequence.qualifiers:
            tmp_dict["synonyms"].append(" ".join(sequence.qualifiers["protein_id"]))
        else:
            tmp_dict["synonyms"].append("")

        if "product" in sequence.qualifiers:
            product = " ".join(sequence.qualifiers["product"])
        else:
            product = ""

        if "note" in sequence.qualifiers:
            note = " ".join(sequence.qualifiers["note"])
        else:
            note = ""

        tmp_dict["product_note"].append(f"| {product} | {note}")

        if tmp_dict["type"][-1] == "CDS" and tmp_dict["status"][-1] == "Valid":
            tmp_dict["sequence_aa"].append(sequence.qualifiers["translation"][0])
            tmp_dict["sequence_nt"].append(str(sequence.extract(replicon).seq))
        else:
            tmp_dict["sequence_aa"].append("")
            tmp_dict["sequence_nt"].append("")

    df = pd.DataFrame(tmp_dict)

    df.to_csv(lst_file, index=False, sep="\t")

    return df


##########################################################################################


def gbk2gen(df_lst, gen_file):

    """
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
    """

    df_lst_Valid_CDS = df_lst[
        (df_lst.type == "CDS") & (df_lst.status == "Valid")
    ].reset_index(drop=True)
    dict_lst_Valid_CDS = df_lst_Valid_CDS.to_dict("records")
    # List of all the sequences to write
    list_sequences = []

    all_descriptions = []
    all_starts = []
    all_stops = []

    for sequence in dict_lst_Valid_CDS:

        start_codon = str(sequence["sequence_nt"][:3])
        stop_codon = str(sequence["sequence_nt"][-3:])

        product_note = sequence["product_note"]

        gene_seq = SeqRecord(Seq(sequence["sequence_nt"]))

        size_gene = len(gene_seq)

        gene_seq.id = sequence["synonyms"]
        gene_seq.name = ""

        list_description = [
            sequence["strand"],
            start_codon,
            stop_codon,
            str(sequence["start"]),
            str(sequence["end"]),
            sequence["gene"],
            str(size_gene),
            sequence["synonyms"],
            sequence["locus_tag"],
            str(sequence["nexons"]),
            sequence["positions"],
            product_note,
        ]

        list_description = [x if x == x else "" for x in list_description]

        gene_seq.description = " ".join(list_description)

        # We add the gene modify in description to the list
        list_sequences.append(gene_seq)

        all_descriptions.append(" ".join(list_description))
        all_starts.append(start_codon)
        all_stops.append(stop_codon)

    df_lst_Valid_CDS["description_gembase"] = all_descriptions
    df_lst_Valid_CDS["start_codon"] = all_starts
    df_lst_Valid_CDS["stop_codon"] = all_stops

    SeqIO.write(list_sequences, gen_file, "fasta")

    return df_lst_Valid_CDS


##########################################################################################


def gbk2gff3(gbk_file, outfile):

    with open(outfile, "wt") as w_file:
        GFF.write(SeqIO.parse(gbk_file, "genbank"), w_file)

    return


##########################################################################################


def gbk2prt(prt_file, df_lst_Valid_CDS):

    """
    Function that will create the .prt file for a specific replicon
    Converting the genbank to multi fasta for each proteins in the replicon

    :params replicon: The actual replicon with all the informations of the genbank
    :type: Bio.SeqRecord.SeqRecord
    :params prt_file: Path of the new created file
    :type: str
    :params df_lst_Valid_CDS: The sub dataframe of the genes CDS and Valid to have the same description
             of .gen file without calculating it again
    :type: pandas.DataFrame
    """

    # Transform the features of the gbk into a dict with the name of the locus_tag as key and feature as value
    # Allow to keep the order of the .lst file but using the translation of the .gbk without calculate it
    dict_lst_Valid_CDS = df_lst_Valid_CDS.to_dict("records")

    # List of all the proteins to write
    all_proteins = []

    for sequence in dict_lst_Valid_CDS:

        protein = SeqRecord(Seq(sequence["sequence_aa"]))

        protein.id = sequence["synonyms"]

        # There is two spaces before the parenthesis in the .prt format so I keep it
        protein.description = f"{sequence['description_gembase']} (translation)"

        all_proteins.append(protein)

    SeqIO.write(all_proteins, prt_file, "fasta")

    return


##########################################################################################


def concat_files(Genomes, Proteins, Genes, taxa):

    folder2concat = {
        Genomes: taxa / "ICTV_all_genomes.fna",
        Proteins: taxa / "ICTV_all_proteins.faa",
        Genes: taxa / "ICTV_all_genes.fna",
    }

    for folder, concat_file in folder2concat.items():
        with open(concat_file, "wt") as w_file:
            files2concat = folder.glob("*")
            num_file = len(list(files2concat))
            for myfile in tqdm(files2concat, colour="GREEN", total=num_file):
                with open(myfile) as r_file:
                    num_line = r_file.read().count("\n")
                    for line in tqdm(
                        r_file, colour="BLUE", leave=False, total=num_line
                    ):
                        w_file.write(line)

    return
