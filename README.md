# ICTV Download and Creation

This script downloads and processes genome data from the International Committee on Taxonomy of Viruses (ICTV) database. It creates separate files for genomes, GenBank files, genes, proteins, Gff, and Lst.

## Installation
To use this script, you will need to have Python 3 installed on your system. You can download Python from the official website: https://www.python.org/downloads/

You will also need to install the following Python packages:

- pandas = 1.4.1
- tqdm = 4.62.3
- python-wget = 3.2
- biopython = 1.81
- bcbio-gff = 0.7.0

You can install these packages using pip:

```
pip install pandas==1.4.1 tqdm==4.62.3 wget==3.2 biopython==1.81 bcbio-gff==0.7.0
```

or mamba:

```
mamba install -c conda-forge -c bioconda pandas=1.4.1 tqdm=4.62.3 python-wget=3.2 biopython=1.81 bcbio-gff=0.7.0
```

## Usage

To run the script, use the following command:

```
python create_ictv.py -o <output_folder> -d <date_stamp> -v <verbosity_level> -t <num_threads> -ictv <ictv_metadata_table>
```

- <output_folder>: The folder where the database will be created (default: current working directory).
- <date_stamp>: The time stamp in the format YYYYMMDD (default: current date).
- <verbosity_level>: The level of logging verbosity (default: 1).
- <num_threads>: The number of threads to use (default: 1).
- <ictv_metadata_table>: The path to the ICTV taxonomy xlsx table that contains the Genbank Id (default: ICTV_metadata.xlsx).

## Output

The script will create a folder named ICTV_database_<date_stamp> in the specified output folder. This folder will contain the following subfolders:

- Genomes: Contains the genome sequences in FASTA format.
- GenBank: Contains the GenBank files.
- Genes: Contains the gene sequences in FASTA format.
- Proteins: Contains the protein sequences in FASTA format.
- Gff: Contains the Gff files.
- Lst: Contains the Lst files.

In addition, the script will create a file named ICTV_metadata.tsv in the ICTV_database_<date_stamp> folder. This file contains the metadata for the downloaded genomes.

# Contributing
If you would like to contribute to this project, please submit a pull request with your changes.

# License
This project is licensed under the MIT License - see the LICENSE file for details.
