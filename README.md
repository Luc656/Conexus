# Conexus
![LogoConexus](https://user-images.githubusercontent.com/94873030/189007177-a9bbac60-7114-4790-a8cf-1dd5f08e2408.jpg)

Conexus (latin for connection) is a bioinformatics pipeline designed to find key co-occuring gene pairs with a link to AMR. The pipline includes a number of different stages using external software tools to highlight key gene pairings likely under the influence of natural selection in a microbial commuity from raw FASTA input. FASTA files are first annotated then scanned for known AMR genes before being clusted together forming an overall pangeneome where gene pairing frequencies can be measured returning those deemed statically significant. The pipeline also stores information on all discovered AMR genes in the population and their corresponding origin files allowing for in-depth future analysis
# Dependancies 
 - [Prokka](https://github.com/tseemann/prokka)
 - [RGI](https://github.com/arpcard/rgi)
 - [Panaroo](https://gtonkinhill.github.io/panaroo/#/gettingstarted/quickstart)
 - [Coinfinder](http://mcinerneylab.com/software/coinfinder/)
 - [Pandas](https://pandas.pydata.org/)
 - [os](https://docs.python.org/3/library/os.html#module-os)
 - [Subprocess](https://docs.python.org/3/library/subprocess.html#module-subprocess)
 - [CSV](https://docs.python.org/3/library/csv.html#module-csv)
# Requirements 
Conexus requires that input files are in FASTA format and end in '.fasta' while each file must also contain a unique ID wether a number or a code seprated by a full stop. The program requires that the position of the ID when the file name is split by full stops is suppied in the -id parameter. 

e.g.
```
This.is.file.number.001.fasta
```
Here the unique ID is at position 4 (counting from 0) when split by the full stops
# Installation & Usage
```
git clone https://github.com/Luc656/Conexus.git
```
Run the Conexus python script in the same directory that contains your sequences, the program will automatically begin processing any file ending in '.fasta'
# Process
##### Annotation
- Conexus uses PROKKA to carry out the annotation step in the hopes of identifying predicted genes and their putative function from the raw genomic input, looping through each of the given files and creating a file specific PROKKA output folder stroring the results
##### Identifying AMR genes 
- In order to filter out the larger number of genes with no known link to AMR Conexus will use the resitance gene identifier (RGI). The software contains a larger database (aka CARD) with genes and specific mutations linked to roles in AMR in peer reviewed papers. Using this information in the published literatture, after identifying predicted resistance genes via a BLASTP search and mutation scanning the software will link the gene to its known resistance mechanisms and drug targets.

- All discovered AMR genes will be added into a pandas dataframe with literature-sourced information of each to allow for future visualisation and analysis
##### Pangenome creation
- Panaroo is used to link all genes stored separatley in different files into an overall pangenome, while additionally clustering highly similar genes into orthologous and paralogous clusters to give an accurate picture of both core and accessory genomes of the bacterial community.
##### Finding key pairings
- Coinfinder takes Panaroo results as input measuring the relative abundance of different gene clusters in the pangenome and calculating expected co-occurances between different pairs, returning those appearing significantly more or less than would be expected by chance.

- Conexus then unpacks these gene clusters, finding their individual genes and the input file they origintated from. After this the program finds the corresponding RGI output for each of the files, collecting a list of all the predicted AMR genes present. Where the list of genes present in the significant co-occuring clusters overlaps with the list of AMR genes present the program considers these genes as key hits retruniung them in hits.csv
# Parameters
|Option      | Description |
|------------|-------------|
|--threads      |Number of threads to run the program (recommended >= 8)     |
|-prefix        |prefix name for output files     |
|-connection    |either associate or disassociate for the type of link between gene pairs    |
|-mode          |Stringency level for which to run panaroo, either strict, moderate or sensitive for details see https://gtonkinhill.github.io/panaroo/#/gettingstarted/params
|-id.           |The postion of the FASTA file's unique identifier starting from 0
# Output
##### Hits.txt
- File conatining all the key gene associations with an identified link to AMR, including the clusters they form part of (as per Panaroo) and the file they originated from
##### Stats.txt
File containing information on:
  - Number of input files
  - Number of files that successfully passed each stage
  - List of files which failed at each particular stage
  - The number of AMR genes dicovered and the number of key pairings
  - The number of key pairings where one was identified in the AMR step
 ##### AMRgenes.csv
 CSV with all predicted AMR genes discovered in the input files, with information on:
  - The originating input file
  - Predicted gene (best hit in the database)
  - Sequence & sequence identity
  - Resistance mechanism
  - Drug target
 ##### Bins.csv
 CSV with overrall information on all of the input bins (files), e.g.:
  - The length of each in bases
  - The number of identified AMR genes in each
  - The number of AMr genes per base 

# Examples
```
python Conexus -prefix Patient1 -threads 10 -connection associate -mode sensitive -id 4
```
