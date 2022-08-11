# Conexus
Conexus is a bioinformatics pipeline designed to find key co-occuring gene pairs with a link to AMR. The pipline includes a number of different stages using external software tools to highlight key gene pairings likely under the influence of natural selection in a microbial commuity from raw FASTA input. FASTA files are first annotated then scanned for known AMR genes before being clusted together forming an overall pangeneome where gene pairing frequencies can be measured returning those deemed statically significant. The pipeline also stores information on all discovered AMR genes in the population and their corresponding origin files allowing for in-depth future analysis
# Dependancies 
 - [Prokka](https://github.com/tseemann/prokka)
 - [RGI](https://github.com/arpcard/rgi)
 - [Panaroo](https://gtonkinhill.github.io/panaroo/#/gettingstarted/quickstart)
 - [Coinfinder](http://mcinerneylab.com/software/coinfinder/)
 - [Pandas](https://pandas.pydata.org/)
 - [os](https://docs.python.org/3/library/os.html#module-os)
 - [Subprocess](https://docs.python.org/3/library/subprocess.html#module-subprocess)
 - [CSV](https://docs.python.org/3/library/csv.html#module-csv)
# Installation
```
git clone https://github.com/Luc656/Conexus.git
cd 
chmod
make conda env?
```
# Process
##### Annotation
- text
