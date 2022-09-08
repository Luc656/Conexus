import os
import subprocess
import csv
import pandas as pd
import argparse

headers = []
dataList = []
extraData = []
noProkka = {}
fastas = [file for file in os.listdir() if file.split('.')[-1] == 'fasta']

parser = argparse.ArgumentParser()
parser.add_argument('-id', type=int, help='The postion of the FASTA files unique identifier starting from 0')
parser.add_argument('-prefix', type=str, help='prefix name for output files')
parser.add_argument('--threads', type=int, help='Number of threads to run the program (recommended = 8)')
parser.add_argument('-connection', type=str, help='either associate or disassociate for the type of link between gene pairs')
parser.add_argument('-mode', type=str, help='Stringency level for which to run panaroo, either strict, moderate or sensitive for details see https://gtonkinhill.github.io/panaroo/#/gettingstarted/params')
args = parser.parse_args()
id = args.id
prefix = args.prefix
threads = args.threads
connection = args.connection
mode = args.mode


def bordered(text):
    lines = text.splitlines()
    width = max(len(s) for s in lines)
    res = ['┌' + '─' * width + '┐']
    for s in lines:
        res.append('│' + (s + ' ' * width)[:width] + '│')
    res.append('└' + '─' * width + '┘')
    return '\n'.join(res)


class Conexus:

    def __init__(self, file):
        self.file = file
        self.prefix = '.'.join(self.file.split('.')[:-1])
        self.bin = self.file.split('.')[id]
        self.current = os.listdir()
        self.noProkka = {}
        self.noRGI = []
        self.prokka = self.prefix + '.prokka'
        self.rgi = self.prefix + '.rgi'
        self.faa = self.prokka + '.faa'
        self.gff = self.prokka + '.gff'
        self.RGItxt = self.rgi + '.txt'
        self.prokkaCount = 0
        self.prokkas = [file for files in os.listdir() if files.split('.')[-1] == 'prokka']
        self.fastas = [file for files in os.listdir() if files.split('.')[-1] == 'fasta']

    # annotate all FASTA files with PROKKA
    def annotate(self):
        # check file wasnt already previously annotated
        if self.prokka not in self.current:
            # PROKKA command
            subprocess.call(['prokka', '--outdir', self.prokka, '--prefix', self.prokka, self.file, '--quiet'])
            # move gff fille from output into new file ready for pangenome creation
            # subprocess.call(['cp', f'{self.prokka}/{self.gff}', f'allGFF'])
            print(bordered(f'finished annotating {self.file}'))
        subprocess.call(['cp', f'{self.prokka}/{self.gff}', f'allGFF'])

    # identify all AMR genes with RGI
    def identifyAMR(self):
        if self.rgi not in self.current:
            subprocess.call(['mkdir', self.rgi])
            # take faa file from PROKKA as input
            subprocess.call(['cp', f'./{self.prokka}/{self.faa}', f'./{self.rgi}'])
            os.chdir(self.rgi)
            # RGI command
            subprocess.call(['rgi', 'main', '-i', self.faa, '-o', self.rgi, '-t', 'protein'])
            os.chdir('..')
            print(bordered(f'finished rgi for {self.file}'))

    def check(self):
        self.prokkaCount += 1
        self.noProkka.update({self.file:''})
        os.chdir(self.file)
        for subFile in os.listdir():
            if subFile.split('.')[-1] == 'faa':
                self.noProkka[self.file] = True
        os.chdir('..')

    # Extract the information on all AMR genes found into dataframe
    def extract(self):
        global dataList
        global extraData
        os.chdir(self.rgi)
        # read the RGI output containing AMR genes
        with open(self.RGItxt,'r') as f1:
            tsvFile = csv.reader(f1, delimiter='\t')
            dataList = list(tsvFile)
            # append genes and their information to a dataframe
            for line in dataList[1:]:
                data = [int(self.bin)] + line
                extraData.append(data)
        os.chdir('..')

class ConexusGroup():

    def pangenome(self):
        os.chdir('allGFF')
        subprocess.call(['mkdir', 'results.panaroo'])
        # panaroo command
        os.system(f"panaroo -i *.gff -o results.panaroo --clean-mode {mode} -t {threads}")
        os.chdir('..')

    def coinfind(self):
        os.chdir('allGFF')
        subprocess.call(['mkdir', 'results.coinfinder'])
        # use the presence_absence.csv from panaroo as input
        subprocess.call(['cp', './results.panaroo/gene_presence_absence.csv', './results.coinfinder'])
        os.chdir('./results.coinfinder')
        # alter the file formatting so confinder can read it
        os.system("""sed -e 's/^/"/g' -e 's/$/"/g' -e 's/,/","/g' gene_presence_absence.csv > gene_presence_absence-withquotes.csv""")
        # coinfinder command
        subprocess.call(
            ['coinfinder', '-i', 'gene_presence_absence-withquotes.csv', '-I', '-o', f'{prefix}', f'--{connection}', '-x',
             f'{threads}'])
        os.chdir('../..')

    # Function to find which genes are contained in the significant associating clusters coinfinder highlighted,to then
    # check if any of those were predicted to be AMR genes in the RGI stage, returning significant gene cluster pairings
    # with a link to AMR
    def extractCoin(self):
        # 1. find significant genes/groups for patient
        with open(f'allGFF/results.coinfinder/{prefix}_pairs.tsv', 'r') as f1, open(f'allGFF/results.coinfinder/gene_presence_absence-withquotes.csv', 'r') as f2:
            csv4 = csv.reader(f1, delimiter='\t')
            data1 = list(csv4)
            sigClusters = set()

            for i in range(1, len(data1)):
                # print(data1[i][0],data1[i][1])
                sigClusters.add(data1[i][0])
                sigClusters.add(data1[i][1])
            print(f'significant clusters = {sigClusters}')

            # 2. find which genes in each group for patient
            csv2 = csv.reader(f2, delimiter=',')
            data = list(csv2)
            groups = {}
            sigClustersUnpacked = {}
            link = {}
            headers = list(data[0][3:])

            # get row number of each group
            for j in range(1, len(data)):
                groups[data[j][0]] = j

            # 3. find genes in sig groups
            # 4. match genes in group to their bin
            for j in sigClusters:
                count = 0
                row = groups[j]
                emptys = list(data[row][3:])
                for i in headers:
                    i = i.replace('-', '.')
                    link[i] = emptys[count]
                    count += 1
                sigClustersUnpacked[j] = dict((k, v) for k, v in link.items() if v)

        # os.chdir('../../')

        print(f'significant clusters unpacked: {sigClustersUnpacked}')

        binRGI = []
        for j in sigClustersUnpacked:
            for k, v in sigClustersUnpacked[j].items():
                file = k + '.rgi'
                os.chdir(file)
                with open(f'{file}.txt', 'r') as f3:
                    txt = csv.reader(f3, delimiter='\t', )
                    data3 = list(txt)
                    for line in data3[1:]:
                        print(line[0].split(' ')[0])
                        binRGI.append(line[0].split(' ')[0])
                os.chdir('..')
        hits = set()
        # print(binRGI)
        values = []

        for i in sigClustersUnpacked.values():
            for j in i.values():
                if ';' in j:
                    j = j.split(';')
                    for item in j:
                        values.append(item)
                else:
                    values.append(j)

        for item in values:
            if item in binRGI:
                hits.add(item)
        # print(values)

        print(f'matching hits are {hits}')
        with open('hits.txt', 'w') as file1:
            file1.write('Significant Clusters:')
            file1.write(str(sigClusters))
            file1.write('\n')
            file1.write('Significant clusters unpacked:')
            file1.write(str(sigClustersUnpacked))
            file1.write('Significant clusters with a link to AMR: \n')
            for items in hits:
                file1.writelines([items])
            file1.close()

    def statement(self):
        global noProkka
        noFAA = []

        for i in noProkka:
            if not noProkka[i]:
                noFAA.append(i)

        print(bordered(f'Of {len(fastas)} fasta files there were {len(fastas)-len(noFAA)} annotated properly'))
        print(bordered(f'Files not properly annotated: {noFAA}'))


if __name__=='__main__':

    subprocess.call(['mkdir', 'allGFF'])

    for file in os.listdir():
        if file.split('.')[-1] == 'fasta':
            individual = Conexus(file)
            individual.annotate()
            individual.identifyAMR()
            individual.extract()

    group = ConexusGroup()
    group.pangenome()
    group.coinfind()
    group.extractCoin()

    # setting up dataframe with AMR genes
    headers = [i for i in dataList[0]]
    columns = ['bin'] + headers
    df = pd.DataFrame(columns=columns)
    for i in extraData:
        df.loc[len(df.index)] = i

    df.to_csv('AMRgenes.csv', index=True)
    print(bordered(f'New csv with all {len(df)} AMR genes found is ready'))

    for file in os.listdir():
        if file.split('.')[-1] == 'prokka':
            individual = Conexus(file)
            individual.check()

    group.statement()
