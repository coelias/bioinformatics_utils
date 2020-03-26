# Carlos Utils! (Yay!)

**bammd5** - Caculates md5 from BAM, having into account only genomic data, not mapping data

```
Usage: bammd5 <file>
```

**fastacmp** - Simple pairwise comparison between two genomes in FASTA format

```
usage: python fastacmp fastafile1 fastafile2 (gzip files allowed)
```

**fastafilter** - Fasta filtering tool and NCBI downloader

```
usage: fastafilter [-h] [-min MINSIZE] [-max MAXSIZE] [-desc DESC]
                   [-extract PATTERN] [-fetch FETCH] [-db DB] [-maxres MAXRES]
                   [-r] [-stats] [-lengths] [-gc] [-o OUTPUT] [--noise NOISE]
                   [--fragment FRAGMENT] [-group GFOLDER]
                   [fastas [fastas ...]]

FastaFilter

positional arguments:
  fastas

optional arguments:
  -h, --help           show this help message and exit
  -min MINSIZE         Minimum Sequence length [Opt]
  -max MAXSIZE         Maximum sequence length [Opt]
  -desc DESC           Filter by description string (regexp) [Opt]
  -extract PATTERN     Extracts a region of the FastA, defined by two
                       subsequences and a number of mismatches (fuzzy regexp).
                       For Example -extract 2,CAGCATA,3,AAAATAGA (Extract
                       subsequences delimited by CAGCATA and AAAATAGA allowing
                       2 and 3 mismatches repectively, if you add a D [DEBUG]
                       (eg. (1,ATG,1,CTG,D) the script will only print the
                       number of matches for every sequence [Opt]
  -fetch FETCH         NCBI query [Opt]
  -db DB               NCBI database [nucleotide]
  -maxres MAXRES       NCBI max number of results [3000]
  -r                   Reverse description filter
  -stats               HTML Statistic information
  -lengths             Display lengths of all sequences
  -gc                  Display GC content for all sequences
  -o OUTPUT            Output file [stdout]
  --noise NOISE        Adds noise to the sequences [muts/base][eg 0.01]
  --fragment FRAGMENT  Return a random fragments [eg: 3,50,10 (three random
                       fragments per read between 50 and 100 bp)]
  -group GFOLDER       Group sequences my -desc regexp and stores them in the
                       provided folder [Opt]
```

**lsum** - du improved, allows to recursively explore dirs and outputs summary by extension

```
$ echo . | lsum 
::::::::::::::::::::::::::::::::::
Total Files: 159
Total size: 231.09 Kb
::::::::::::::::::::::::::::::::::
extension summary:
pyc: 9.30 Kb
fq: 11.99 Kb
py: 44.68 Kb
NOEXTENSION: 163.56 Kb

1 of files (%1 smallest): 1.55 Kb
::::::::::::::::::::::::::::::::::
```

**mapar** - parallel ported to python (parallel sometimes stops!! this doesn't)
```
$ mapar --help
usage: mapar [-h] [-t THREADS] [-o STDOUT] [-e STDERR] [-s]

Carlos parallelizer

optional arguments:
  -h, --help  show this help message and exit
  -t THREADS  Number of processes (defaut #cpu processors)
  -o STDOUT   redirect stdout to a file [and prints progress]
  -e STDERR   redirect stderr to a file [and prints progress]
  -s          supress stdout and stderr and show progress
```

**mappinginfo.py** - tool to generate mapping sattistics from BAM or VCF (see help)

```
$ mappinginfo.py --help
usage: mappinginfo.py [-h] [-s SAMPLEID] [-b BAM] [-v VCF] [-r REFERENCE]
                      [-o HEATOUT] [-cs COVSTATS] [-vs VARIANTSTATS]

Mapping Statistics

optional arguments:
  -h, --help        show this help message and exit
  -s SAMPLEID       SampleId
  -b BAM            Bam file
  -v VCF            Vcf file
  -r REFERENCE      fastafile (if you want to get mutations)
  -o HEATOUT        HeatMaps output file
  -cs COVSTATS      Coverage Stats (bam or VCF)
  -vs VARIANTSTATS  Variant statistics (only vcf)
```

**mgrep** - multiple grep

```
$ mgrep --help
usage: mgrep [-h] -i INPUT [INPUT ...] -r REGEXPS [REGEXPS ...] [-v] [-n]

Multiple regular expression search

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT [INPUT ...]  input files (or stdin -)
  -r REGEXPS [REGEXPS ...]
                        regexp file (or stdin -)
  -v                    Output non matching regexps
  -n                    non regexps, just presence
```

**n50** - quick fasta info (n50 statistic)

```
$ n50 ~/ANALYSIS/centrifugeAnalysys/mockAlignment/high_confidence/Zymo.fa
</Users/cdelojo/ANALYSIS/centrifugeAnalysys/mockAlignment/high_confidence/Zymo.fa>
N50: 218942
Genome Size: 69788472
Contigs: 13586
Largest contig: 2951805
-------------------------
```

**acc2taxid.py** - Tool to transform ncbi accession numbers into taxonomic ids

```
$ acc2taxid.py  --help
usage: acc2taxid.py [-h] [-i INPUT] [-o OUTFILE] [-n]

NCBI Accessions to taxid

optional arguments:
  -h, --help  show this help message and exit
  -i INPUT    File input (one accession per line, default stdin)
  -o OUTFILE  file output (default stdout)
  -n          Use only if massive number of accession numbers, will download
              ncbidbs
```

**lols** - Amazing tool, naive replacement for ls that can list folder with gazillion of files in no time, even before you press enter. ('make lols' will create the executable, usage: 'lols DIRTOLIST')

**clim** - Replacement for column command, works much better. (basically a tablify, but super awesome) 

```
pip install clim
```
**cross** - Tool for crossing data sets like TSV's or CSV. Do not write that python script again, cross does that for you

```
pip install cross
```