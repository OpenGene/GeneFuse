# GeneFuse
A tool to detect and visualize target gene fusions by scanning FASTQ files directly. This tool accepts FASTQ files and reference genome as input, and outputs detected fusion results in TEXT, JSON and HTML formats.

# Take a quick glance of the informative report
* Sample HTML report: http://opengene.org/GeneFuse/report.html
* Sample JSON report: http://opengene.org/GeneFuse/report.json
* Dataset for testing: http://opengene.org/dataset.html

# Get genefuse program
## download binary
This binary is only for Linux systems, http://opengene.org/GeneFuse/genefuse
```shell
# this binary was compiled on CentOS, and tested on CentOS/Ubuntu
wget http://opengene.org/GeneFuse/genefuse
chmod a+x ./genefuse
```
## or compile from source
```shell
# get source (you can also use browser to download from master or releases)
git clone https://github.com/OpenGene/genefuse.git

# build
cd genefuse
make

# Install
sudo make install
```

# Usage
You should provide following arguments to run genefuse
* the reference genome fasta file, specified by `-r` or `--ref=`
* the fusion setting file, specified by `-f` or `--fusion=`
* the fastq file(s), specified by `-1` or `--read1=` for single-end data. If dealing with pair-end data, specify the read2 file by `-2` or `--read2=`
* use `-h` or `--html=` to specify the file name of HTML report
* the plain text result is directly printed to STDOUT, you can pipe it to a file using a `>`

## Example
```shell
genefuse -r hg19.fasta -f genes/druggable.hg19.csv -1 R1.fq.gz -2 R2.fq.gz -h report.html > result
```

## Reference genome
The reference genome should be a single whole FASTA file containg all chromosome data. This file shouldn't be compressed. For human data, typicall `hg19/GRch37` or `hg38/GRch38` assembly is used, which can be downloaded from following sites:
* `hg19/GRch37`: ftp://ftp.ncbi.nlm.nih.gov/sra/reports/Assembly/GRCh37-HG19_Broad_variant/Homo_sapiens_assembly19.fasta
* `hg38/GRch38`: http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz  Remember to decompress hg38.fa.gz since it is gzipped and is not supported currently.

## Fusion file
The fusion file is a list of coordinated target genes together with their exons. A sample is:
```CSV
>EML4_ENST00000318522.5,chr2:42396490-42559688
1,42396490,42396776
2,42472645,42472827
3,42483641,42483770
4,42488261,42488434
5,42490318,42490446
...

>ALK_ENST00000389048.3,chr2:29415640-30144432
1,30142859,30144432
2,29940444,29940563
3,29917716,29917880
4,29754781,29754982
5,29606598,29606725
...
```
The coordination system should be consistent with the reference genome.  
### Fusion files provided in this package
Four fusion files are provided with `genefuse`:
1. `genes/druggable.hg19.csv`: all druggable fusion genes based on `hg19/GRch37` reference assembly.
2. `genes/druggable.hg38.csv`: all druggable fusion genes based on `hg38/GRch38` reference assembly.
3. `genes/cancer.hg19.csv`: all COSMIC curated fusion genes (http://cancer.sanger.ac.uk/cosmic/fusion) based on `hg19/GRch37` reference assembly.
4. `genes/cancer.hg38.csv`: all COSMIC curated fusion genes (http://cancer.sanger.ac.uk/cosmic/fusion) based on `hg38/GRch38` reference assembly.

Notes:
* `genefuse` runs almost ~5X faster with `druggable` genes than `cancer` genes, since `druggable` genes are only a small subset of `cancer` genes. Use this one if you only care about the fusion related personalized medicine for cancers.
* The `cancer` genes should be enough for most cancer related studies, since all COSMIC curated fusion genes are included.
* If you want to create a custom gene list, please follow the instructions given on next section.
### Create a fusion file based on hg19 or hg38
If you'd like to create a custom fusion file, you can use `scripts/gen_fusion_file.jl`, which is based on the Julia library `OpenGene.jl` to generate the fusion file you want.   
You should prepare a file containing all genes you want, seperated by `space` or `line break`. Please note that `comma` is not supported. Each gene should be the HGNC standard name.  
By default, the primary transcript (named as GENE_001) will be used. But you can specify the transcript by add `_TranscriptId` to the gene. For example: use `CD74_ENST00000009530` to specify the transcript of `CD74`.   
When the gene list file (`genes.txt`) is prepared, you can used following command to generate a fusion file (`fusion.csv`):
```shell
julia scripts/gen_fusion_file.jl -r hg19 -g genes.txt -f fusion.csv
```
The reference genome is specified by `-r` option, which can be hg19/GRch37/GRch38.

# HTML report
GeneFuse can generate very informative and interactive HTML pages to visualize the fusions with following information:
* the fusion genes, along with their transcripts.
* the inferred break point with reference genome coordinations.
* the inferred fusion protein, with all exons and the transcription direction.
* the supporting reads, with all bases colorized according to their quality scores.
* the number of supporting reads, and how many of them are unique (the rest may be duplications)
## A HTML report example
![image](http://www.opengene.org/GeneFuse/eml4alk.png)  
See the HTML page of this picture: http://opengene.org/GeneFuse/report.html

# All options
```
options:
  -1, --read1       read1 file name (string)
  -2, --read2       read2 file name (string [=])
  -f, --fusion      fusion file name, in CSV format (string)
  -r, --ref         reference fasta file name (string)
  -u, --unique      specify the least supporting read number is required to report a fusion, default is 2 (int [=2])
  -d, --deletion    specify the least deletion length of a intra-gene deletion to report, default is 50 (int [=50])
  -h, --html        file name to store HTML report, default is genefuse.html (string [=genefuse.html])
  -j, --json        file name to store JSON report, default is genefuse.json (string [=genefuse.json])
  -t, --thread      worker thread number, default is 4 (int [=4])
  -?, --help        print this message
```
