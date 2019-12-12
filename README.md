# Teleogryllus_occipitalis

## Genome Assembly


### GenomeScope
```
cat T_occipitalis_1P.fq T_occipitalis_2P.fq > list
jellyfish count -C -m 21 -s 1000000000 -t 40 /home/kataoka/work/DATA/shortread/toccipitalis/list -o reads.jf
jellyfish histo -t 40 reads.jf > reads.histogram
Rscript genomescope.R reads.histogram 31 301 output_dir_kmer_31
```

### MaSuRCA assembler
```
masurca config.txt
./assemble.sh
```

```
DATA
PE= pe 200 116 T_occipitalis_R1.fastq T_occipitalis_R2.fastq
NANOPORE=nanopore.fasta
END
PARAMETERS
EXTEND_JUMP_READS=0
GRAPH_KMER_SIZE = auto
USE_LINKING_MATES = 0
MEGA_READS_ONE_PASS=0
LIMIT_JUMP_COVERAGE = 300
CA_PARAMETERS =  cgwErrorRate=0.15
KMER_COUNT_THRESHOLD = 2
CLOSE_GAPS=1
NUM_THREADS = 40
JF_SIZE = 20000000000
SOAP_ASSEMBLY=0
END
```
### HaploMerger2
```
windowmasker -checkdup true -mk_counts -in /work/korogi/genome/toccipitalis.hybrid.flye.3.fasta -out masking_library.ustat -mem 6500
windowmasker -ustat masking_library.ustat -in /work/korogi/genome/toccipitalis.hybrid.flye.3.fasta  -out  /work/korogi/genome/toccipitalis.hybrid.flye.3.softmasked.fasta -outfmt fasta -dust true
PATH={$PATH}:/home/kataoka/work/HaploMerger2/HaploMerger2_20180603/bin
run_all.batch >run_all.log 2>&1
```

### QUAST
```
quast.py /work/korogi/genome/toccipitalis_genome_final.fa
```
### Trimommatic
```
java -jar trimommatic-0.38.jar PE -threads 40 -phred33 -summary summary.txt -validatePairs T_occipitalus_R1.fastq T_occipitalus_R2.fastq T_occipitalus_1P.fastq T_occipitalus_1U.fastq T_occipitalus_2P.fastq T_occipitalus_2U.fastq LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
```

### NOVOplasty
```
NOVOPlasty3.7.1.pl -c config.txt
```

```
Project:
----------------------
Project name          = Toccipitalis_mitogenome
Type                  = mito
Genome range          = 12000-22000
K-mer                 = 39
Max memory            =
Extended log          = 0
Save assembled reads  = no
Seed Input            = /home/kataoka/work/NOVOPlasty/toccipitalis_seed.fasta
Reference sequence    =
Variance detection    =
Chloroplast sequence  =

Dataset 1:
----------------------
Read Length           = 151
Insert size           =
Platform              = illumina
Single/Paired         = PE
Combined reads        =
Forward reads         = /work/korogi/genome/T_occipitalis_1P_trimmed.fq
Reverse reads         = /work/korogi/genome/T_occipitalis_2P_trimmed.fq

Heteroplasmy:
-----------------------
Heteroplasmy          =
HP exclude list       =
PCR-free              =

Optional:
----------------------
Insert size auto      = yes
Insert range          = 1.9
Insert range strict   = 1.3
Use Quality Scores    =
```
## Genome Annotation

### RepeatModeler, RepeatMasker
```
BuildDatabase -name toccipitalis_repeaet_database -engine ncbi /work/korogi/genome/toccipitalis_genome_final.fa
RepeatModeler -engine ncbi -pa 30 -database toccipitalis_repeaet_database
RepeatMasker /work/korogi/genome/toccipitalis_genome_final.fa -pa 30 -lib ./RM_354200.ThuOct102327132019/consensi.fa.classified -xsmall -s -html
```

### Trinity
```
PATH={$PATH}:/home/kataoka/work/samtools
PATH={$PATH}:/home/kataoka/work/masurca/MaSuRCA-3.3.0/bin
PATH={$PATH}:/home/kataoka/work/bowtie2-2.2.6
PATH={$PATH}:/home/kataoka/work/salmon-latest_linux_x86_64/bin
Trinity --seqType fq --max_memory 100G --left /home/kataoka/work/DATA/transcriptome/Muscle_1.fastq --right /home/kataoka/work/DATA/transcriptome/Muscle_2.fastq --CPU 8
```

### Transdecorder
```
TransDecoder.LongOrfs -t /home/kataoka/work/trinityrnaseq/190806_T.occipitalis/trinity_out_dir/Trinity.fasta
```

### Extraction of complete gene from output of Transdecoder
```
#!/usr/bin/env python
# -*- coding: utf-8 -*-

#PBS -S /bin/env python
#PBS -l nodes=1
#PBS -q yuri_all
#PBS -N compgene_finder
#PBS -j eo
#PBS -m ae

import sys
from Bio import SeqIO
import csv


fasta_in = sys.argv[1]
for record in SeqIO.parse(fasta_in, 'fasta'):
        id_part = record.id
        desc_part = record.description
        seq = record.seq
        if 'type:complete' in desc_part:
                fasta_seq = '>' +  desc_part + '\n' + seq
                print(fasta_seq)
```

### Braker2
```
braker.pl --genome=toccipitalis_genome_final_masked_HeaderEdit.fa --bam=toccipitalis.sort.bam,toccipitalis_temma.sort.bam,toccipitalis_toceanicus.sort.bam --species=191206_1 --cores=8 --AUGUSTUS_ab_initio --softmasking --gff3 --workingdir=/home/kataoka/work/BRAKER/191206_1 --GENEMARK_PATH=gm_et_linux_64 --min_contig=5000 --prot_seq=compgene_translated.fa --prg=gth --ALIGNMENT_TOOL_PATH=gth-1.7.1-Linux_x86_64-64bit/bin --trainFromGth
```

### diamond
```
db=/UNIPROT/swissprot
diamond blastp -d ${db} -q Teleogryllus_occipitalis_UnassignedOG.fasta -p 20 --outfmt 6 -o Unassigned.m6
diamond blastp -d ${db} -q Teleogryllus_occipitalis_UniqOG.fasta -p 20 --outfmt 6 -o Uniq.m6
```

### prosite
```
ps_scan.pl -d prosite.dat Teleogryllus_occipitalis_UniqOG.fasta > Uniq_prosite.txt
ps_scan.pl -d prosite.dat Teleogryllus_occipitalis_UnassignedOG.fasta > Unassigned_prosite.txt
```
