# Teleogryllus_occipitalis

## Genome Assembly


### GenomeScope
```

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

```

### QUAST
```

```
### NOVOplasty
```

```
## Genome Annotation

### RepeatModeler, RepeatMasker
```

```
### Trinity, Transdecorder
```

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
