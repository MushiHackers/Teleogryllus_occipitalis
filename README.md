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
PE= pe 200 116 /work/korogi/genome/T_occipitalis_R1.fastq /work/korogi/genome/T_occipitalis_R2.fastq
NANOPORE=/work/korogi/MinION/nanopore.fasta
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

```
