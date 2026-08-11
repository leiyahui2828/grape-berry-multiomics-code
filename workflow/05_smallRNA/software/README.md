# External software used by the historical small-RNA workflow

The supplied analysis references:

```text
FastQC
Cutadapt
Trim Galore
SeqKit
SeqTK
BEDTools
Bowtie
Bowtie2
miRDeep2 0.1.3 (mapper.pl / quantifier.pl path)
miRDP2 pipeline v1.1.4
NCBI BLAST+
psRobot v1.2
TargetFinder wrapper
MUSCLE
trimAl
TreeBest
FastTree
IQ-TREE
MEME
DESeq2
```

Do not copy large third-party software installations into GitHub.

Instead, before the final `v1.0.0` release, create a software-version table and
record the exact executable/package versions used.

Expected local paths used by the organized scripts:

```text
software/mirdeep2-0.1.3/bin/mapper.pl
software/mirdeep2-0.1.3/bin/quantifier.pl
software/miRDP2-v1.1.4/miRDP2-v1.1.4_pipeline.bash
software/psRobot_v1.2/bin/psRobot_tar
software/targetfinder_threads.pl
```
