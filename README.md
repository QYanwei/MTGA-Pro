### **MTGA-Pro**


Desciption: 
---
This software is aim to analyze the mitochondiral genome with ultra-high sequencing depth(5000X) data by BGISEQ500 platform. 
The main function contain the rare heteroplasmy variantion calling, haplopgroup Calling & genome assembly for mitochondrial genome.

Development history:
---
MTGA-Pro - MiTochondrial Genome Analysis Pro (2017 Oct 15, compiled Thu 17:10:58 Aug 10 2017)

Auther:
---
Yanwei QI [qiyanwei1@genomics.cn]

Usage:
---

* Argument:

```
-i,--input # File fold stored Paired-end fastq gzip files.
-o,--output # Output file fold.
-s,--sample # Sample name.
-f,--isfilter # If the data has be cleaned by QC, please set with 1. 
-m,--mapping # Mapping step with BWA.
-c,--iscircle # If the D-loop region needs detection the variation carefully, please set with 1.
-t,--isheter # run the step set with 1 and filter the low frequency varation with <5, please combine set with like 1,5.
-a,--isassembly # run the step set with 1 and set the k-mer with a number, default with 1,45.
-p,--ishaplo # run the step set with 1 and set the k-mer with a number, default with 1,45.
-h
```
* Example:

```
Python2.7 MTGA-Pro.py \
--input=/PE_FASTQ/Store/to/PATH/ \
--output=/PATH/to/OUTPUT_DIR/ \
--sample=CL100010855_L01_25 \
--isfilter=0 --mapping=1 --iscircle=1,0 --isheter=0,1 --isassembly=0,45 --ishaplo=0,45
```

Dependencies:
---

+ Python: Version >2.7.0
+ Third-party software: SOAPnuke, bwa, samtools, Haplo, SOAPdenovoTrans.

License:
---
* GNU General Public License v3.0
