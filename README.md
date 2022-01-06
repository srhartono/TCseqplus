### 0. Synopsis

Pipeline downloaded from: https://github.com/robinmeyers/transloc_pipeline

```
TranslocPreprocess.pl meta.txt preprocess/ --indir ./
TranslocWrapper.pl meta.txt preprocess/ results/
TranslocFilter.pl results/SL2003D1/SL2003D1.tlx results/SL2003D1/SL2003D1_result.tlx --filters "f.unaligned=1 f.baitonly=1 f.uncut=G10 f.misprimed=L10 f.freqcut=1 f.largegap=G30 f.mapqual=1 f.breaksite=1 f.sequential=1 f.repeatseq=1 f.duplicate=1"
```

---
  
### 1. Run TranslocPreprocess.pl which preprocess the fq files

TranslocPreprocess.pl <meta.txt> <preprocess output folder> --indir <input folder with fastq files>

---
  
### 2. IIRC TranslocWrapper.pl is the mapping part of the software which will put results in the results/ directory

TranslocWrapper.pl <meta.txt> <preprocess output folder from (1)> <results output folder>

---
  
### 3. Then do final filter. Please see TCseq pipeline for more info (https://github.com/robinmeyers/transloc_pipeline)

For example, this is for SL2003D1

`TranslocFilter.pl <resulting tlx file from (2)> <filtered tlx output file> --filters "<filters>"`

---
  
### 4. meta.txt (tab delimited):

```
Library Assembly    Chr Start   End Strand  MID Primer  Adapter Description
SL2003W1    mm9 chr12   114664888   114664910   -   AGTCAA  CACACAAAGACTCTGGACCTC   CCACGCGTGCTCTACA    Smu_as_bait_site
SL2003S1    mm9 chr12   114664888   114664910   -   GTCCGC  CACACAAAGACTCTGGACCTC   CCACGCGTGCTCTACA    Smu_as_bait_site
SL2003R1    mm9 chr12   114664888   114664910   -   CCGTCC  CACACAAAGACTCTGGACCTC   CCACGCGTGCTCTACA    Smu_as_bait_site
SL2003W2    mm9 chr12   114664888   114664910   -   ATCACG  CACACAAAGACTCTGGACCTC   CCACGCGTGCTCTACA    Smu_as_bait_cut
SL2003S2    mm9 chr12   114664888   114664910   -   CGATGT  CACACAAAGACTCTGGACCTC   CCACGCGTGCTCTACA    Smu_as_bait_cut 
SL2003R2    mm9 chr12   114664888   114664910   -   TTAGGC  CACACAAAGACTCTGGACCTC   CCACGCGTGCTCTACA    Smu_as_bait_cut
SL2003W3    mm9 chr12   114664888   114664910   -   ACAGTG  CACACAAAGACTCTGGACCTC   CCACGCGTGCTCTACA    Smu_as_bait_site
SL2003S3    mm9 chr12   114664888   114664910   -   GCCAAT  CACACAAAGACTCTGGACCTC   CCACGCGTGCTCTACA    Smu_as_bait_site
SL2003R3    mm9 chr12   114664888   114664910   -   CAGATC  CACACAAAGACTCTGGACCTC   CCACGCGTGCTCTACA    Smu_as_bait_site<run TCseq pipeline until fliter.tlx>
```

---
  
### 5. To post-process and create figure of % read with mutation, junction distance etc:


- 1st run Folder: `/group/stella/Work/Data/Fastq/200323_SLIMS0323_BC_mm_TCseq_JackieBarlowTCseq/2_Peak/results`
- 2nd run Folder: `/group/stella/Work/Data/Fastq/200903_SLIMS4168_BC_mm_TCseq_JackieBarlowTCseq/2_Peak/results`

-m 1 means use the metafile in the 1st run fastq folder `/group/stella/Work/Data/Fastq/200323_SLIMS0323_BC_mm_TCseq_JackieBarlowTCseq/0_Fastq/0_META.TXT`
-m 2 means use the metaFile in the 2nd run fastq folder `/group/stella/Work/Data/Fastq/200903_SLIMS4168_BC_mm_TCseq_JackieBarlowTCseq/0_Fastq/0_META.TXT`

There are 2 methods that should give similar result, which was created to "double check" each other.

Synopsis of run:

```
cd <folder_of_results>

# A. parse filter.tlx file and create a table of read and their mutations
3_Mutation.pl -m <1 or 2> -i <filter.tlx file>

# B. calculate distance for each read
5_Dist.pl -i <filter.tlx file>

# C. parse table from 3_Mutation.pl and rearrange for R script

# C.1. Method 1 (This also create junc distance figure)
4_FixStats.pl -i . -e "*.final.tsv.all"
run_Rscript.pl graph_and_juncdist.R

# C.2. Method 2:
both.pl -i ./
run_Rscript.pl Jackie_JuncDist.R
run_Rscript.pl JackieTCseq_BigTable.R  
run_Rscript.pl JackieTCseq_NEW_combined.R  

# results_*.pdf = barplots of mutations
# results3_*.pdf = junc distance figure and barplots
```
