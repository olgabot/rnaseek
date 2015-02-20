# rnaseek
Scripts and modules for parsing and aggregating output files from 
RNA-sequencing data. This was initially created for personal use and there
is no warranty or guarantee.

## Installation

Clone the repository, navigate to that directory, and use `pip` to install.

```
git clone https://github.com/olgabot/rnaseek.git
cd rnaseek
pip install .
```

## Usage

This is the functionality so far.

### Combine [Sailfish](http://www.cs.cmu.edu/~ckingsf/software/sailfish/) output files

Run this from the directory where you have folders containing Sailfish output.
To keep the user informed, this program outputs progress in addition to 
writing files.
For simplicity, this assumes a few things:

- You created your Sailfish index using Gencode genes and ERCC spike ins:

```
wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.pc_transcripts.fa.gz
wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.lncRNA_transcripts.fa.gz
wget http://tools.lifetechnologies.com/downloads/ERCC92.fa
zcat gencode.v19.pc_transcripts.fa.gz gencode.v19.lncRNA_transcripts.fa.gz | cat - ERCC92.fa > gencode.v19.pc_lncRNA_transcripts.ercc_spikein.fa
sailfish index --kmerSize 31 --threads 8 --transcripts gencode.v19.pc_lncRNA_transcripts.ercc_spikein.fa --out gencode.v19.pc_lncRNA_transcripts.ercc_spikein.fa_sailfish_index_k31
```
 
- Your samples' folders are named: `<sample_id>.anythingelse.sailfish`
- You want Transcripts Per Million (TPM) (See 
  [this blog post](https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/)
  for a great review of RNA-seq expression units)
- You want to use the `quant_bias_corrected.sf` file, not `quant.sf`
- You want to collapse the expression values for all the transcripts of a gene,
  into gene-level expression


```
combine_sailfish_output.py
```

Get help/documentation using

```
$ combine_sailfish_output.py --help
usage: combine_sailfish_output.py [-h] [-g GLOB_COMMAND] [-o OUT_DIR]
                                  [-n N_PROGRESS]

Combine sailfish output files

optional arguments:
  -h, --help            show this help message and exit
  -g GLOB_COMMAND, --glob-command GLOB_COMMAND
                        Where to find sailfish output directories. Default is
                        folders in the current directory whose names end with
                        "sailfish"
  -o OUT_DIR, --out-dir OUT_DIR
                        Where to output the combined matrices. Does not need
                        to exist already. Default is to create a folder called
                        "combined_input"
  -n N_PROGRESS, --n-progress N_PROGRESS
                        Number of files to show per iterative progress, e.g.
                        10/58 files completed, 20/58 files completed. Can
                        increase this ifyou have thousands of files, or
                        decrease to1 if you only have a few
```


#### Example output

```
$ combine_sailfish_output.py
Reading 246 of sailfish's quant_bias_corrected.sf files ...
        10/246 files read
        20/246 files read
        30/246 files read
        40/246 files read
        50/246 files read
        60/246 files read
        70/246 files read
        80/246 files read
        90/246 files read
        100/246 files read
        110/246 files read
        120/246 files read
        130/246 files read
        140/246 files read
        150/246 files read
        160/246 files read
        170/246 files read
        180/246 files read
        190/246 files read
        200/246 files read
        210/246 files read
        220/246 files read
        230/246 files read
        240/246 files read
        Done.
Separating out spike-ins from regular genes ...
        Done.
Summing TPM expression of all transcripts in a gene ...
        Done.
Writing output files ...
        Wrote /oasis/tscc/scratch/obotvinnik/projects/singlecell_pnms/analysis/singlecell_pnms_pe_v6/combined_output/tpm_genes.csv
        Wrote /oasis/tscc/scratch/obotvinnik/projects/singlecell_pnms/analysis/singlecell_pnms_pe_v6/combined_output/tpm_spikein.csv
Done, son.
```

### Combine [MISO](http://genes.mit.edu/burgelab/miso/) output files

Run this in the directory containing a folder called `miso` which has all your
MISO output. This function outputs progress to keep the user informed on how
many files the program has processed, whether the file only had a header and
no output, splicing event filtering, making a [pivot table](http://en.wikipedia.org/wiki/Pivot_table)
of the psi scores, and writing the files.

For simplicity, this assumes a few things:

- Your miso summary files are located,

```
miso/<sample_id>/<splice_type>/summary/<splice_type>.summary
```

- You want to filter your MISO summary files for only events with confidence
  intervals < 0.5 and at least 10 reads counted by MISO to be unique to one of
  the two possible isoforms.

Run it using

```
combine_miso_output.py
```

Get help/documentation using

```
$ combine_miso_output.py --help
usage: combine_miso_output.py [-h] [-g GLOB_COMMAND] [-o OUT_DIR]
                              [-n N_PROGRESS] [--ci-max CI_MAX]
                              [--per-isoform-reads-min PER_ISOFORM_READS_MIN]
                              [--downsampled]

Combine sailfish output files

optional arguments:
  -h, --help            show this help message and exit
  -g GLOB_COMMAND, --glob-command GLOB_COMMAND
                        Where to find sailfish output directories. The default
                        assumes you have your MISO summary files are located
                        as so: <current_directory>/miso/<sample_id>/<splice_ty
                        pe>/summary/<splice_type>.miso_summary
  -o OUT_DIR, --out-dir OUT_DIR
                        Where to output the combined matrices. Does not need
                        to exist already. Default is to create a folder called
                        "combined_input"
  -n N_PROGRESS, --n-progress N_PROGRESS
                        Number of files to show per iterative progress, e.g.
                        10/58 files completed, 20/58 files completed. Can
                        increase this ifyou have thousands of files, or
                        decrease to 1 if you only have a few
  --ci-max CI_MAX       Used for filtering. Maximum size of the confidence
                        interval for the 'percent-spliced-in' aka PSI score of
                        a splicing event
  --per-isoform-reads-min PER_ISOFORM_READS_MIN
                        Used for filtering. Minimum number of readsthat MISO
                        deemed unique to one isoform, for that splicing event
                        to pass filtering.
  --downsampled         If given, then assumed that the samples were
                        downsampled at certain numbers of reads, and assumes
                        your summary files are located in,<current_directory>/
                        miso/<sample_id>_prob<p>_iter<i>/<splice_type>/summary
                        /<splice_type>.miso_summary
```

#### Example output

```
$ combine_miso_output.py -n 100
Reading 1972 MISO summary files ...
        100/1972 files attempted to read
        200/1972 files attempted to read
        300/1972 files attempted to read
        400/1972 files attempted to read
        Only found header and an empty table for ./miso/N2_06_R1/TANDEMUTR/summary/TANDEMUTR.miso_summary
        Only found header and an empty table for ./miso/N2_06_R1/ALE/summary/ALE.miso_summary
        Only found header and an empty table for ./miso/N2_06_R1/AFE/summary/AFE.miso_summary
        Only found header and an empty table for ./miso/N2_06_R1/A3SS/summary/A3SS.miso_summary
        Only found header and an empty table for ./miso/N2_06_R1/MXE/summary/MXE.miso_summary
        Only found header and an empty table for ./miso/N2_06_R1/SE/summary/SE.miso_summary
        Only found header and an empty table for ./miso/N2_06_R1/A5SS/summary/A5SS.miso_summary
        Only found header and an empty table for ./miso/N2_06_R1/RI/summary/RI.miso_summary
        Only found header and an empty table for ./miso/N2_03_R1/TANDEMUTR/summary/TANDEMUTR.miso_summary
        Only found header and an empty table for ./miso/N2_03_R1/ALE/summary/ALE.miso_summary
        Only found header and an empty table for ./miso/N2_03_R1/AFE/summary/AFE.miso_summary
        Only found header and an empty table for ./miso/N2_03_R1/A3SS/summary/A3SS.miso_summary
        Only found header and an empty table for ./miso/N2_03_R1/MXE/summary/MXE.miso_summary
        Only found header and an empty table for ./miso/N2_03_R1/SE/summary/SE.miso_summary
        Only found header and an empty table for ./miso/N2_03_R1/A5SS/summary/A5SS.miso_summary
        Only found header and an empty table for ./miso/N2_03_R1/RI/summary/RI.miso_summary
        500/1972 files attempted to read
        600/1972 files attempted to read
        700/1972 files attempted to read
        800/1972 files attempted to read
        Only found header and an empty table for ./miso/P8_09_R1/TANDEMUTR/summary/TANDEMUTR.miso_summary
        Only found header and an empty table for ./miso/P8_09_R1/ALE/summary/ALE.miso_summary
        Only found header and an empty table for ./miso/P8_09_R1/AFE/summary/AFE.miso_summary
        Only found header and an empty table for ./miso/P8_09_R1/A3SS/summary/A3SS.miso_summary
        Only found header and an empty table for ./miso/P8_09_R1/MXE/summary/MXE.miso_summary
        Only found header and an empty table for ./miso/P8_09_R1/SE/summary/SE.miso_summary
        Only found header and an empty table for ./miso/P8_09_R1/A5SS/summary/A5SS.miso_summary
        Only found header and an empty table for ./miso/P8_09_R1/RI/summary/RI.miso_summary
        900/1972 files attempted to read
        Only found header and an empty table for ./miso/M3_10_R1/TANDEMUTR/summary/TANDEMUTR.miso_summary
        Only found header and an empty table for ./miso/M3_10_R1/ALE/summary/ALE.miso_summary
        Only found header and an empty table for ./miso/M3_10_R1/AFE/summary/AFE.miso_summary
        Only found header and an empty table for ./miso/M3_10_R1/A3SS/summary/A3SS.miso_summary
        Only found header and an empty table for ./miso/M3_10_R1/MXE/summary/MXE.miso_summary
        Only found header and an empty table for ./miso/M3_10_R1/SE/summary/SE.miso_summary
        Only found header and an empty table for ./miso/M3_10_R1/A5SS/summary/A5SS.miso_summary
        Only found header and an empty table for ./miso/M3_10_R1/RI/summary/RI.miso_summary
        1000/1972 files attempted to read
        1100/1972 files attempted to read
        Only found header and an empty table for ./miso/P8_07_R1/TANDEMUTR/summary/TANDEMUTR.miso_summary
        Only found header and an empty table for ./miso/P8_07_R1/ALE/summary/ALE.miso_summary
        Only found header and an empty table for ./miso/P8_07_R1/AFE/summary/AFE.miso_summary
        Only found header and an empty table for ./miso/P8_07_R1/A3SS/summary/A3SS.miso_summary
        Only found header and an empty table for ./miso/P8_07_R1/MXE/summary/MXE.miso_summary
        Only found header and an empty table for ./miso/P8_07_R1/SE/summary/SE.miso_summary
        Only found header and an empty table for ./miso/P8_07_R1/A5SS/summary/A5SS.miso_summary
        Only found header and an empty table for ./miso/P8_07_R1/RI/summary/RI.miso_summary
        1200/1972 files attempted to read
        1300/1972 files attempted to read
        Only found header and an empty table for ./miso/MSA_25_R1/TANDEMUTR/summary/TANDEMUTR.miso_summary
        Only found header and an empty table for ./miso/MSA_25_R1/ALE/summary/ALE.miso_summary
        Only found header and an empty table for ./miso/MSA_25_R1/AFE/summary/AFE.miso_summary
        Only found header and an empty table for ./miso/MSA_25_R1/A3SS/summary/A3SS.miso_summary
        Only found header and an empty table for ./miso/MSA_25_R1/MXE/summary/MXE.miso_summary
        Only found header and an empty table for ./miso/MSA_25_R1/SE/summary/SE.miso_summary
        Only found header and an empty table for ./miso/MSA_25_R1/A5SS/summary/A5SS.miso_summary
        Only found header and an empty table for ./miso/MSA_25_R1/RI/summary/RI.miso_summary
        Only found header and an empty table for ./miso/M6_04_R1/TANDEMUTR/summary/TANDEMUTR.miso_summary
        Only found header and an empty table for ./miso/M6_04_R1/ALE/summary/ALE.miso_summary
        Only found header and an empty table for ./miso/M6_04_R1/AFE/summary/AFE.miso_summary
        Only found header and an empty table for ./miso/M6_04_R1/A3SS/summary/A3SS.miso_summary
        Only found header and an empty table for ./miso/M6_04_R1/MXE/summary/MXE.miso_summary
        Only found header and an empty table for ./miso/M6_04_R1/SE/summary/SE.miso_summary
        Only found header and an empty table for ./miso/M6_04_R1/A5SS/summary/A5SS.miso_summary
        Only found header and an empty table for ./miso/M6_04_R1/RI/summary/RI.miso_summary
        1400/1972 files attempted to read
        1500/1972 files attempted to read
        1600/1972 files attempted to read
        Only found header and an empty table for ./miso/M2nd_31_R1/TANDEMUTR/summary/TANDEMUTR.miso_summary
        Only found header and an empty table for ./miso/M2nd_31_R1/ALE/summary/ALE.miso_summary
        Only found header and an empty table for ./miso/M2nd_31_R1/AFE/summary/AFE.miso_summary
        Only found header and an empty table for ./miso/M2nd_31_R1/A3SS/summary/A3SS.miso_summary
        Only found header and an empty table for ./miso/M2nd_31_R1/MXE/summary/MXE.miso_summary
        Only found header and an empty table for ./miso/M2nd_31_R1/SE/summary/SE.miso_summary
        Only found header and an empty table for ./miso/M2nd_31_R1/A5SS/summary/A5SS.miso_summary
        Only found header and an empty table for ./miso/M2nd_31_R1/RI/summary/RI.miso_summary
        Only found header and an empty table for ./miso/P8_12_R1/TANDEMUTR/summary/TANDEMUTR.miso_summary
        Only found header and an empty table for ./miso/P8_12_R1/ALE/summary/ALE.miso_summary
        Only found header and an empty table for ./miso/P8_12_R1/AFE/summary/AFE.miso_summary
        Only found header and an empty table for ./miso/P8_12_R1/A3SS/summary/A3SS.miso_summary
        Only found header and an empty table for ./miso/P8_12_R1/MXE/summary/MXE.miso_summary
        Only found header and an empty table for ./miso/P8_12_R1/SE/summary/SE.miso_summary
        Only found header and an empty table for ./miso/P8_12_R1/A5SS/summary/A5SS.miso_summary
        Only found header and an empty table for ./miso/P8_12_R1/RI/summary/RI.miso_summary
        1800/1972 files attempted to read
        Only found header and an empty table for ./miso/N3_08_R1/TANDEMUTR/summary/TANDEMUTR.miso_summary
        Only found header and an empty table for ./miso/N3_08_R1/ALE/summary/ALE.miso_summary
        Only found header and an empty table for ./miso/N3_08_R1/AFE/summary/AFE.miso_summary
        Only found header and an empty table for ./miso/N3_08_R1/A3SS/summary/A3SS.miso_summary
        Only found header and an empty table for ./miso/N3_08_R1/MXE/summary/MXE.miso_summary
        Only found header and an empty table for ./miso/N3_08_R1/SE/summary/SE.miso_summary
        Only found header and an empty table for ./miso/N3_08_R1/A5SS/summary/A5SS.miso_summary
        Only found header and an empty table for ./miso/N3_08_R1/RI/summary/RI.miso_summary
        1900/1972 files attempted to read
        Done.
Merging all 1892 MISO summaries into a gigantic one ...
        Done.
Writing raw MISO summary files ...
        Wrote /oasis/tscc/scratch/obotvinnik/projects/singlecell_pnms/analysis/singlecell_pnms_pe_v6/combined_output/miso_summary_raw.csv
Filtering MISO summaries with ci_max=0.5, per_isoform_counts=10 ...
 -906864 events removed with poor confidence (ci >0.50)
 -787610 events removed with too few reads are unique to individual isoforms (n < 10)
        Done.
Writing filtered MISO summary files ...
        Wrote /oasis/tscc/scratch/obotvinnik/projects/singlecell_pnms/analysis/singlecell_pnms_pe_v6/combined_output/miso_summary_filtered.csv
Creating ((event_name, splice_type), samples) PSI matrix ...
        Wrote /oasis/tscc/scratch/obotvinnik/projects/singlecell_pnms/analysis/singlecell_pnms_pe_v6/combined_output/miso_summary_filtered.csv
```