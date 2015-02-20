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
combine_sailfish_output.py --help
```

### Combine [MISO](http://genes.mit.edu/burgelab/miso/) output files

Run this in the directory containing a folder called `miso` which has all your
MISO output.

This assumes a few things:

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
combine_miso_output.py --help
```