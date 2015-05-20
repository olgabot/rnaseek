# API plan for `rnaseek`

The goal for this package is to collect scripts for RNA-seq processing, and to
annotate alternative splicing events. As of now (2015-04-28), the plan is to 
only support skipped exon (SE) and mutually exclusive exon (MXE) splicing 
events, as thes are the largest splice type groups in mammalian systems.

Here are all the splice types that would hopefully be eventually supported: 

[![Types of alternative splicing](http://upload.wikimedia.org/wikipedia/commons/a/ab/Alt_splicing_bestiary2.jpg)](http://commons.wikimedia.org/wiki/File%3AAlt_splicing_bestiary2.jpg)


## Alternative event annotations


### Inputs

Besides the IDs, we also need the splice type, the genome fasta file, and the
genome name. It **should** be easy to grab these given a genome name (i.e. 
`wget` them from ENSEMBL/gencode in the background), but will require some 
slick engineering.

#### For now

- MISO ids, e.g. `'chr1:100:200:+@chr1:300:400:+@chr1:500:600:+'` for an SE 
  event on the positive strand, or `'chr2:700:800:-@chr2:500:600:-@chr2:300:400:-@chr2:100:200:-'`
  for an MXE event on the negative strand

#### Future plans

- Alternatively, the user should be able to supply a `.bed` file, so all these 
  things could be computed for any arbitrary bed file.

### Outputs

- **Can be calculated given raw ids or bed files**
    - Exon length
    - Intron length
- **Can be calculated given genome version (and sequence?)**
    - Exon GC content
    - Exon conservation
    - Conservation of "constitutive," flanking exons
        - "Constitutive" in quotes because they not be truly constitutive relative
           to the gene, but are constitutive relative to the alternative exon.
    - Exon sequence
    - **Require external annotations/**
        - *MaxEntScan*
            - Splice site strength
        - *gffutils* `FeatureDB` made from a gtf file
            - Reading frame
            - Translation
            - Codon adaptation index
        - Overlap with repetitive elements
            - RepBase
        - `RNAhybrid`: microRNA target sites
        - mirbase: microRNA birth sites in the introns
        - Experimental annotation
            - CLIP/ChIP-Seq binding
         - RNA binding protein motifs
            
            
Given a translation, I want to be able to compute things similar to 
[Carvunis et al, Nature (2012)](http://www.ncbi.nlm.nih.gov/pubmed/22722833)

From the supplementary:

> **Sequence properties related to translation**
> The TMHMM46 program was used to predict putative transmembrane regions. The 
> DISOPRED247 program was used to predict disordered regions after removing 
> predicted transmembrane regions and coiled coil regions predicted with COILS48,
> as previously described49. The general average hydropathicity (or GRAVY score)
> was calculated as the arithmetic mean of the sum of the hydropathic indices 
> of each amino acid50, as provided by SGD16. Codon adaptation index was 
> evaluated using the original methodology51. To assess the AUG context 
> optimality, we considered the presence or absence of an adenine in 
> position -3 relative to the AUG start codon52.


- Average hydropathicity
- Predict ...
    - Protein domains
        - *Pfam-A `hmmscan`*
            - From Pfam-A, then to aggregate the domain types:
                - `pfam2go` to get gene ontology annotation of the domains
                - Use the Clans from pfam
            - **For splicing only**:
                - Aggregate changes in predicted protein domains, i.e. 
                  "domain disruption"
            - License?
         - DomPred
    - Disordered protein regions
        - [DISOPRED](http://bioinf.cs.ucl.ac.uk/software_downloads/)
            - License: GNU GPL
        - Disprot? 
            - Use `water` to align to sequences in Disprot fasta file
            - License?
    - Coiled-coil regions
        - [COILS](http://www.ch.embnet.org/software/COILS_form.html)
    - Transmembrane regions
        - [TMHMM](http://www.cbs.dtu.dk/services/TMHMM-2.0/)
        
### API Design constraints

From these descriptions, I see three main possible inputs:

- MISO ids
- Bed files
- Fasta files
    - raw DNA/mRNA
    - protein fasta

As a result, the objects that operate on miso IDs should be independent of the
objects that operate on bed files, which are independent of those that operate
on fasta files. This will allow for maximum portability, and allow the user
to start the annotation process at any point.
