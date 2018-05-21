## Reproducible Scripts for the Publication

>Jungnam Cho, Matthias Benoit, Marco Catoni, Hajk-Georg Drost, Anna Brestovitsky, Matthijs Oosterbeek and Jerzy Paszkowski. __Sensitive detection of pre-integration intermediates of LTR retrotransposons in crop plants__. _bioRxiv_ (2018). [doi: http://dx.doi.org/10.1101/317479](https://www.biorxiv.org/content/biorxiv/early/2018/05/15/317479.full.pdf)   

For the computational reproducibility of _de novo_ annotations of LTR retrotransposons we implemented the pipeline `LTRpred`.
`LTRpred` calls the command line tools `suffixerator`, `LTRharvest` (Ellinghaus et al., 2008), and `LTRdigest` (Steinbiss et al., 2009), which are part of the [GenomeTools library](http://genometools.org/) (Gremme et al., 2013) to screen for repeated LTRs, specific sequence motifs such as primer binding sites (PBS), polypurine tract motifs (PPT), and target site duplications (TSD) and for conserved protein domains such as reverse transcriptase (gag), integrase DNA binding domain, integrase Zinc binding domain, RNase H, and the integrase core domain. Subsequently, `LTRpred` implements customized parser functions to import `LTRdigest` output and structures the data in a `tidy data format` (Wickham, 2014) which subsequently enables automation of false positive curation. In a second step, open reading frame (ORF) prediction is performed by a customized wrapper function that runs the command line tool `usearch` (Edgar, 2010). This automated step allows to automatically filter out RTEs that might have conserved protein domains such as an integrase or a reverse transcriptase, but fail to have ORFs and thus are not expressed. In a third step, RTE family clustering is performed using the command line tool `vsearch` (Rognes et al., 2016) which defines family members by >90% sequence homology of the full element to each other. In a fourth step, an automated `hmmer search` (Finn et al., 2011) against the `Dfam` database (Hubley et al., 2016) is performed to assign super-family associations such as Copia or Gypsy by comparing the protein domains of de novo predicted RTEs with already annotated RTEs in the `Dfam` database. For each step, `LTRpred` implements customized parser functions to import `usearch`, `vsearch`, and `Dfam` output and transforms this output in `tidy data format` for subsequent automated false positive curation. In a fifth step, for each predicted element the count and proportion (count divided by element length) of CHH, CHG, CG, and NNN motifs are quantified for the entire element, the 3’ LTR and 5’ LTR separately. In a sixth step, automated false positive curation is performed by the `LTRpred` function `quality.filter()` to conservatively reduce false positive predictions.  


Please install the following R packages before running the reproducible scripts:

```r
install.packages("dplyr")
install.packages("ggplot2")
install.packages("readr")
install.packages("readxl")

source("http://bioconductor.org/biocLite.R")
biocLite('biomartr')

biocLite("devtools")
biocLite("HajkD/LTRpred")
```

Please also make sure that you follow the [INSTALLATION instructions](https://hajkd.github.io/LTRpred/articles/Introduction.html#installation)
of the `LTRpred` package to install all __command line tools__ that `LTRpred` depends on.

## Resources for Annotation

### tRNAs

We retrieved tRNA sequences in `*.fasta` format from the following databases:

- [GtRNAdb](http://gtrnadb2009.ucsc.edu/download.html)
- [plantRNA database](http://plantrna.ibmp.cnrs.fr/plantrna/search/;jsessionid=14635D3979E56DA4F076CE252D4E2078)

We combined tRNA sequences from both databases to have a comprehensive collection of tRNA sequences
specific for each kingdom of life.

### HMM Models

We retrieved the HMM models for protein domain annotation of the region
between de novo predicted LTRs from [Pfam](http://pfam.xfam.org):

  - RNA dependent RNA polymerase: [Overview](http://pfam.xfam.org/clan/CL0027)
      - [RdRP_1](http://pfam.xfam.org/family/PF00680#tabview=tab6)
      - [RdRP_2](http://pfam.xfam.org/family/PF00978#tabview=tab6)
      - [RdRP_3](http://pfam.xfam.org/family/PF00998#tabview=tab6)
      - [RdRP_4](http://pfam.xfam.org/family/PF02123#tabview=tab6)
      - [RVT_1](http://pfam.xfam.org/family/PF00078#tabview=tab6)
      - [RVT_2](http://pfam.xfam.org/family/PF07727#tabview=tab6)
      - [Integrase DNA binding domain](http://pfam.xfam.org/family/PF00552#tabview=tab6)
      - [Integrase Zinc binding domain](http://pfam.xfam.org/family/PF02022#tabview=tab6)
      - [Retrotrans_gag](http://pfam.xfam.org/family/PF03732#tabview=tab6)
      - [RNase H](http://pfam.xfam.org/family/PF00075#tabview=tab6)
      - [Integrase core domain](http://pfam.xfam.org/family/PF00665#tabview=tab6)

## Running `LTRpred`

The following code can be run on a computer with 4 cores. Please be aware that 
computation times might correspond to days due to the genome sizes of the respective species.

For further details about `LTRpred` please consult the [LTRpred: Introduction Vignette](https://hajkd.github.io/LTRpred/articles/Introduction.html).

### Generating _de novo_ LTR retrotransposon annotation for `Arabidopsis thaliana`

```r
library(LTRpred)
# de novo LTR transposon prediction of 'A. thaliana'
LTRpred(
      genome.file = "Athaliana.fa",
      cluster     = TRUE,
      cores       = 4,
      copy.number.est = FALSE,
      minlenltr   = 100,
      maxlenltr   = 5000,
      mindistltr  = 4000,
      maxdistltr  = 30000,
      mintsd      = 3,
      maxtsd      = 20,
      vic         = 80,
      overlaps    = "no",
      xdrop        = 7,
      motifmis    = 1,
      pbsradius   = 60,
      pbsalilen   = c(8,40),
      pbsoffset   = c(0,10),
      quality.filter = TRUE,
      n.orfs      = 0
      )
# import LTRpred output
Athaliana_LTRpred <- read.ltrpred("Athaliana_ltrpred/Athaliana_LTRpred_DataSheet.tsv")
```

### Generating _de novo_ LTR retrotransposon annotation for `Solanum lycopersicum`

```r
library(LTRpred)
# de novo LTR transposon prediction of 'S. lycopersicum'
LTRpred(
      genome.file = "Slycopersicum.fa",
      cluster     = TRUE,
      cores       = 4,
      copy.number.est = FALSE,
      minlenltr   = 100,
      maxlenltr   = 5000,
      mindistltr  = 4000,
      maxdistltr  = 30000,
      mintsd      = 3,
      maxtsd      = 20,
      vic         = 80,
      overlaps    = "no",
      xdrop        = 7,
      motifmis    = 1,
      pbsradius   = 60,
      pbsalilen   = c(8,40),
      pbsoffset   = c(0,10),
      quality.filter = TRUE,
      n.orfs      = 0
      )
# import LTRpred output
Slycopersicum_LTRpred <- read.ltrpred("Slycopersicum_ltrpred/Slycopersicum_LTRpred_DataSheet.tsv")
```

### Generating _de novo_ LTR retrotransposon annotation for `Oryza sativa`

```r
library(LTRpred)
# de novo LTR transposon prediction of 'O. sativa'
LTRpred(
      genome.file = "Osativa.fa",
      cluster     = TRUE,
      cores       = 4,
      copy.number.est = FALSE,
      minlenltr   = 100,
      maxlenltr   = 5000,
      mindistltr  = 4000,
      maxdistltr  = 30000,
      mintsd      = 3,
      maxtsd      = 20,
      vic         = 80,
      overlaps    = "no",
      xdrop        = 7,
      motifmis    = 1,
      pbsradius   = 60,
      pbsalilen   = c(8,40),
      pbsoffset   = c(0,10),
      quality.filter = TRUE,
      n.orfs      = 0
      )
# import LTRpred output
Osativa_LTRpred <- read.ltrpred("Osativa_ltrpred/Osativa_LTRpred_DataSheet.tsv")
```

