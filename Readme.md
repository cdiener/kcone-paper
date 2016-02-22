![paper status](https://img.shields.io/badge/paper-under_review-yellow.svg)
[![protocol status](https://img.shields.io/badge/protocol-online-green.svg)](https://cdiener.github.io/kcone-paper)

## Materials and Protocol for the k-cone paper

Full citation still pending. This repository contains all materials in order
to reproduce *all* results from our publication. You can either create a detailed
computational protocol from an [Rmarkdown](http://rmarkdown.rstudio.com/) file
or just run an [R](http://r-project.org) script to reproduce the results.

The output can be found at https://cdiener.github.io/kcone-paper.

The publication is currently under review, and the results will be citeable upon
acceptance. 

For the source code of the `dycone` package please see 
https://github.com/cdiener/dycone.

### Overview of the provided files

This repository provides the following input files:

- metabolome measurements for HaCaT and HeLa (metabolome.csv)
- model for central carbon metabolism with annotations (reactions.csv) 
- a mapping of the model metabolites to [KEGG](http://www.genome.jp/kegg/) 
  compound and [HMDB](http://hmdb.ca) IDs (id_map.csv)
- a list of 58 microrray samples from the [GEO](http://www.ncbi.nlm.nih.gov/geo/)
  database annotated with the corresponding cell line (ge_samples.csv)
- the model in dycone's CSV format (reactions.csv)
- the model in SBML format (cemet.xml)
  
Additionally, upon first run we will generate some intermediate data files that 
include results for the most time consuming steps in the analysis. This way the 
first run will take about 5h on a machine with 6+ cores, but repeated runs will 
finish in about 4 minutes. The intermediate files cover the following protocol 
steps:

- web scraping of missing metabolite concentrations from HMDB (scraped_concs.Rd)
- calculation of the polytope basis for the model (basis.Rd)
- stability analysis for all of the 6 basis (stab.Rd)

For the sake of completeness we also provide the pathway image we use in the
paper as [Inkscape](https://inkscape.org) SVGs in the images folder. 
**This image is licensed under the [CC-BY 4.0](https://creativecommons.org/licenses/by/4.0/).**

### Reproducing the paper

In order create HTML and PDF outputs of the protocol follow the initial 
[installation intructions](http://cdiener.github.io/kcone-paper) (section 
"Installation")and render the entire Rmarkdown protocol with

```bash
Rscript -e "rmarkdown::render('protocol.Rmd', output_format='all')"
```

HTML output is the same as found on https://cdiener.github.io/kcone-paper.

In case you want to run the analysis in a "quiet" manner without producing the
HTML and PDF reports but yielding all other output (tables and figure files) you
can also run the script version via

```bash
Rscript -e "source('protocol.R')"
```

This is useful if you want to perform small changes and see how this alters the
results.

