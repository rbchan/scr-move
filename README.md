# README

Data, code, and other files associated with the manuscript "Modeling abundance, distribution, movement, and space  use with camera and telemetry data" hopefully to appear in Ecology as part of the special issue on joint spatial capture-recapture movement models.

## Contents

- [camera_telemetry.tex](camera_telemetry.tex) Manuscript

- [ecology.bst](ecology.bst) Bibliography style file for Ecology

- [mybib.bib](mybib.bib) BibTeX bibliography

- [R/](R/) Directory with R scripts, shell scripts, and data
  * [deer_scr_telem.RData](R/deer_scr_telem.RData) Compressed R data file with deer camera and telemetry data used in the analysis
  * [cluster](R/cluster) Directory with R scripts and shell files for execution on the [GACRC cluster](https://wiki.gacrc.uga.edu/wiki/Georgia_Advanced_Computing_Resource_Center)
    + [run_scr_DA.sh](R/cluster/run_scr_DA.sh) Shell script to run `scr_DA.R` on the GACRC cluster
    + [scr_DA.R](R/cluster/scr_DA.R) R script to fit basic SCR model without movement
	+ Additional shell and R scripts are similar to the two above, but for the models with movement
  * [R/sims](R/sims) Directory with R scripts and shell files for running the simulation study
  
- [figs/](figs/) Directory with figures used in the manuscript

- [scrmove/](scrmove/) R package with MCMC samplers

- [supp/](supp/) Directory with supplemental information
  * [Appendix-S1](supp/Appendix-S1.pdf) 
  * [Appendix-S2](supp/Appendix-S2.pdf)


# DOI

[![DOI](https://zenodo.org/badge/355622082.svg)](https://zenodo.org/badge/latestdoi/355622082)

