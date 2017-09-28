# expands current version: 1.7.2; use following command to run:
# /gsc/software/linux-x86_64-centos6/R-3.1.1/bin/Rscript {{rscript}}
setwd("{{wkdir}}")
dir.create("{{patient_status}}",
           showWarnings = TRUE,
           recursive = FALSE,
           mode ="0777")

setwd("{{wkdir}}/{{patient_status}}")

library(expands)

runExPANdS("{{snv_input}}",
           "{{cnv_input}}",
           maxScore=2.5,
           max_PM=6,
           min_CellFreq=0.1,
           precision=NA,
           plotF=2,
           snvF="{{patient_status}}.expands",
           maxN=8000,
           region=NA)
setwd("{{wkdir}}")
