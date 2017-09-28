
require(rCharts)
## for some reason Rscript this script report an error
## run this script in rstudio under R 3.1.1, where rCharts is installed
## export RSTUDIO_WHICH_R=/gsc/software/linux-x86_64-centos6/R-3.1.1/bin/R
## rstudio

#require(json)
setwd('{{wkdir}}')
## pdf(file.path('{{wkdir}}', "sankey.pdf"), height=100, width=160)

file <-"{{sankey_input}}"
df.sps <- read.table(file, header=TRUE, sep="\t", quote="",
                     as.is=TRUE, stringsAsFactors = FALSE)
workingdata=df.sps
colnames(workingdata)=c('Sector','source','target','value')
mut.num <- nrow(workingdata)

sankeyPlot=function(df){
sankeyPlot <- rCharts$new()
sankeyPlot$setLib('~/projects/development/cfanalyzer/rCharts_d3_sankey-gh-pages/')
sankeyPlot$setTemplate(script = "~/projects/development/cfanalyzer/rCharts_d3_sankey-gh-pages/layouts/chart.html")

sankeyPlot$set(
  data = df,
  nodeWidth = 15,
  nodePadding = 10,
  layout = 32,
  width = 750,
  height = 12*mut.num,
  labelFormat = ".1%"
)
sankeyPlot
}
agg.pri.sps.2.mutations=aggregate(value ~ {{status1sp}} + mutation, df.sps, sum)
colnames(agg.pri.sps.2.mutations)=c('source','target','value')

## We can now generate a single data file combing all source and target data
mutation.2.ref.sps=subset(workingdata,select=c('source','target','value'))
all.source.target.df=rbind(mutation.2.ref.sps, agg.pri.sps.2.mutations)
sankeyPlot(all.source.target.df)

## dev.off()
