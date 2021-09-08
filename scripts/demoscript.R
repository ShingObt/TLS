library(FORTLS)
library(vroom)
library(lifecycle)
root="/Users/shingo/Dropbox/1-PDProje/TLS"
setwd(root)
las=readLAS("demoData/Block1_plot1.las")
plot(las)
las.norm=normalize("demoData/Block1_plot1.las",save.result = TRUE)
id.tree=treeDetection(las.norm)
temp=metrics.variables(id.tree,
                       plot.parameters=data.frame(radius=30),
                       dir.data = "/Users/shingo/Dropbox/1-PDProje/TLS")
dim(id.tree)
