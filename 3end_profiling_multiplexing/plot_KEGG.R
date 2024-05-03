#Chai Ruochen
#2021-05-19
#This script is to visualize gene expression of RAINBOW-seq data in interested KEGG pathways
#input: gene expression matrix and path csv
#Output: png XML
#For multiple mapping node, z-score and log2 scale methods use different parameters to calculate node summary. Z-score: max.abs, log2: mean 

#To execute:
#Rscript plot_KEGG.R --expression <expression.csv> --poi <path.csv> --output <output directory>
#For help:
#Rscript plot_KEGG.R --help



library(argparser, quietly = T)
options(connectionObserver = NULL)
library(pathview)
library(org.EcK12.eg.db)

p <- arg_parser("plot_KEGG")
#add command line arguments
p <- add_argument(p,"--expression",help ="gene expression matrix file path",
                  type="character",default = NULL)
p <- add_argument(p,"--poi",help ="path of interest, a csv file with column1 for KEGGID, column2 for KEGGNAME",
                  type="character",default = NULL)
p <- add_argument(p,"--output",help ="ouput directory", 
                  type="character",default = '.')
p <- add_argument(p,"--cutoff",help ="a expression threshold, if a gene with expression less than cutoff, expression is set to zero.", 
                  type="numeric",default = 0)
p <- add_argument(p,"--log",help ="output log2scaled expression on pathview,--log & zscore at least one should be set", 
                  flag=TRUE)
p <- add_argument(p,"--zscore",help ="output per gen z-scaled expression on pathview,--log & zscore at least one should be set", 
                  flag=TRUE)




argv <- parse_args(p)

expression.path = argv$expression
poi = argv$poi
OUTDIR = argv$output
c = argv$cutoff
l = argv$log
z = argv$zscore

if(l+z == 0){
  stop("please choose at least one mode to scale data, --log, --zscore")
}

cutexpression = function(x){ifelse(x>=c,x,0)}


#load expression matrix
expression = read.csv(expression.path)
row.names(expression) = expression[,1]
expression = expression[,-1]
sample = colnames(expression)
expression = apply(expression,c(1,2),cutexpression)

if(l+z == 0){
  stop("please choose at least one mode to scale data, --log, --zscore")
}


#load expression matrix
expression = read.csv(expression.path)
row.names(expression) = expression[,1]
expression = expression[,-1]
sample = colnames(expression)
#z-score
if(z == T){
  expression.z = as.data.frame(t(apply(expression,1,scale)))
  colnames(expression) = sample}

if(l==T){
  expression.l = log2(expression+1)
}




#load path
path = read.table(poi,sep = ",", colClasses = "character",header = T)


#plot
setwd(OUTDIR)
for(i in 1:nrow(path)){
  id = path$id[i]
  name = path$name[i]
  if(z == T){
    #for z-scored expression, if a node(rectangle in the output graph) is mapped to multiple genes, the gene with max absulote expression 
    p.1 <- pathview(gene.data =  as.matrix(expression.z),pathway.id = id, 
                    species = "eco",gene.idtype = "symbol", 
                    kegg.native = T,same.layer=F,out.suffix = paste("zscore",name),
                    limit = list(gene=3),bins = list(gene = 30),
                    high = "#710320",mid = "#FFFFFF",low = "#072F63",
                    new.signature = F,res = 600,cex = 0.12,plot.col.key = F,node.sum = "max.abs")
    }
  if(l == T){
    #for log2 scaled data, the average exxpression is used to summary multi mapping nodes
    p.1 <- pathview(gene.data =  as.matrix(expression.l),pathway.id = id, 
                    species = "eco",gene.idtype = "symbol", 
                    kegg.native = T,same.layer=F,out.suffix = paste("log2",name),
                    limit = list(gene=c(0,13)),bins = list(gene = 30),
                    high = "#641E16",mid = "#CD6155",low = "#FFFFFF",
                    new.signature = F,res = 600,cex = 0.12,plot.col.key = T,node.sum = "mean")
  }
}




