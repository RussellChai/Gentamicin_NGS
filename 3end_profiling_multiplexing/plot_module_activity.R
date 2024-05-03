#ChaiRuochen
#20210506
#This script is to save heatmap and bubble plot of regulator activities
#The output plots are [0,1] normalized or show as non-normalized data.
#To execute:
#Rscript plot_activity.R --sigma_activity <sigma.activity.file.xlsx> --TF_activity <TF.activity.file.xlsx> --output <output directory>


library(tidyverse, quietly = T)
library(argparser, quietly = T)
library(pheatmap, quietly = T)
library(Rmisc, quietly = T)

p <- arg_parser("plot activity")
#add command line arguments
p <- add_argument(p,"--sigma_activity",help ="sigma factor activity file, derived from python script",
                  type="character",default = NULL)
p <- add_argument(p,"--TF_activity",help ="TF activity file, derived from python script",
                  type="character",default = NULL)
p <- add_argument(p,"--output",help ="ouput directory", 
                  type="character",default = '.')

argv <- parse_args(p)
sigma.activity.file <- argv$sigma_activity
TF.activity.file <- argv$TF_activity
OUTDIR <- argv$output


######define minor functions

normalize01 = function(x){
  y = (x-min(x))/(max(x)-min(x))
  return(y)
}

plot_sigma = function(sigma){
  sigma.level <- c("Sigma70","Sigma54","Sigma28","Sigma38","Sigma32","Sigma24")
  
  df.sigma <- left_join(sigma %>% pivot_longer(cols = starts_with("lib"),
                                               names_to = "sample",
                                               values_to = "activity"),
                        sigma.quantile %>% pivot_longer(cols = starts_with("lib"),
                                                        names_to = "sample",
                                                        values_to = "quantile"),
                        by = c("sigma_factor","sample")) %>% 
    mutate(sigma_factor = factor(sigma_factor,levels = sigma.level))
  
  
  
  p.b <- ggplot(df.sigma)+
    geom_point(aes(x = sample,y = sigma_factor,color = activity,size = quantile))+
    coord_flip()+
    theme(panel.background=element_rect(fill='transparent', color='gray'),
          axis.text.x = element_text(hjust = 1,angle = 65))+
    scale_color_gradient(low = "#FFF5F0",high = "#670005")+
    theme(axis.title.x =element_text(size=12), 
          axis.title.y=element_text(size=12),
          axis.text.x =element_text(size=14), 
          axis.text.y=element_text(size=14))
  
  
  sigma <- as.data.frame(sigma)
  row.names(sigma) = sigma$sigma_factor
  sigma = sigma[sigma.level,-1]
  sigma = sigma[,ncol(sigma):1]
  p.h <- pheatmap(t(sigma),cluster_rows = F,cluster_cols = F,
                  color = colorRampPalette(c("#FFF5F0", "#670005"))(50),
                  fontsize = 10.5,angle_col = 45,border_color = "black",
                  cellwidth = 25,cellheight = 25)
  plot = list(bubble = p.b,heatmap=p.h)
  return(plot)
}

plot_tf = function(tf){
  c2 = nrow(tf)
  c1 = round(c2/2)
  p.b.1 <- left_join(tf[1:c1,] %>% pivot_longer(cols = starts_with("lib"),
                                         names_to = "sample",
                                         values_to = "activity"),
                     tf.quantile %>% pivot_longer(cols = starts_with("lib"),
                                                  names_to = "sample",
                                                  values_to = "quantile"),
                     by = c("TF_name","sample")) %>% 
    #mutate(sample = sample)
    ggplot()+
    geom_point(aes(x = sample,y = TF_name,color = activity,size = quantile))+
    coord_flip()+
    theme(panel.background=element_rect(fill='transparent', color='gray'),
          axis.text.x = element_text(hjust = 1,angle = 65))+
    scale_color_gradient(low = "#FFF5F0",high = "#670005")+
    theme(axis.title.x =element_text(size=12), 
          axis.title.y=element_text(size=12),
          axis.text.x =element_text(size=14), 
          axis.text.y=element_text(size=14))
  
  p.b.2 <- left_join(tf[(c1+1):c2,] %>% pivot_longer(cols = starts_with("lib"),
                                                names_to = "sample",
                                                values_to = "activity"),
                     tf.quantile %>% pivot_longer(cols = starts_with("lib"),
                                                  names_to = "sample",
                                                  values_to = "quantile"),
                     by = c("TF_name","sample")) %>% 
    #mutate(sample = sample)
    ggplot()+
    geom_point(aes(x = sample,y = TF_name,color = activity,size = quantile))+
    coord_flip()+
    theme(panel.background=element_rect(fill='transparent', color='gray'),
          axis.text.x = element_text(hjust = 1,angle = 65))+
    scale_color_gradient(low = "#FFF5F0",high = "#670005")+
    theme(axis.title.x =element_text(size=12), 
          axis.title.y=element_text(size=12),
          axis.text.x =element_text(size=14), 
          axis.text.y=element_text(size=14))
  
  
  tf<- as.data.frame(tf)
  row.names(tf) = tf$TF_name
  tf = tf[,-1]
  tf = tf[,ncol(tf):1]
  c2 = nrow(tf)
  c1 = round(c2/2)
  p.h.1 <- pheatmap(t(tf) %>% .[,1:c1],cluster_rows = F,cluster_cols = F,
                  color = colorRampPalette(c("#FFF5F0", "#670005"))(50),
                  fontsize = 10.5,angle_col = 45,border_color = "black",
                  cellwidth = 25,cellheight = 25)
  p.h.2 <- pheatmap(t(tf) %>% .[,(c1+1):c2],cluster_rows = F,cluster_cols = F,
                    color = colorRampPalette(c("#FFF5F0", "#670005"))(50),
                    fontsize = 10.5,angle_col = 45,border_color = "black",
                    cellwidth = 25,cellheight = 25)
  plot = list(bubble.1 = p.b.1,bubble.2=p.b.2,heatmap.1=p.h.1,heatmap.2=p.h.2)
  return(plot)
}

####bubble plot & heatmap of sigma activity
#####read files
#
#sigma.activity.file = "sigma_activity_split.xlsx"
  sigma.activity <- readxl::read_excel(sigma.activity.file,
                                       sheet = "sigma_activity_mean")
  sigma.activity.norm <- as.data.frame(t(apply(sigma.activity[,-1],1,normalize01)))
  colnames(sigma.activity.norm) <- colnames(sigma.activity)[-1]
  sigma.activity.norm$sigma_factor <- sigma.activity$sigma_factor
  sigma.activity.norm <- select(sigma.activity.norm,sigma_factor,!contains("sigma_factor"))
  sigma.quantile <- readxl::read_excel(sigma.activity.file,
                                       sheet = "quantile_mean")
  plot = plot_sigma(sigma = sigma.activity)
  
  pdf(file = paste(OUTDIR,"sigma_activity_bubble.pdf",sep = "/"),width = 4, height =4)
  plot$bubble
  dev.off()
  pdf(file = paste(OUTDIR,"sigma_activity_heatmap.pdf",sep = "/"),width = 4, height =4)
  plot$heatmap
  dev.off()
  
  plot = plot_sigma(sigma = sigma.activity.norm)
  
  pdf(file = paste(OUTDIR,"sigma_activity_norm_bubble.pdf",sep = "/"),width = 4, height =4)
  plot$bubble
  dev.off()
  pdf(file = paste(OUTDIR,"sigma_activity_norm_heatmap.pdf",sep = "/"),width = 4, height =4)
  plot$heatmap
  dev.off()



####bubble plot & heatmap of TF activity
#####read files

  tf.quantile <- readxl::read_excel(TF.activity.file,
                                    sheet = "quantile_mean")
  tf.activity <- readxl::read_excel(TF.activity.file,
                                    sheet = "tf_activity_mean")
  tf.activity.norm <- as.data.frame(t(apply(tf.activity[,-1],1,normalize01)))
  colnames(tf.activity.norm) <- colnames(tf.activity)[-1]
  tf.activity.norm$TF_name <- tf.activity$TF_name
  tf.activity.norm <- select(tf.activity.norm,TF_name,!contains("TF_name"))
  
  plot = plot_tf(tf = tf.activity)
  
  pdf(file = paste(OUTDIR,"TF_activity_bubble.pdf",sep = "/"),width = 40, height =8)
  multiplot(plot$bubble.1,plot$bubble.2)
  dev.off()
  
  pdf(file = paste("TF_activity_heatmap_1.pdf",sep = "/"),width = 60, height =10)
  plot$heatmap.1
  dev.off()
  pdf(file = paste("TF_activity_heatmap_2.pdf",sep = "/"),width = 60, height =10)
  plot$heatmap.2
  dev.off()
  
  plot = plot_tf(tf = tf.activity.norm)
  
  pdf(file = paste(OUTDIR,"TF_activity_norm_bubble.pdf",sep = "/"),width = 40, height =8)
  multiplot(plot$bubble.1,plot$bubble.2)
  dev.off()
  
  pdf(file = paste("TF_activity_norm_heatmap_1.pdf",sep = "/"),width = 60, height =10)
  plot$heatmap.1
  dev.off()
  pdf(file = paste("TF_activity_norm_heatmap_2.pdf",sep = "/"),width = 60, height =10)
  plot$heatmap.2
  dev.off()




