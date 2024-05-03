# This script is to merge depth files, output enrichment file and 3' end freq.
# Usage: Rscript end_enrichment_from_bam.R -b <bam file path> -p <forward depth> -n <reverse depth> -o <output directory> 
# Note that this script is dedicated for library structure same as 3'end-seq of DR-seq.
# Written by Ruochen Chai @ Tsinghua
# Last updated, Dec 01, 2022

library(jsonlite)
library(argparser, quietly = T)
library(tidyverse)
library(DescTools)
library(parallel)
library(aplot)
library(gggenes)
library(Rsamtools,quietly=TRUE)
library(reshape2,quietly=TRUE)
library(snowfall,quietly=TRUE)
library(APERO,quietly=TRUE)
library(data.table)


p <- arg_parser("end_enrichment_from_bam")
#add command line arguments
p <- add_argument(p,"--bam_path",help ="path of input bam file",
                  type="character",default = NULL)
#p <- add_argument(p,"--invert_strand",help ="invert fragments strand",
#                  flag = TRUE)
p <- add_argument(p,"--out_path",help ="output directory",
                  type="character",default = NULL)
p <- add_argument(p,"--pos_depth",help ="positive depth file",
                  type="character",default = NULL)
p <- add_argument(p,"--neg_depth",help ="negative depth file",
                  type="character",default = NULL)

argv <- parse_args(p)
my_bam_path <- argv$bam_path
output_dir <- argv$out_path
invertstr <- TRUE
forward_depth <- argv$pos_depth
reverse_depth <- argv$neg_depth



#generate coverage file
merge_depth_file = function(pos_path,neg_path,out_path,invert_str,total_reads){
  if(invert_str == T){
    pos_depth = read_tsv(neg_path,
                         col_names = c("chrom","position","depth")) %>% 
      mutate(str = "+")
    
    neg_depth = read_tsv(pos_path,
                         col_names = c("chrom","position","depth")) %>% 
      mutate(str = "-")
    
    write.csv(rbind(pos_depth,neg_depth) %>% 
                select(!chrom) %>% 
                mutate(depth = depth*1000000/total_reads),
              out_path,row.names = F)
  }else{
    pos_depth = read_tsv(pos_path,
                         col_names = c("chrom","position","depth")) %>% 
      mutate(str = "+")
    
    neg_depth = read_tsv(neg_path,
                         col_names = c("chrom","position","depth")) %>% 
      mutate(str = "-")
    
    write.csv(rbind(pos_depth,neg_depth) %>% 
                select(!chrom) %>% 
                mutate(depth = depth*1000000/total_reads),
              out_path,row.names = F)
  }

}
enrich_end_parallel = function(end3_freq,w = 1,f = 50,genome_size){
  colnames(end3_freq)[1] = "end"
  
  end3_freq$sum_freq = as.numeric(as.character(end3_freq$freq))
  dp = end3_freq[which(end3_freq$str == "+"), ]
  dn = end3_freq[which(end3_freq$str == "-"), ]
  M = data.frame(P = rep(0, genome_size), N = rep(0, genome_size))
  M$P[dp$end] = dp$freq
  M$N[dn$end] = dn$freq
  
  
  enrichment = data.frame(P=c(),N=c())
  
  genome_position = 1:genome_size
  
  get_enrich = function(g,w,f,e){
    window = c((g-w):(g+w)) %>% genome_circulazition(genome_size)
    window_freq = sum(e[window])
    flank = c((g-w-f):(g-w-1),(g+w+1):(g+w+f)) %>% genome_circulazition(genome_size)
    flank_freq = sum(e[flank])
    flank_freq = ifelse(flank_freq==0,1,flank_freq)
    enrich = (window_freq/(2*w+1))/(flank_freq/(2*f))
    
    return(c(window_freq,flank_freq,enrich))
  }
  
  
  enrichment_P = lapply(genome_position,get_enrich,w,f,M$P)
  enrichment_P = data.frame(t(matrix(unlist(enrichment_P),ncol = genome_size)))
  colnames(enrichment_P) = c("window_freq","flank_freq","enrichment")
  enrichment_P$str = "+"
  enrichment_P$position = c(1:genome_size)
  enrichment_P$end_freq = M$P
  enrichment_N = lapply(genome_position,get_enrich,w,f,M$N)
  enrichment_N = data.frame(t(matrix(unlist(enrichment_N),ncol = genome_size)))
  colnames(enrichment_N) = c("window_freq","flank_freq","enrichment")
  enrichment_N$str = "-"
  enrichment_N$position = c(1:genome_size)
  enrichment_N$end_freq = M$N
  enrichment = rbindlist(list(enrichment_P,enrichment_N))
  enrichment$w = w
  enrichment$f = f
  
  enrichment = enrichment %>% 
    mutate(enrichment = ifelse(is.nan(enrichment)==T,
                               0,enrichment))
  
  #enrichment = unlist(enrichment)
  return(enrichment)
}

draw_ecdf_enrichment = function(enrichment){
  ecdf_enrichment = ecdf(enrichment$enrichment)
  p = enrichment %>% 
    mutate(cdf = ecdf_enrichment(enrichment)) %>% 
    ggplot(aes(x = enrichment,y = 1-cdf))+
    geom_line(color = "#922B21",size = 1.2)+
    scale_y_log10(limits = c(10^-6,10^0),
                  breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x)))+
    theme_bw()+
    theme(plot.margin = unit(rep(3,4),"mm"))+
    scale_x_continuous(expand = c(0,0),limits = c(0,1000))
  
  return(p)
}

draw_ecdf_end_CPM = function(enrichment){
  ecdf_endCPM = ecdf(enrichment$end_CPM)
  p = enrichment %>% 
    mutate(cdf = ecdf_endCPM(end_CPM)) %>% 
    ggplot(aes(x = end_CPM,y = 1-cdf))+
    geom_line(color = "#922B21",size = 1.2)+
    scale_y_log10()+
    theme_bw()+
    theme(plot.margin = unit(rep(3,4),"mm"))+
    scale_x_continuous(expand = c(0,0),limits = c(0,1000))
  
  return(p)
}

get_fragment_from_pe_bam = function(PE_bam_file){
  a=scanBam(PE_bam_file,
            param=ScanBamParam(what=c("pos","flag","mpos","qname","strand","cigar"), 
                               flag=scanBamFlag(isUnmappedQuery=FALSE, 
                                                hasUnmappedMate=FALSE)))
  pos=a[[1]]$pos
  flag=a[[1]]$flag
  mpos=a[[1]]$mpos
  qname=a[[1]]$qname
  str=a[[1]]$strand
  cigar=a[[1]]$cigar
  bam=data.frame(pos=pos,mpos=mpos,qname=qname,str=str,flag=flag,cigar=cigar)
  rm(a,pos,flag,mpos,qname,str,cigar)
  
  FLAG=c(99,147,83,163)
  bam=bam[which(bam$flag%in%FLAG),]
  rm(FLAG)
  
  bam$qwidth=apply(bam,1,APERO::qw)
  
  ap = bam[which(bam$str == "+"), ]
  an = bam[which(bam$str == "-"), ]
  #rm(bam)
  ap$right_end_read = ap$pos + ap$qwidth - 1
  an$right_end_read = an$pos + an$qwidth - 1
  p = match(ap$qname, an$qname, nomatch = 0)
  ap$right_end_mate = an$right_end_read[p]
  n = match(an$qname, ap$qname, nomatch = 0)
  an$right_end_mate = ap$right_end_read[n]
  r1 = rbind(ap, an)
  rm(ap, an, p, n)
  r1p = r1[which(r1$str == "+"), ]
  r1n = r1[which(r1$str == "-"), ]
  r1p$lg = r1p$right_end_mate - r1p$pos + 1
  r1n$lg = r1n$right_end_read - r1n$mpos + 1
  colnames(r1p) = c("5end", "mate3end", "qname", "str", "flag", 
                    "cigar", "qwidth", "3end", "mate5end", "TLEN")
  colnames(r1n) = c("3end", "mate5end", "qname", "str", "flag", 
                    "cigar", "qwidth", "5end", "mate3end", "TLEN")
  r1n = r1n[, c(8, 9, 3, 4, 5, 6, 7, 1, 2, 10)]
  r1 = rbind(r1p, r1n)
  r1 = r1[which(r1$flag < 128), ]
  rm(r1p, r1n)
  r1 = r1 %>% 
    mutate(left = ifelse(str == "+",`5end`,mate5end),right = ifelse(str == "+",mate5end,`5end`)) %>% 
    select(qname,str,left,right)
  
  return(r1)
} 

get_3end_freq_from_fragments = function(input_fragment){
  
  input_fragment = input_fragment %>% mutate("3end" = ifelse(str == "+",right,left))
  
  input_fragment = input_fragment[,c("3end","str")]
  input_fragment$agreg=paste(input_fragment$'3end',input_fragment$str,sep=".")
  
  c = data.frame(table(input_fragment$agreg))
  c = data.frame(colsplit(c$Var1,"\\.",c("3end","str")),sum_freq=c$Freq)
  c=c[,c(1,3,2)]
  colnames(c)=c("end","freq","str")
  
  c= c %>% arrange(end)
  
  return(c)
  
}

genome_circulazition = function(x,gs){
  x = ifelse(x<=0,gs+x,x)
  x = ifelse(x>gs,x-gs,x)
}
#merge depth file

e_coli_genome_size = 4641652
my_fragments = get_fragment_from_pe_bam(my_bam_path)
tr = dim(my_fragments)[1]

merge_depth_file(pos_path=forward_depth,
                 neg_path=reverse_depth,
                 out_path=str_c(output_dir,"/coverage.csv"),
                 invert_str=invertstr,
                 total_reads=tr)


#count functional reads and all reads
stop_bed = './genomeFile/gene_stop.bed'
stop_out = str_c(output_dir,'/functional_count.txt')
count_cmd = str_c('/usr/bin/bedtools2/bin/bedtools multicov -bams ',my_bam_path,' -bed ',stop_bed, ' -s > ',stop_out)
system(count_cmd)
all_bed = './genomeFile/gene_extension.bed'
all_out = str_c(output_dir,'/all_count.txt')
count_cmd = str_c('/usr/bin/bedtools2/bin/bedtools multicov -bams ',my_bam_path,' -bed ',all_bed, ' -s > ',all_out)
system(count_cmd)
all_reads=read.csv(all_out,sep='\t',header = F)
colnames(all_reads)[4]='gene'
colnames(all_reads)[7]='all'
functional_reads=read.csv(stop_out,sep='\t',header = F)
colnames(functional_reads)[4]='gene'
colnames(functional_reads)[7]='functional'
reads = left_join(functional_reads %>% select(gene,functional),
             all_reads %>% select(gene,all),by='gene')
file.remove(all_out)
file.remove(stop_out)                  
write.csv(reads,str_c(output_dir,'functional_count.csv'),row.names=F)
print('functionally counted')


if(invertstr == T){
  my_fragments$str = ifelse(my_fragments$str == "+","-","+")
}

my_fragments = my_fragments %>% mutate(end3 = ifelse(str == "+",right,left),
                                       end5 = ifelse(str == "+",left,right))


my_end_freq = get_3end_freq_from_fragments(my_fragments)
my_enrichment = enrich_end_parallel(my_end_freq,w = 1,f = 50,genome_size = e_coli_genome_size)
my_enrichment$end_CPM = my_enrichment$end_freq*1000000/sum(my_enrichment$end_freq)


my_list = my_fragments %>% mutate(position = end3) %>% 
  left_join(my_enrichment) %>% 
  mutate(other_end = end5) %>% 
  select(position,str,end_CPM,enrichment,other_end)%>% 
  dplyr::group_by(str,position,end_CPM,enrichment) %>%
  dplyr::summarise(other_end = str_c(other_end,collapse = ",")) %>% 
  as.data.frame() %>% as.list() 
my_list$other_end = my_list$other_end %>% str_split(pattern = ",") %>% lapply(as.numeric)
 


my_list$w = 1
my_list$f = 50
json_output = str_c(output_dir,"/enrichment_reads_location.json")
cat(toJSON(my_list,pretty = T),file = json_output,labels = NULL,append = F,fill = F)

my_ecdf_enrichment = ecdf(my_enrichment$enrichment)
my_ecdf_df = my_enrichment %>% 
  mutate(cdf = my_ecdf_enrichment(enrichment)) %>% 
  select(enrichment,cdf)
ecdf_output = str_c(output_dir,"/cdf_enrichment.csv")
write.csv(my_ecdf_df,ecdf_output,row.names = F)

my_read_length = my_fragments %>% 
  mutate(length = abs(left-right)+1) %>% .$length %>% 
  table() %>% as.data.frame() %>%dplyr::rename("length"=".","count" = "Freq")
p=ggplot(my_read_length,aes(x=as.numeric(length),y=as.numeric(count)))+geom_col()+theme_bw()
pdf(str_c(output_dir,"/read_length_distribution.pdf"),height = 4,width = 4)
print(p)
dev.off()
length_output = str_c(output_dir,"/reads_length.csv")
write.csv(my_read_length,length_output,row.names = F)


p_ecdf = draw_ecdf_enrichment(my_enrichment)
pdf(str_c(output_dir,"/ecdf_enrichment.pdf"),height = 4,width = 4)
print(p_ecdf)
dev.off()



