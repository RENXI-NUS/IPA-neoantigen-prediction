library(dplyr)
library(caroline)

###
# Rscript run.r for expression table of IPA RI and SM
###
my_args = commandArgs(trailingOnly = TRUE)
read_length = my_args[3] %>% as.numeric()
sample = my_args[2]
outdir = my_args[1]
datadir = paste0(paste(outdir, sample, sep = "/"), ".")


### prepare the input
#cat("reading input files ...\n")
cluster = read.table(paste(datadir,"intron.clusterID.output",sep = ""))
upsacount = read.table(paste(datadir,"softclip.polya.ups100.SAF.all.count",sep = ""),head = T)
upsjcount = read.table(paste(datadir,"softclip.polya.ups100.SAF.junc.count",sep = ""),head = T)
dnsacount = read.table(paste(datadir,"softclip.polya.dns100.SAF.all.count",sep = ""),head = T)
dnsjcount = read.table(paste(datadir,"softclip.polya.dns100.SAF.junc.count",sep = ""),head = T)
totalcount = read.table(paste(datadir,"total.aligned.all.count",sep = ""))
subcluster = cluster[cluster$V4 %in% upsacount$Geneid,]
ss5inf = read.table(paste(datadir,"softclip.polya.ups24.ss5.wo",sep = ""))
ss3inf = read.table(paste(datadir,"softclip.polya.dns24.ss3.wo",sep = ""))
### checking
if(sum(as.character(subcluster[,4]) == as.character(upsacount[,1])) == nrow(subcluster)){
    cat("there are",nrow(subcluster),"PAS cluster ...\n")
}else{
    cat("error:cluster ID is not consistent across different input files\n")
}
# upstream and downstram coverage
tmp1 = rep(0,nrow(upsjcount))
tmp1[upsjcount[,1] %in% ss5inf[,4]] = upsjcount[upsjcount[,1] %in% ss5inf[,4],7]
upsrpm = (upsacount[,7] - tmp1)/unlist(totalcount[1]/1000000)
tmp2 = rep(0,nrow(dnsjcount))
tmp2[dnsjcount[,1] %in% ss3inf[,4]] = dnsjcount[dnsjcount[,1] %in% ss3inf[,4],7]
dnsrpm = (dnsacount[,7] - tmp2)/unlist(totalcount[1]/1000000)

# merge to one table
IPA_table = cbind(upsrpm,dnsrpm,subcluster$V5) 
colnames(IPA_table) = c("upsrpm","dnsrpm","soft_clip_reads")
IPA_table = IPA_table %>% as.data.frame() %>% mutate(geneid = as.character(subcluster$V4)) %>% mutate(soft_clipped_rpm=(soft_clip_reads/(read_length-3))*100/unlist(totalcount[1]/1000000)) %>% select(geneid,upsrpm,dnsrpm,soft_clipped_rpm) %>% mutate(IPA_transcript_expression = pmax(soft_clipped_rpm,(upsrpm-dnsrpm)))
IPA_table = IPA_table[IPA_table$IPA_transcript_expression>0.1,]
colnames(cluster) = c("chr", "start", "end", "geneid", "soft_clip_counts", "strand")
IPA_merged = merge(cluster, IPA_table) %>% select(chr,start,end,geneid,IPA_transcript_expression,strand)
write.delim(IPA_merged, paste0(datadir,"ipa.table"), col.names = T)

