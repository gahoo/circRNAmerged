library(UpSetR)
library(plyr)
library(dplyr)
library(tidyr)
library(ggplot2)

ciri_files<-dir('CIRI')
ciri_lst<-lapply(ciri_files, function(ciri_file){
  read.table(paste0('CIRI/', ciri_file), sep='\t',header=T) %>%
    mutate(sample=gsub('.CIRI.merged','',ciri_file)) ->
    data
  colnames(data)[11:16]<-c('junction.Normal', 'junction.Tumor',
                           'non_junction.Normal', 'non_junction.Tumor',
                           'ratio.Normal', 'ratio.Tumor')
  data
})

ciri<-lapply(ciri_lst, function(ci){ci %>% select(circRNA_ID, p.values)})

ciri_rbind<-do.call(rbind, ciri_lst)

ciri_rbind %>%
  group_by(circRNA_ID) %>%
  summarise(count=n()) %>%
  merge(ciri_rbind)->
  ciri_rbind_count

ciri_rbind_count %>%
  select(circRNA_ID, sample, p.values) %>%
  mutate(p.values = 1) %>%
  spread(sample, p.values, fill=0) ->
  ciri_df

# names(ciri)<-gsub('.CIRI.merged','',ciri_files)
# ciri_df<-ldply(ciri,.id='Somatic') %>%
#   mutate(p.values=1) %>%
#   #rename(Identifier = circRNA_ID) %>%
#   spread(Somatic, p.values, fill=0)
#   
# 
# caclIntersect<-function(n){
#   sum(sapply(1:n, function(i){choose(n, i)}))
# }

caclIntersect2<-function(n){
  2^n-1
}

upset(ciri_df, nsets=10, 
      #nintersects=caclIntersect2(length(ciri_df)-1),
      nintersects=50,
      order.by='freq',
      number.angles =45,
      decreasing=c(T,F)
      )

ciri_df_t<-ciri_df
row.names(ciri_df_t)<-ciri_df$circRNA_ID
ciri_df_t<-as.data.frame(t(ciri_df_t[,-1]))
upset(ciri_df_t, sets = colnames(ciri_df_t)[1:5])
stable_candidate<-c('chr10:103557737|103560157','chr10:103557737|103570071','chr10:103565802|103567658','chr10:105777918|105778666','chr10:112640991|112647644','chr10:112723883|112724819','chr10:116247717|116251637','chr10:120832402|120833449','chr10:121275021|121286936','chr10:12277047|12280484','chr10:123629521|123683844','chr10:123670421|123683844','chr10:123718839|123721032','chr10:126097111|126097534','chr10:128795012|128798571','chr10:128859932|128926028','chr10:14562976|14572514','chr10:22019829|22024164','chr10:27047991|27059274','chr10:27403450|27410377')
upset(ciri_df_t, sets = stable_candidate[1:5], nsets = 20, nintersects = 50)
idx<-colnames(ciri_df_t) %in% stable_candidate[1:5]
ciri_df_t[,idx]

ciri_df$count<-rowSums(ciri_df[,2:11])

ciri_df$count %>%
  table %>%
  data.frame %>%
  ggplot(aes(x=., y=Freq)) +
  geom_bar(stat='identity') +
  xlab('occurrence count')


ciri_rbind %>%
  select(circRNA_ID,sample,junction_reads_ratio) %>%
  spread(sample, junction_reads_ratio) ->
  ciri_ratio
ciri_ratio$sd<-apply(100*ciri_ratio[,2:11], 1 , sd, na.rm=T)

getColSD<-function(df, column, sd_name){
  df %>%
    select_("circRNA_ID", "sample", column) %>%
    spread_("sample", column) ->
    df_sd
  df_sd[[sd_name]]<-apply(100*df_sd[,2:ncol(df_sd)], 1 , sd, na.rm=T)
  df_sd[,c('circRNA_ID', sd_name)]
}

merge.all<-function(x, y){merge(x, y, all=T)}

lapply(c('junction_reads_ratio', 'ratio.Normal', 'ratio.Tumor', 'p.values'),
       function(column){
         getColSD(ciri_rbind, column, paste0(column,'.sd'))
       }) %>%
  Reduce(f=merge.all) ->
  ciri_sd

ciri_rbind_count %>%
  mutate(ratio.diff=ratio.Tumor-ratio.Normal) %>%
  left_join(ciri_sd, by='circRNA_ID') ->
  ciri_rbind_count_sd

write.table(ciri_rbind_count_sd, file='ciri_rbind_sd.txt', row.names=F, quote=F, sep='\t')

ciri_rbind_count_sd %>%
  filter(symbol=='ABL2', count>4, circRNA_ID != 'chr1:179087722|179091002') %>%
  select(ratio.Normal, ratio.Tumor, sample, circRNA_ID, p.values) %>%
  gather(type,ratio,ratio.Normal,ratio.Tumor) %>%
  mutate(type=gsub('ratio.','',type),
         logP=-log10(p.values)) %>%
  ggplot(aes(x=sample, y=ratio, group=type, color=type)) + 
  geom_line() +
  geom_point(aes(size=logP))+
  facet_grid(circRNA_ID~.) +
  ggtitle('ABL2') +
  theme(axis.text.x=element_text(angle=90))
