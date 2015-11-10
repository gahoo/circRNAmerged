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
      nintersects=100,
      order.by='freq',
      number.angles =45,
      decreasing=c(T,F)
      )


ciri_df$count<-rowSums(ciri_df[,2:11])

ciri_df$count %>%
  table %>%
  data.frame %>%
  ggplot(aes(x=., y=Freq)) +
  geom_bar(stat='identity')


ciri_rbind %>%
  select(circRNA_ID,sample,junction_reads_ratio) %>%
  spread(sample, junction_reads_ratio) ->
  ciri_ratio
ciri_ratio$sd<-apply(100*ciri_ratio[,2:11], 1 , sd, na.rm=T)

ciri_rbind_count_ratio_sd<-merge(ciri_rbind_count, ciri_ratio[,c('circRNA_ID','sd')])
