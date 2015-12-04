library(dplyr)
library(GenomicRanges)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

GeneID_SYMBOL<-AnnotationDbi::select(org.Hs.eg.db, 
                      keys = keys(org.Hs.eg.db),
                      columns = c("ENTREZID", "SYMBOL")) %>%
  rename_(GENEID="ENTREZID")

AnnotationDbi::select(TxDb.Hsapiens.UCSC.hg19.knownGene,
       keys = keys(TxDb.Hsapiens.UCSC.hg19.knownGene, keytype='GENEID'),
       columns = c('GENEID', 'TXID', 'TXCHROM', 'TXSTRAND', 'TXSTART', 'TXEND'),
       keytype = 'GENEID') %>%
  group_by(GENEID, TXCHROM, TXSTRAND) %>%
  summarise(start=min(TXSTART),
            end = max(TXEND))%>%
  merge(GeneID_SYMBOL, by='GENEID', all.x=T) %>%
  makeGRangesFromDataFrame(seqnames.field="TXCHROM",
                           strand.field="TXSTRAND",
                           keep.extra.columns=T) ->
  GeneRanges

ciri_rbind %>%
  #filter(circRNA_type == 'intergenic region') %>%
  filter(!is.na(circRNA_type)) %>%
  dplyr::select(chr, circRNA_start, circRNA_end, circRNA_ID) %>%
  unique %>%
  makeGRangesFromDataFrame(seqnames.field='chr',
                           start.field='circRNA_start',
                           end.field='circRNA_end',
                           keep.extra.columns=T) ->
  intergenicRanges

hits<-findOverlaps(intergenicRanges, GeneRanges)
anno <- cbind(mcols(intergenicRanges[queryHits(hits)]),
      mcols(GeneRanges[subjectHits(hits)])) %>%
  unique %>%
  as.data.frame %>%
  group_by(circRNA_ID) %>%
  summarise(region_gene_id=paste0(GENEID, collapse = ','),
            region_symbol=paste0(SYMBOL, collapse = ','),
            region_gene_cnt=n())

ciri_rbind %>%
  mutate(depth.total = junction.Normal+junction.Tumor+non_junction.Normal+non_junction.Tumor,
         depth.Normal = junction.Normal+non_junction.Normal,
         depth.Tumor = junction.Tumor+non_junction.Tumor) %>%
  left_join(anno) ->
  ciri_rbind_anno

abs_diff<-function(x){sqrt(sum(x^2)/length(x))}

write.table(ciri_rbind_anno, file='ciri_rbind_region_check.txt', row.names=F, quote=F, sep='\t')

ciri_rbind_anno %>%
  dplyr::select(circRNA_ID, ratio.Normal, ratio.Tumor) %>%
  mutate(ratio.Normal = 100 * ratio.Normal,
         ratio.Tumor = 100 * ratio.Tumor,
         ratio.diff = ratio.Tumor - ratio.Normal) %>%
  group_by(circRNA_ID) %>%
  summarise(occurrence = n(),
            ratio.Normal.sd = sd(ratio.Normal, na.rm=T),
            ratio.Tumor.sd = sd(ratio.Tumor, na.rm=T),
            ratio.Diff.sd = sd(ratio.diff, na.rm=T),
            #ratio.Normal.sd = ifelse(is.na(ratio.Normal.sd),1,ratio.Normal.sd),
            #ratio.Tumor.sd = ifelse(is.na(ratio.Tumor.sd),1,ratio.Tumor.sd),
            #ratio.Diff.sd = ifelse(is.na(ratio.Diff.sd),1,ratio.Diff.sd),
            ratio.abs_diff = abs_diff(ratio.diff),
            ratio.rank = ratio.abs_diff/(ratio.Normal.sd * ratio.Tumor.sd * ratio.Diff.sd),
            ratio.rank2 = ratio.abs_diff/(ratio.Normal.sd * ratio.Tumor.sd),
            ratio.rank3 = ratio.abs_diff/ratio.Diff.sd
  ) ->
  ciri_rbind_anno_rank

write.table(ciri_rbind_anno_rank, file='ciri_rbind_rank.txt', row.names=F, quote=F, sep='\t')

countSample<-function(df){
  df %>%
  dplyr::select(circRNA_ID) %>%
    group_by(circRNA_ID) %>%
    summarise(sample_cnt=n())
}

ciri_rbind_anno %>%
  filter(junction.Normal <= 1,
         junction.Tumor >= 3,
         depth.total >= 10,
         p.values <=0.05) ->
  ciri_rbind_anno_more_tumor

  
ciri_rbind_anno %>%
  filter(junction.Normal >= 3,
         junction.Tumor <= 1,
         depth.total >= 10,
         p.values <=0.05) ->
  ciri_rbind_anno_more_normal

ciri_rbind_anno_more_normal<-merge(ciri_rbind_anno_more_normal,
                                   countSample(ciri_rbind_anno_more_normal),
                                   all=T)
ciri_rbind_anno_more_tumor<-merge(ciri_rbind_anno_more_tumor,
                                  countSample(ciri_rbind_anno_more_tumor),
                                  all=T)

write.table(ciri_rbind_anno_more_normal, file='ciri_rbind_region_check_more_normal.txt', row.names=F, quote=F, sep='\t')
write.table(ciri_rbind_anno_more_tumor, file='ciri_rbind_region_check_more_tumor.txt', row.names=F, quote=F, sep='\t')

getSetDf<-function(df){
  df %>%
    select(circRNA_ID, sample, p.values) %>%
    mutate(p.values = 1) %>%
    spread(sample, p.values, fill=0)
}

plotSampleSets<-function(df, ...){
  setCnt<-length(unique(df$sample))
   df %>%
    getSetDf %>%
    upset(nsets=setCnt, 
        #nintersects=50,
        order.by='freq',
        #number.angles =45,
        decreasing=c(T,F),
        ...
        )
}

filterCnt<-function(df, cnt){
  df %>%
    filter(sample_cnt >= cnt) %>%
    "$"('circRNA_ID') %>%
    as.character %>%
    unique
}

plotCircSets<-function(df, cnt=0, ...){
  candidates<-filterCnt(df, cnt)
  df_set<-getSetDf(df)
  df_set_t<-df_set
  df_set %>%
    left_join(unique(df[c('circRNA_ID', 'symbol', 'region_symbol')]), by='circRNA_ID') %>%
    mutate(anno = sprintf("%s.%s", 
                          ifelse(is.na(symbol), region_symbol, as.character(symbol)),
                          circRNA_ID)) -> df_set
  row.names(df_set_t)<-unique(df_set$anno)
  df_set %>%
    filter(circRNA_ID %in% candidates) %>%
    "$"("anno") -> candidates
  #row.names(df_set_t)<-df_set$circRNA_ID
  df_set_t<-as.data.frame(t(df_set_t[,-1]))
  upset(df_set_t, sets = candidates, nsets = length(candidates), ...)
}

# plotCircSets<-function(df, cnt=0, ...){
#   candidates<-filterCnt(df, cnt)
#   df_set<-getSetDf(df)
#   df_set_t<-df_set
#   row.names(df_set_t)<-df_set$circRNA_ID
#   df_set_t<-as.data.frame(t(df_set_t[,-1]))
#   upset(df_set_t, sets = candidates, nsets = length(candidates), ...)
# }

pdf('plots/ciri_rbind_anno_more_normal.samples.pdf', width=8, height=6)
plotSampleSets(ciri_rbind_anno_more_normal,nintersects=20)
dev.off()
pdf('plots/ciri_rbind_anno_more_tumor.samples.pdf', width=8, height=6)
plotSampleSets(ciri_rbind_anno_more_tumor)
dev.off()
pdf('plots/ciri_rbind_anno_more_tumor.circRNA_Sets.inCGC.sampleCnt_ge2.pdf', width=8, height=6)
plotCircSets(ciri_rbind_anno_more_tumor %>% filter(inCGC==T),2)
dev.off()
pdf('plots/ciri_rbind_anno_more_tumor.circRNA_Sets.ge3.pdf', width=8, height=6)
plotCircSets(ciri_rbind_anno_more_tumor, 3)
dev.off()
pdf('plots/ciri_rbind_anno_more_normal.circRNA_Sets.inCGC.sampleCnt_ge3.pdf', width=8, height=6)
plotCircSets(ciri_rbind_anno_more_normal %>% filter(inCGC==T), 3)
dev.off()
pdf('plots/ciri_rbind_anno_more_normal.circRNA_Sets.sampleCnt_ge6.pdf', width=8, height=6)
plotCircSets(ciri_rbind_anno_more_normal, 6)
dev.off()

plotRelExpPattern<-function(df, circRNA_IDs){
  df %>%
    filter(circRNA_ID %in% circRNA_IDs) %>%
    dplyr::select(ratio.Normal, ratio.Tumor, sample, circRNA_ID, p.values, symbol, region_symbol) %>%
    gather(type,ratio,ratio.Normal,ratio.Tumor) %>%
    mutate(type=gsub('ratio.','',type),
           significant=p.values<=0.05,
           log10P=log10(p.values),
           anno = sprintf("%s\n%s", 
                          ifelse(is.na(symbol), region_symbol, as.character(symbol)),
                          circRNA_ID)
           ) %>%
    ggplot(aes(x=sample, y=ratio, group=type, color=type)) + 
    geom_line() +
    geom_point(aes(size=-log10P, alpha=significant)) +
    #scale_size_continuous(range=c(2,8)) +
    facet_grid(anno~.) +
    ylab('Relative Expression Ratio') +
    theme(axis.text.x=element_text(angle=90))
}

pdf('plots/relExpPattern.pdf', height=12, width=12)
plotRelExpPattern(ciri_rbind_anno, 
                  circRNA_IDs=filterCnt(ciri_rbind_anno_more_tumor %>% filter(inCGC==T), 2) ) +
                  ggtitle("Relative Expression Pattern of\n Tumor circRNA in CGC and more than 1 pairs")
plotRelExpPattern(ciri_rbind_anno, 
                  circRNA_IDs=filterCnt(ciri_rbind_anno_more_tumor, 3) ) +
                  ggtitle("Relative Expression Pattern of\n Tumor circRNA in more than 2 pairs")
plotRelExpPattern(ciri_rbind_anno, 
                  circRNA_IDs=filterCnt(ciri_rbind_anno_more_normal %>% filter(inCGC==T), 3) ) +
                  ggtitle("Relative Expression Pattern of\n Normal circRNA in CGC and more than 2 pairs")
plotRelExpPattern(ciri_rbind_anno, 
                  circRNA_IDs=filterCnt(ciri_rbind_anno_more_normal, 6) ) +
                  ggtitle("Relative Expression Pattern of\n Normal circRNA in more than 5 pairs")
dev.off()