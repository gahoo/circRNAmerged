library(UpSetR)
library(GenomicRanges)
library(dplyr)

loadCIRI<-function(ciri_files){
  lapply(ciri_files, function(ciri_file){
    read.table(ciri_file, sep='\t',header=T,
               nrow=1000
               ) %>%
      mutate(sample=gsub('.*/','',gsub('.CIRI.merged','',ciri_file))) ->
      data
    colnames(data)[11:16]<-c('junction.Normal', 'junction.Tumor',
                             'non_junction.Normal', 'non_junction.Tumor',
                             'ratio.Normal', 'ratio.Tumor')
    data
  })
}

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



plotRelExpPattern<-function(df, circRNA_IDs){
  df %>%
    filter(circRNA_ID %in% circRNA_IDs) %>%
    select(ratio.Normal, ratio.Tumor, sample, circRNA_ID, p.values, symbol, region_symbol) %>%
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

getGeneArc<-function(df, symbols){
  prepare4rbind<-function(df, type){
    df %>% 
      select(circRNA_ID, chr, circRNA_start, circRNA_end,
             circRNA_type, p.values, sample, ends_with(type)) %>%
      rename_(junction = sprintf("junction.%s", type),
              non_junction = sprintf("non_junction.%s", type),
              relExp = sprintf("ratio.%s", type)
      ) %>%
      mutate(type=type,
             log10P=log10(p.values))
  }
  
  df %>%
    filter(symbol %in% symbols|region_symbol %in% symbols) %>%
    select(circRNA_ID, chr, circRNA_start, circRNA_end, circRNA_type,
                  junction.Normal, junction.Tumor, 
                  non_junction.Normal,non_junction.Tumor,
                  ratio.Normal, ratio.Tumor,p.values, sample) ->
    df
  rbind(
    df %>% prepare4rbind(type='Normal'),
    df %>% prepare4rbind(type='Tumor')
  )
}

plotArc<-function(arc){
  if(length(arc)==0){
    message("empty")
    NULL
  }else{
    message(length(arc))
    ggbio() +
      geom_arch(
        data=arc,
        ylab='Relative Expression Ratio',
        aes(fill = type, color = type, height = relExp, size = -log10P),
        alpha = 0.4,
        facets = sample ~ .) +
      scale_size(range = c(0, 6))
    
  }
}


plotTrack<-function(df, symbol, reads_cnt=0){
  df %>%
    getGeneArc(symbols=symbol) %>%
    filter(junction >= reads_cnt) %>%
    makeGRangesFromDataFrame(keep.extra.columns=T) %>%
    plotArc ->
    arc
  
  which_gene<-genesymbol[symbol]
  transcripts <- ggbio() + geom_alignment(data=txdb, which = which_gene )
  track_list<-list(
    arc=arc,
    transcripts=transcripts
  )
  
  track_list<-track_list[!sapply(track_list, is.null)]
  tracks(track_list, title=symbol, heights = c(3,1))
}

df2GRanges<-function(df){
  df %>%
    select(chr, circRNA_start, circRNA_end, circRNA_ID) %>%
    filter(!is.na(chr)) %>%
    unique %>%
    makeGRangesFromDataFrame(seqnames.field='chr',
                             start.field='circRNA_start',
                             end.field='circRNA_end',
                             keep.extra.columns=T)
}