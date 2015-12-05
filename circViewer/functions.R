library(UpSetR)
library(GenomicRanges)
library(dplyr)

loadCIRI<-function(ciri_files){
  withProgress(message = 'Loading CIRI.merge files',
               detail = 'This may take a while...', value = 0, {
    n<-length(ciri_files)
    step_size<-1/n
    lapply(ciri_files, function(ciri_file){
      incProgress(step_size, detail=ciri_file)
      read.table(ciri_file, sep='\t',header=T,
                 nrow=100
                 ) %>%
        mutate(sample=gsub('.*/','',gsub('.CIRI.merged','',ciri_file))) ->
        data
      colnames(data)[11:16]<-c('junction.Normal', 'junction.Tumor',
                               'non_junction.Normal', 'non_junction.Tumor',
                               'ratio.Normal', 'ratio.Tumor')
      data
    })
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

loadGeneRanges<-function(){
  GeneID_SYMBOL<-AnnotationDbi::select(org.Hs.eg.db, 
                                       keys = keys(org.Hs.eg.db),
                                       columns = c("ENTREZID", "SYMBOL", "GENENAME")) %>%
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
                             keep.extra.columns=T)
}


annotateRmsk<-function(dfRanges){
  hits<-findOverlaps(dfRanges, rmsk)
  cbind(mcols(dfRanges[queryHits(hits)]),
        as.data.frame(rmsk[subjectHits(hits)])) %>%
    select(circRNA_ID, strand, name) %>%
    unique %>%
    as.data.frame %>%
    group_by(circRNA_ID, name) %>%
    summarise(cnt=n()) %>%
    filter(cnt==2) %>%
    select(circRNA_ID, name) %>%
    summarise(
      region_repeat_type_cnt=n(),
      region_repeat_name=paste0(name, collapse = ',')
    )
}

annotateGene<-function(dfRanges, GeneRanges){
  hits<-findOverlaps(dfRanges, GeneRanges)
  cbind(mcols(dfRanges[queryHits(hits)]),
        mcols(GeneRanges[subjectHits(hits)])) %>%
    unique %>%
    as.data.frame %>%
    group_by(circRNA_ID) %>%
    summarise(
      region_gene_id=paste0(GENEID, collapse = ','),
      region_symbol=paste0(SYMBOL, collapse = ','),
      region_gene_name=paste0(GENENAME, collapse = ';'),
      region_gene_cnt=n()
    )
}

extend<-function(gr, width, start=T, end=T){
  if(start){
    start(gr)<-start(gr)-width
  }
  if(end){
    end(gr)<-end(gr)+width
  }
  gr
}

flank_both<-function(gr, width){
  c(flank(gr, start=T, width),
    flank(gr, start=F, width))
}

rankCircRNA<-function(df) {
  abs_diff<-function(x){sqrt(sum(x^2)/length(x))}
  
  df %>%
    select(circRNA_ID, ratio.Normal, ratio.Tumor) %>%
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
    )
}