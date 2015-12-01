library(ggbio)
library(dplyr)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)

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
    dplyr::select(circRNA_ID, chr, circRNA_start, circRNA_end, circRNA_type,
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

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
data(genesymbol, package = "biovizBase")

pdf('plots/arc.filtered.noZero.pdf', height=12, width=12)
for (symbol in c('BRAF', 'FOXO1', 'EP300', 'ATRX')){
  message(symbol)
  #ciri_rbind_anno %>%
  rbind(ciri_rbind_anno_more_normal,ciri_rbind_anno_more_tumor) %>%
    plotTrack(symbol, reads_cnt=1) %>%
    print
}
dev.off()
