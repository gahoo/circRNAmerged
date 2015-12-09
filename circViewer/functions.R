library(UpSetR)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(dplyr)
library(tidyr)

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
data(genesymbol, package = "biovizBase")
load('rmsk_0.0.1.RData')
load('rmsk.family.RData')

progressTip<-function(session, message='Processing', detail='This may take a while...'){
  progress <- shiny::Progress$new(session)
  on.exit(progress$close())
  progress$set(message = message,
               detail = detail)
}

datatable_template<-function(data, ...){
  datatable(data, filter='top',
            extensions = c('TableTools','FixedHeader','ColReorder','ColVis'),
            options=list(
              pageLength = 10,
              stateSave = FALSE,
              dom = 'CRT<"clear">lfrtip',
              tableTools = list(sSwfPath = copySWF())
            ),
            ...
  )
}

loadCIRI<-function(ciri_files){
  withProgress(message = 'Loading CIRI.merge files',
               detail = 'This may take a while...', value = 0, {
    n<-length(ciri_files)
    step_size<-1/n
    lapply(ciri_files, function(ciri_file){
      incProgress(step_size, detail=ciri_file)
      read.table(ciri_file, sep='\t',header=T,
                 nrow=10000
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
    select(circRNA_ID, sample) %>%
    mutate(value = 1) %>%
    spread(sample, value, fill=0)
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
    filter(occurrence >= cnt) %>%
    "$"('circRNA_ID') %>%
    as.character %>%
    unique
}


addSymbolAnno<-function(df){
  df %>%
    select(circRNA_ID, symbol, region_symbol) %>%
    unique %>%
  mutate_each(funs(as.character)) %>%
  mutate(
    region_symbol = ifelse(is.na(region_symbol), '', region_symbol),
    anno = sprintf(
      "%s\n%s", circRNA_ID,
      ifelse(is.na(symbol), region_symbol, symbol))
    ) %>%
  select(circRNA_ID, anno)
}

plotCircSets<-function(df, cnt=0, ...){
  candidates<-filterCnt(df, cnt)
  df_set<-getSetDf(df)
  df_set_t<-df_set
  df_set %>%
    left_join(addSymbolAnno(df)) ->
    df_set
  row.names(df_set_t)<-unique(df_set$anno)
  df_set %>%
    filter(circRNA_ID %in% candidates) %>%
    "$"("anno") -> candidates
  #row.names(df_set_t)<-df_set$circRNA_ID
  df_set_t<-as.data.frame(t(df_set_t[,-1]))
  upset(df_set_t, sets = candidates, nsets = length(candidates), ...)
}



plotRelExpPattern<-function(df, facet=T, significant2alpha=T, line=T, p.value=0.05){
  df %>%
    select(ratio.Normal, ratio.Tumor, sample, circRNA_ID, p.values) %>%
    gather(type,ratio,ratio.Normal,ratio.Tumor) %>%
    mutate(type=gsub('ratio.','',type),
           significant=p.values<=p.value,
           log10P=log10(p.values) ) %>%
    left_join(addSymbolAnno(df)) %>%
    ggplot(aes(x=sample, y=ratio, group=type)) + 
    ylab('Relative Expression Ratio') +
    theme(axis.text.x=element_text(angle=90)) ->
    p
  
  if(significant2alpha&facet){
    point_aes<-aes(size=-log10P, alpha=significant, color=type, shape=type)
    line_aes<-aes(color=type)
    p <- p + facet_grid(anno~.)
  }else if(significant2alpha&!facet){
    point_aes<-aes(size=-log10P, alpha=significant, color=anno, shape=type)
    line_aes<-aes(color=circRNA_ID, group=factor(paste0(circRNA_ID,type)))
  }else if(!significant2alpha&facet){
    point_aes<-aes(size=-log10P, color=type, shape=type)
    line_aes<-aes(color=type)
    p <- p + facet_grid(anno~.)
  }else if(!significant2alpha&!facet){
    point_aes<-aes(size=-log10P, color=anno, shape=type)
    line_aes<-aes(color=circRNA_ID, group=factor(paste0(circRNA_ID,type)))
  }
  
  if(line){
    p<-p + geom_line(line_aes)
  }
  
  p + geom_point(point_aes)
}

prepareHeatmapRatio<-function(df, diff=F, color_fix=T, rownames_fix=T){
  getRatioDf<-function(df){
    df %>%
      select(circRNA_ID, sample, ratio.Normal, ratio.Tumor) %>%
      gather(type, value, starts_with("ratio")) %>% 
      mutate(type=gsub('ratio.','',type)) %>% 
      unite(sample_type, sample, type, sep='.') %>% 
      spread(sample_type, value, fill=0)
  }
  
  getRatioDiff<-function(df){
    df %>%
      select(circRNA_ID, sample, ratio.Diff) %>%
      spread(sample, ratio.Diff, fill=0)
  }
  
  anno<-addSymbolAnno(df)
  if(rownames_fix){
    anno <- anno %>% mutate(anno=gsub('\n','\t',anno))
  }
  
  if(diff){
    df <- getRatioDiff(df) 
    if(color_fix){
      df <- df %>% mutate(color_fix_min=-1,color_fix_max=1)
    }
  }else{
    df<-getRatioDf(df)
  }
  
  df %>%
    left_join(anno) ->
    df
  
  rownames(df)<-df$anno
  df %>% select(-circRNA_ID, -anno)
}

plotRelExpHeatmap<-function(df){
  d3heatmap(df, scale = "column", colors = "Spectral")
}

prepareArc<-function(df){
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
    select(circRNA_ID, chr, circRNA_start, circRNA_end, circRNA_type,
           junction.Normal, junction.Tumor, 
           non_junction.Normal,non_junction.Tumor,
           ratio.Normal, ratio.Tumor,p.values, sample) ->
    df
  rbind(
    df %>% prepare4rbind(type='Normal'),
    df %>% prepare4rbind(type='Tumor')
  ) %>%
    filter(junction > 0) %>%
    makeGRangesFromDataFrame(keep.extra.columns=T)
}

plotArc<-function(arc){
  if(length(arc)==0){
    message("empty")
    NULL
  }else{
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

fixSymbol<-function(df){
  df %>%
    mutate(
      symbol = as.character(symbol),
      region_symbol = as.character(region_symbol),
      symbol = ifelse(is.na(symbol), region_symbol, symbol)
    )
}

prepareRepeat<-function(df, extend_size=0, flank_only=F,
                        reverse_complement_only=T, 
                        repeat_column='name'){
  df2extendRanges<- function(df){
    if(flank_only){
      change_range<-flank_both
    }else{
      change_range<-extend
    }
    
    df %>%
      df2GRanges %>%
      change_range(extend_size)->
      which
  }
  
  reverseComplementFilter<-function(repeats){
    repeats %>%
      as.data.frame %>%
      select_('strand', repeat_column) %>%
      unique %>%
      "[["(repeat_column) ->
      candidates
    rev_comp_name<-candidates[duplicated(candidates)]
    
    mcols(repeats) %>%
      "[["(repeat_column) %>%
      "%in%"(rev_comp_name) ->
      rev_comp_idx
    
    repeats[rev_comp_idx]
  }
  
  which <- df2extendRanges(df)
  repeats <- subsetByOverlaps(rmsk, which, ignore.strand=T)
  mcols(repeats) %>%
    as.data.frame %>%
    left_join(rmsk.family) ->
    mcols(repeats)
  
  if(reverse_complement_only){
    reverseComplementFilter(repeats)
  }else{
    repeats
  }
  
}

plotRepeat<-function(repeatGranges, repeat.y, repeat.fill = 'strand'){
  mcols(repeatGranges)[[repeat.y]]<-as.character(mcols(repeatGranges)[[repeat.y]])
  if(length(repeatGranges)!=0){
    ggbio() + 
      geom_alignment(data = repeatGranges,
                     aes_string(group = repeat.y, fill = repeat.fill),
                     alpha = 0.7,
                     label = T)
  }else{
    NULL
  }
}

plotTrack<-function(df, plot.transcript=T, plot.repeats=T, ...,
                    repeat.y=repeat_column, repeat.fill = 'strand'){
  
  getSymbol<-function(df){
    df %>%
      select(symbol, region_symbol) %>%
      unique %>%
      fixSymbol %>%
      "$"("symbol") %>%
      strsplit(split=',') %>%
      unlist %>%
      unique
    
  }
  
  checkError<-function(df){
    n_circRNA<-length(unique(df$circRNA_ID))
    n_chr<-length(unique(df$chr))
    if(n_circRNA > 200){
      return(stop('too many circRNA to plot: ', n_circRNA))
    }
    if(n_chr != 1){
      return(
        stop('circRNA not in the same chromosome: ',
             paste0(unique(df$chr), collapse = ', ') )
      )
    }
  }
  
  error <- checkError(df)
  if(!is.null(error)){
    return(error)
  }
  
  symbol <- getSymbol(df)
  symbol <- symbol[!is.na(symbol)]
  
  df %>% prepareArc %>% plotArc -> arc
  
  if(plot.repeats){
    df %>%
      prepareRepeat(...) %>%
      plotRepeat(repeat.y=repeat.y,
                 repeat.fill) ->
      repeats
  }else{
    repeats <- NULL
  }
  
  if(plot.transcript & length(symbol) > 0 ){
    str(symbol)
    which_gene<-genesymbol[symbol]
    transcripts <- ggbio() + 
      geom_alignment(data=txdb, which = which_gene)    
    track_title <- paste0(symbol, collapse = "|")
  }else{
    transcripts <- NULL
    track_title <- paste0(unique(df$circRNA_ID), collapse = "; ")
  }
  
  track_list<-list(
    arc=arc,
    repeats=repeats,
    transcripts=transcripts
  )
  
  track_height<-c(arc=6, repeats=2, transcripts=2)
  track_list<-track_list[!sapply(track_list, is.null)]
  tracks(track_list, title=track_title, heights=track_height[names(track_list)])
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

countCircRNA<-function(df, columnName){
  df %>%
    select(circRNA_ID) %>%
    group_by(circRNA_ID) %>%
    summarise(Cnts=n()) %>%
    plyr::rename(c(Cnts=columnName))
}

addFamily<-function(df){
  df %>%
    select(circRNA_ID, strand, name) %>%
    unique %>%
    as.data.frame %>%
    left_join(rmsk.family)
}

countRepat<-function(df, columnName){
  df %>%
    group_by_('circRNA_ID', columnName) %>%
    summarise(cnt=n()) %>%
    filter(cnt>1) %>%
    select_('circRNA_ID', columnName) %>%
    summarise_(
      type_cnt="n()",
      region_repeat=sprintf("paste0(%s, collapse = ',')", columnName)
      ) %>%
    plyr::rename(c(
      type_cnt=paste0('repeat.',columnName,'Cnt'),
      region_repeat=paste0('repeat.',columnName))
      )
}

annotateRmsk<-function(dfRanges){
  hits<-findOverlaps(dfRanges, rmsk)
  overlapRmsk<-cbind(
    mcols(dfRanges[queryHits(hits)]),
    as.data.frame(rmsk[subjectHits(hits)])
    )
    
  repeatCnts<-countCircRNA(overlapRmsk, 'repeatCnt')
  overlapRmskFamily<-addFamily(overlapRmsk)
  repeatCnts %>%
    left_join(countRepat(overlapRmskFamily, 'name')) %>%
    left_join(countRepat(overlapRmskFamily, 'class')) %>%
    left_join(countRepat(overlapRmskFamily, 'family'))
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
  remove_na_zero<-function(x, v=1){ifelse(is.na(x)|x==0,v,x)}
  
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
              ratio.Normal.sd = remove_na_zero(ratio.Normal.sd),
              ratio.Tumor.sd = remove_na_zero(ratio.Tumor.sd),
              ratio.Diff.sd = remove_na_zero(ratio.Diff.sd),
              ratio.abs_diff = abs_diff(ratio.diff),
              ratio.rank = ratio.abs_diff/(ratio.Normal.sd * ratio.Tumor.sd * ratio.Diff.sd),
              ratio.rank2 = ratio.abs_diff/(ratio.Normal.sd * ratio.Tumor.sd),
              ratio.rank3 = ratio.abs_diff/ratio.Diff.sd
    )
}

rows2df<-function(df){
  lapply(1:nrow(df), function(i){
    df[i,] %>%
      as.list %>%
      lapply(as.character) %>%
      unlist %>%
      as.data.frame
  }) %>%
    do.call(what=cbind)
}

doNothing<-function(df){
  df
}

transformer<-function(transformTable){
  if(transformTable){
    rows2df
  }else{
    doNothing
  }
}

build_a<-function(x){
  ifelse(is.na(x), NA, sprintf('<a href="http://www.ncbi.nlm.nih.gov/gene/%s" target="_blank">%s</a>', x, x))
}

addCriteria<-function(column_name, filter_string, df_class, criteria){
  if(column_name!='free'){
    column_class<-df_class[[column_name]]
    if(column_class %in% c('factor', 'character')){
      filter_string<-sprintf('grepl("%s", as.character(%s))',
                             filter_string, column_name)
    }else if(column_class %in% c('numeric', 'integer')){
      filter_string<-sprintf('%s%s', column_name, filter_string)
    }else if(column_class %in% c('logical')){
      filter_string<-sprintf('%s%s', column_name, filter_string)
    }
  }
  
  if(filter_string %in% criteria){
    criteria
  }else{
    n<-length(criteria)+1
    criteria[[n]]<-filter_string
  }
  
  criteria
}