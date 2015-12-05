library(shiny)
library(DT)
library(magrittr)
library(ggbio)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(dplyr)

source('functions.R')
GeneRanges<-loadGeneRanges()
load('rmsk_0.0.1.RData')
load('rmsk.family.RData')
c('circRNA_ID', 'sample', 'X.junction_reads', 'SM_MS_SMS', 'X.non_junction_reads', 
  'junction_reads_ratio', 'junction.Normal', 'junction.Tumor',
  'non_junction.Normal', 'non_junction.Tumor', 'ratio.Normal',
  'ratio.Tumor', 'p.values', 'fdr') ->
  sample_columns

shinyServer(function(input, output, session) {
  output$ciri_files<-DT::renderDataTable({
    ciri_files<-dir(input$ciri_path)
    data.frame(files=paste(input$ciri_path,ciri_files,sep='/')) %>%
      datatable
  })
  
  ciri_list<-reactive({
    ciri_files<-dir(input$ciri_path)
    paste(input$ciri_path,ciri_files,sep='/') %>%
      loadCIRI
  })
  
  ciri_rbind<-reactive({
    df<-do.call(rbind, ciri_list()) %>%
      mutate(length = circRNA_end - circRNA_start)
  })
  
  #ciri_rbind_filter? could speed up table rebuild when data is large or changing extend_size
  #ui could use dynamic ui with add remove filter by column class
  
  ciri_rbind_gr<-reactive({
    ciri_rbind() %>% df2GRanges
  })
  
  ciri_rbind_anno<-reactive({
    annotateGene(ciri_rbind_gr(), GeneRanges)
  })
  
  ciri_rbind_rmsk<-reactive({
    progress <- shiny::Progress$new(session)
    on.exit(progress$close())
    progress$set(message = 'Annotate Repeats',
                 detail = 'This may take a while...')
    
    if(input$anno_repeat){
      circRmsk<-annotateRmsk(ciri_rbind_gr())
      flankRmsk<-annotateRmsk(flank_both(ciri_rbind_gr(), input$extend_size))
      
      colnames(circRmsk)<-gsub('repeat', 'circRepeat', colnames(circRmsk))
      colnames(flankRmsk)<-gsub('repeat', 'flankRepeat', colnames(flankRmsk))
      
      merge(circRmsk, flankRmsk, by='circRNA_ID', all=T)
    }else{
      data.frame(circRNA_ID=factor())
    }
  })
  
  ciri_rbind_rank<-reactive({
    ciri_rbind() %>% rankCircRNA
  })
  
  ciri_merged<-reactive({
    ciri_rbind() %>%
      #left_join(ciri_rbind_anno()) %>%
      #left_join(ciri_rbind_rmsk()) %>%
      left_join(ciri_rbind_rank())
  })
  
  output$ciri_datatable<-DT::renderDataTable({
    ciri_merged() %>%
      datatable
  })
  
  ciri_selected<-reactive({
    row_id<-input$ciri_datatable_rows_selected
    ciri_merged()[row_id,]
  })
  
  ciri_selected_sample<-reactive({
    ciri_selected()[sample_columns] %>%
      rows2df
  })
  
  output$rows_sample_table<-DT::renderDataTable({
    ciri_selected_sample() %>%
      datatable
  })
  
  ciri_selected_circRNA<-reactive({
    selected<-ciri_selected()
    circRNA_columns<-c('circRNA_ID', setdiff(names(selected), sample_columns))
    selected[circRNA_columns] %>%
      unique %>%
      rows2df
  })
  
  output$rows_circRNA_table<-DT::renderDataTable({
    ciri_selected_circRNA() %>%
      datatable
  })
  
  output$helper<-renderText({
    #names(ciri_rbind())
  })
})