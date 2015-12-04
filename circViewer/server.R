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
  
  ciri_rbind_gr<-reactive({
    ciri_rbind() %>% df2GRanges
  })
  
  ciri_rbind_anno<-reactive({
    annotateDf(ciri_rbind_gr(), GeneRanges)
  })
  
  ciri_rbind_rank<-reactive({
    ciri_rbind() %>% rankCircRNA
  })
  
  output$ciri_table<-DT::renderDataTable({
    ciri_rbind() %>%
      left_join(ciri_rbind_anno()) %>%
      left_join(ciri_rbind_rank()) %>%
      datatable
  })
  
  output$helper<-renderText({
    #names(ciri_rbind())
  })
})