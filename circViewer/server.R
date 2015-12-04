library(shiny)
library(DT)
library(magrittr)
library(ggbio)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(dplyr)

source('functions.R')
GeneRanges<-loadGeneRanges()

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
    df<-do.call(rbind, ciri_list())
    df %>%
      left_join(annotateDf(df2GRanges(df), GeneRanges)) %>%
      left_join(rankCircRNA(df))
  })
  
  output$ciri_table<-DT::renderDataTable({
    ciri_rbind() %>%
      datatable
  })
  
  output$helper<-renderText({
    names(ciri_rbind())
  })
})