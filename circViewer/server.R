library(shiny)
library(DT)
library(magrittr)
library(ggbio)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(dplyr)

source('functions.R')

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
    do.call(rbind, ciri_list())
  })
  
  output$helper<-renderText({
    names(ciri_rbind())
  })
})