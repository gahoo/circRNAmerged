library(shiny)
library(DT)

collapsibleDiv<-function(id, ..., label='Show/Hide', .func=actionButton,
                         collapse = FALSE, class=NULL, icon=NULL, width=NULL){
  
  collapse_status<-ifelse(collapse, "on", "in")
  
  list(
    .func(
      sprintf("b_%s",id), label=label, icon=icon, class=class, width=width,
      "data-toggle"='collapse', "data-target"=sprintf('#%s',id)
    ),
    div(
      id = id, class = sprintf("collapse %s", collapse_status),
      ...
    )
  )
}

shinyUI(fluidPage(
  titlePanel("circRNA BAM Viewer"),
  textOutput('helper'),
  collapsibleDiv(id='settings', collapse = T,
                 label = 'Settings',
                 class = 'btn-info btn-xs',
                 icon = icon('info-sign',lib='glyphicon'),
                 textInput('ciri_path', 'CIRI ouput path', value = '../CIRI'),
                 checkboxInput('anno_repeat', 'annotate repeat', value = T),
                 numericInput('extend_size', 'extend size for repeat annotation:',
                              value = '2000', min = 0, max = 100000),
                 DT::dataTableOutput('ciri_files')
  ),
  collapsibleDiv(id='ciri_table', collapse = F,
                 label = 'CIRI',
                 class = 'btn-info btn-xs',
                 DT::dataTableOutput('ciri_datatable')
  ),
  collapsibleDiv(id='ciri_sample_table', collapse = F,
                 label = 'circRNA',
                 class = 'btn-info btn-xs',
                 DT::dataTableOutput('rows_circRNA_table'),
                 DT::dataTableOutput('rows_sample_table')
  )
  
))