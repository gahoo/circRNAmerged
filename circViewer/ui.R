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
  collapsibleDiv(id='info', collapse = T,
                 label = 'Settings',
                 class = 'btn-info btn-xs',
                 icon = icon('info-sign',lib='glyphicon'),
                 textInput('ciri_path', 'CIRI ouput path', value = '../CIRI'),
                 DT::dataTableOutput('ciri_files')
  )
))