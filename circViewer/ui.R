library(shiny)
library(DT)
library(d3heatmap)


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
                 DT::dataTableOutput('ciri_datatable'),
                 selectInput('showBy', 'by:', choices = c('circRNA_ID', 'symbol'),
                             selected='circRNA_ID'),
                 checkboxInput('col2row', 'long table', value = F)
  ),
  collapsibleDiv(id='selected_rows_circRNA_table', collapse = T,
                 label = 'circRNA',
                 class = 'btn-info btn-xs',
                 DT::dataTableOutput('rows_circRNA_table')
  ),
  collapsibleDiv(id='selected_rows_sample_table', collapse = T,
                 label = 'Samples',
                 class = 'btn-info btn-xs',
                 DT::dataTableOutput('rows_sample_table')
  ),
  collapsibleDiv(id='selected_rows_set_samples', collapse = T,
                 label = 'SetsSamples',
                 class = 'btn-info btn-xs',
                 plotOutput('upset_samples')
  ),
  collapsibleDiv(id='selected_rows_set_circRNA', collapse = T,
                 label = 'SetsCircRNA',
                 class = 'btn-info btn-xs',
                 plotOutput('upset_circRNA')
  ),
  collapsibleDiv(id='selected_rows_heatmap', collapse = F,
                 label = 'ratioHeatmap',
                 class = 'btn-info btn-xs',
                 checkboxInput('diff_ratio', 'Diff Ratio', value = F),
                 checkboxInput('d3heatmap_symm', 'symm', value = F),
                 selectInput('d3heatmap_dendrogram', 'dendrogram',
                             choices=c('none', 'row', 'column', 'both'),
                             selected='both'),
                 selectInput('d3heatmap_scale', 'scale',
                             choices=c('none', 'row', 'column'),
                             selected='none'),
                 d3heatmapOutput('ratio_heatmap', height='700px')
  )
  
))