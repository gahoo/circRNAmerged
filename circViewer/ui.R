library(shiny)
library(DT)
library(shinyAce)
library(d3heatmap)

preset_filter<-'p.values<0.05 & depth.Total > 10 & 
(depth.Normal <=1 & depth.Tumor >= 5) &
(depth.Normal >=5 & depth.Tumor <= 1)'
inlineDiv<-function(...){div(style="display:inline-block", ...)}

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

absoluteCollaspablePanel<-function(
  ..., draggable=T, fixed=F, left=NULL, right=NULL, top=NULL, bottom=NULL){
  absolutePanel(
    collapsibleDiv(...), fixed=fixed,
    draggable=draggable, left=left, right=right, top=top, bottom=bottom)
}

fixedCollaspablePanel<-function(
  ..., draggable=T, fixed=T, left=NULL, right=5, top=10, bottom=NULL,
  style='background: rgba(255, 255, 255, 0.6);'){
  fixedPanel(
    collapsibleDiv(..., style=style),
    draggable=draggable, left=left, right=right, top=top, bottom=bottom)
}

collapsibleDiv(id='selected_rows_summary_table', collapse = T,
               label = 'summary',
               class = 'btn-info btn-xs',
               uiOutput('summary_columns'),
               DT::dataTableOutput('rows_summary_table')) ->
  summary_table

collapsibleDiv(id='config', collapse = T,
               label = 'config',
               class = 'btn-danger btn-xs pull-right',
               icon = icon('info-sign',lib='glyphicon'),
               textInput('ciri_path', 'CIRI ouput path', value = '../CIRI'),
               DT::dataTableOutput('ciri_files')) ->
  config

collapsibleDiv(
  id='filters', collapse = T,
  label = 'filters',
  class = 'btn-warning btn-xs pull-right',
  style='background: rgba(255, 255, 255, 0.9);',
  icon = icon('info-sign',lib='glyphicon'),
  numericInput('load_nrow', 'load rows:',
               min = -1, max = 1/0, value=10),
  checkboxInput('preview_table', 'preview', value = T),
  numericInput('preview_nrow', 'preview rows:',
               min = -1, max = 1/0, value=1000),
  checkboxInput('anno_repeat', 'annotate repeat', value = F),
  numericInput('extend_size', 'extend size for repeat annotation:',
               value = '2000', min = 0, max = 100000),
  uiOutput('ciri_filtering_column'),
  textInput('ciri_filter_string', '', preset_filter),
  actionButton('add_filter','Add'),
  actionButton('remove_filter','Remove'),
  actionButton('clear_filter','Clear'),
  downloadButton('downloadTableData', 'Download'),
  verbatimTextOutput('filter_nrow'),
  verbatimTextOutput('criteria')
)->filtering

collapsibleDiv(
  id='subsetting', collapse = T,
  label = 'subsetting',
  class = 'btn-success btn-xs pull-right',
  style='background: rgba(255, 255, 255, 0.9);',
  icon = icon('info-sign',lib='glyphicon'),
  selectInput('subsettingBy', 'subsetting by:',
              choices = c('rows', 'pages', 'none'),
              selected = 'rows'),
  selectInput('showBy', 'show by:',
              choices = c('circRNA_ID', 'symbol'),
              selected='circRNA_ID'),
  checkboxInput('filter_only', 'filtered data only', value = T),
  checkboxInput('col2row', 'colums 2 rows', value = F)
  #actionButton('clear_selection','Clear Selection'),
) -> subsetting

collapsibleDiv(id='ciri_table', collapse = F,
               label = 'CIRI',
               class = 'btn-info btn-xs',
               DT::dataTableOutput('ciri_datatable'),
               fixedPanel(
                 draggable=T, top=40, right=5,
                 filtering,
                 subsetting
               )) ->
  ciri_table

collapsibleDiv(id='selected_rows_circRNA_table', collapse = T,
               label = 'circRNA',
               style='background: #FFFFFF;',
               class = 'btn-primary btn-xs',
               DT::dataTableOutput('rows_circRNA_table', width='800px')
) ->
  circRNA_table

collapsibleDiv(id='selected_rows_sample_table', collapse = T,
               label = 'Samples',
               style='background: #FFFFFF;',
               class = 'btn-primary btn-xs',
               DT::dataTableOutput('rows_sample_table', width='800px')
) ->
  sample_table

collapsibleDiv(id='selected_rows_symbol_hpa', collapse = T,
               label = 'HPA_Cancer',
               style='background: #FFFFFF;',
               class = 'btn-primary btn-xs',
               selectInput('hpa_position', 'position',
                           choices=c('fill', 'stack', 'dodge'),
                           selected='fill'),
               plotOutput('hpa_cancer_symbols')
) ->
  hpa_cancer

collapsibleDiv(id='selected_rows_plot_table', collapse = T,
               label = 'tablePlot',
               class = 'btn-info btn-xs',
               fixedCollaspablePanel(
                 id='plot_table_controls', collapse = F,
                 label = 'plotTblControls',
                 class = 'btn-success btn-xs pull-right',
                 uiOutput('tbl_plot_ctrls')
               ),
               plotOutput('table_plot')
) ->
  table_plot

collapsibleDiv(id='selected_rows_set_samples', collapse = T,
               label = 'SetsSamples',
               class = 'btn-info btn-xs',
               plotOutput('upset_samples')
)->
  upset_sample

collapsibleDiv(id='selected_rows_set_circRNA', collapse = T,
               label = 'SetsCircRNA',
               class = 'btn-info btn-xs',
               plotOutput('upset_circRNA', )
)->
  upset_circRNA

collapsibleDiv(id='selected_rows_ratio_pattern', collapse = T,
               label = 'ratioPattern',
               class = 'btn-info btn-xs',
               fixedCollaspablePanel(
                 id='pattern_controls', collapse = F,
                 label = 'ratioPatternControls',
                 class = 'btn-success btn-xs pull-right',
                 icon = icon('info-sign',lib='glyphicon'),
                 checkboxInput('ratio_pattern_facet','facet',value=T),
                 checkboxInput('ratio_pattern_map_significant',
                               'significant dots only',value=T),
                 numericInput('ratio_pattern_p', 'P value <= ',
                              value = 0.05, step = 0.005,
                              min = 0, max = 1),
                 checkboxInput('ratio_pattern_line',
                               'connect dots',value=T)
               ),
               plotOutput('ratio_pattern', )
)->
  ratio_pattern

collapsibleDiv(id='selected_rows_heatmap', collapse = T,
               label = 'ratioHeatmap',
               class = 'btn-info btn-xs',
               fixedCollaspablePanel(
                 id='heatmap_controls', collapse = F,
                 label = 'heatmapControls',
                 class = 'btn-success btn-xs pull-right',
                 icon = icon('info-sign',lib='glyphicon'),
                 checkboxInput('diff_ratio', 'Diff Ratio', value = F),
                 checkboxInput('d3heatmap_symm', 'symm', value = F),
                 selectInput('d3heatmap_dendrogram', 'dendrogram',
                             choices=c('none', 'row', 'column', 'both'),
                             selected='both'),
                 selectInput('d3heatmap_scale', 'scale',
                             choices=c('none', 'row', 'column'),
                             selected='none'),
                 checkboxInput('pheatmap_cluster_cols',
                               'cluster columns', value = F),
                 checkboxInput('pheatmap_cluster_rows',
                               'cluster rows', value = F)
               ),
               tabsetPanel(
                 tabPanel('d3heatmap',
                          style='background: rgba(255, 255, 255, 0.9);',
                          d3heatmapOutput('ratio_d3heatmap')
                 ),
                 tabPanel('pheatmap',
                          plotOutput('ratio_pheatmap', )
                 )
               )
)->ratio_heatmap

collapsibleDiv(id='selected_rows_arc', collapse = T,
               label = 'arcPlot',
               class = 'btn-info btn-xs',
               fixedCollaspablePanel(
                 id='arc_controls', collapse = F,
                 label = 'arcControls',
                 class = 'btn-success btn-xs pull-right',
                 icon = icon('info-sign',lib='glyphicon'),
                 checkboxInput('track_transcript', 'plot transcript', value=F),
                 checkboxInput('track_repeats', 'plot repeats', value=F),
                 checkboxInput('track_repeats_flank_only', 'flank only', value=F),
                 numericInput('track_extend_size', 'extend size for repeat:',
                              value = '2000', min = 0, max = 100000),
                 selectInput('track_repeats_y', 'repeat.y', 
                             choices = c('name', 'class', 'family'),
                             selected = 'name'),
                 selectInput('track_repeats_column', 'repeat reverse complement by', 
                             choices = c('name', 'class', 'family'),
                             selected = 'name'),
                 selectInput('track_repeats_fill', 'repeat fill color by', 
                             choices = c('name', 'score', 'class', 'family', 'strand'),
                             selected = 'strand')
               ),
               plotOutput('arc_plot')
)->
  arc_plot

collapsibleDiv(id='batch_mode', collapse = T,
               label = 'batchMode',
               class = 'btn-info btn-xs',
               checkboxInput('batch_one_by_one', 'one by one', T),
               selectizeInput('batch_plots', 'plot:', multiple=T,
                              choices=c(
                                HPA_cancer = 'plotHPA',
                                tablePlot='plotTable',
                                SetsSamples='plotSampleSets',
                                SetsCircRNA='plotCircSets',
                                ratioPattern='plotRelExpPattern',
                                ratioHeatmap='plotRelExpPheatmap',
                                arcPlot='plotTrack'),
                              selected = c('plotRelExpPattern', 'plotTrack')
                              ),
               textInput('sequence_fa', 'indexed fasta:', value = '../chrALL.fa'),
               aceEditor("batch_ids", mode='txt', value="", height="200px"),
               downloadButton('downloadPlotData', 'Download PDF'),
               downloadButton('downloadFa', 'Download Fa')
)->
  batch_mode

collapsibleDiv(id='qqPvalues', collapse = T,
               label = 'qqPvalues',
               class = 'btn-info btn-xs',
               fixedCollaspablePanel(
                 id='qq_controls', collapse = F,
                 label = 'qqControls',
                 class = 'btn-success btn-xs pull-right',
                 icon = icon('info-sign',lib='glyphicon'),
                 selectInput('qq_column', 'column:',
                             choices = c('p.values', 'fdr', 'group.p.values', 'group.fdr'),
                             selected = 'group.p.values'),
                 checkboxInput('qq_facet', 'facet by sample', value = F)
                 ),
               plotOutput('qqP')
)->
  qqPvalues

shinyUI(fluidPage(
  titlePanel("circRNA Viewer"),
  tags$head(
    tags$script(src="jquery-ui.min.js"),
    tags$script(src="resize.js"),
    tags$link(rel = "stylesheet", type = "text/css", href = "jquery-ui.css")
    ),
  textOutput('helper'),
  config,
  ciri_table,
  fixedPanel(
    draggable=T, top=63, left=58, width='auto', height='auto', class='resizable',
    style='border: 1px solid grey;',
    table_plot,
    qqPvalues,
    upset_sample,
    upset_circRNA,
    ratio_pattern,
    ratio_heatmap,
    arc_plot
  ),
  fixedPanel(
    draggable=T, top=0, left=15, width='auto', height='auto', class='resizable',
    style='border: 1px solid grey;',
    circRNA_table,
    sample_table,
    hpa_cancer
  ),
  summary_table,
  batch_mode
))