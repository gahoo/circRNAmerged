library(shiny)
library(DT)
library(magrittr)
library(ggbio)
library(org.Hs.eg.db)
library(dplyr)
library(lazyeval)
library(d3heatmap)
library(RColorBrewer)


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
      mutate(length = circRNA_end - circRNA_start,
             ratio.Diff=ratio.Tumor-ratio.Normal)
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
      left_join(ciri_rbind_anno()) %>%
      #left_join(ciri_rbind_rmsk()) %>%
      #left_join(ciri_rbind_db()) %>%
      left_join(ciri_rbind_rank()) %>%
      filter_(.dots = filtering[['criteria']] )
  })
  
  output$ciri_datatable<-DT::renderDataTable({
    ciri_merged() %>%
      datatable(filter='top', options=list(stateSave = TRUE) )
  })
  
  ciri_column_class<-reactive({
    ciri_merged() %>%
      sapply(class) %>%
      as.list
  })
  
  output$ciri_filtering_column<-renderUI({
    selectInput('ciri_filtering_column', 'Filter',
                choices = c('free', names(ciri_merged())),
                selected = 'circRNA_ID')
  })
  
  filtering<-reactiveValues(criteria=list())
  
  observeEvent(input$add_filter,{
    addCriteria(input$ciri_filtering_column, 
                input$ciri_filter_string,
                ciri_column_class(),
                filtering[['criteria']]) ->
      filtering[['criteria']]
  })
  
  observeEvent(input$remove_filter,{
    n<-length(filtering[['criteria']])
    filtering[['criteria']][[n]]<-NULL
  })
  
  observeEvent(input$clear_filter,{
    filtering[['criteria']]<-list()
  })
  
  output$criteria<-renderText({
    filtering[['criteria']] %>% paste(collapse = ', ')
  })
  
  ciri_selected<-reactive({
    if(input$subsettingBy=='rows'){
      row_id<-input$ciri_datatable_rows_selected
    }else if(input$subsettingBy=='pages'){
      row_id<-input$ciri_datatable_rows_current
    }
    
    if(input$subsettingBy=='none'){
      ciri_merged()
    }else{
      if(is.null(row_id)){
        row_id<-1
      }
      
      selected<-ciri_merged()[row_id,] %>% fixSymbol
      
      columnName<-input$showBy
      filterValues<-unique(selected[[columnName]])
      filter_criteria <- interp(~ columnName %in% filterValues,
                                columnName=as.name(columnName))
      ciri_merged() %>%
        fixSymbol %>%
        filter_(filter_criteria)
    }
  })
  
  selected_rows_circRNA<-reactive({
    transformerFunc<-transformer(input$col2row)
    selected<-ciri_selected()
    circRNA_columns<-c('circRNA_ID', setdiff(names(selected), sample_columns))
    selected[circRNA_columns] %>%
      unique %>%
      mutate(gene_id=build_a(gene_id)) %>%
      transformerFunc
  })
  
  output$rows_circRNA_table<-DT::renderDataTable({
    selected_rows_circRNA() %>%
      datatable(escape = F)
  })
  
  selected_rows_samples<-reactive({
    transformerFunc<-transformer(input$col2row)
    ciri_selected()[sample_columns] %>%
      transformerFunc
  })
  
  output$rows_sample_table<-DT::renderDataTable({
    selected_rows_samples() %>%
      datatable
  })
  
  output$upset_samples<-renderPlot({
    ciri_selected() %>%
      plotSampleSets
  })
  
  output$upset_circRNA<-renderPlot({
    ciri_selected() %>%
      plotCircSets
  })
  
  output$ratio_pattern<-renderPlot({
    ciri_selected() %>%
      plotRelExpPattern(
        facet=input$ratio_pattern_facet,
        significant2alpha=input$ratio_pattern_map_significant,
        line=input$ratio_pattern_line,
        p.value=input$ratio_pattern_p)
  })
  
  output$ratio_heatmap<-renderD3heatmap({
    if(input$diff_ratio){
      colors_scheme = rev(brewer.pal(3,"RdYlGn"))
    }else{
      colors_scheme = 'Blues'
    }
    
    ciri_selected() %>%
      prepareHeatmapRatio(diff=input$diff_ratio, color_fix=T) %>%
      d3heatmap(colors = colors_scheme,
                dendrogram = input$d3heatmap_dendrogram,
                scale = input$d3heatmap_scale,
                symm = input$d3heatmap_symm,
                xaxis_height=150, yaxis_width=350)
  })
  
  output$ratio_heatmap2<-renderPlot({
    if(input$diff_ratio){
      colors_scheme = function(n){colorpanel(n, "green", "yellow", "red")}
      isSymbreaks <- T
    }else{
      colors_scheme = function(n){colorpanel(n, "white", "blue")}
      isSymbreaks <- F
    }
    
    ciri_selected() %>%
      prepareHeatmapRatio(
        diff=input$diff_ratio,
        color_fix=F) %>%
      as.matrix %>%
      heatmap.2(col = colors_scheme,
                symbreaks = isSymbreaks,
                dendrogram = input$d3heatmap_dendrogram,
                scale = input$d3heatmap_scale,
                symm = input$d3heatmap_symm,
                margins = c(input$heatmap2_height, input$heatmap2_width),
                trace = 'none')
  })
  
  output$arc_plot<-renderPlot({
    ciri_selected() %>%
      plotTrack(
        plot.transcript = input$track_transcript,
        plot.repeats = input$track_repeats,
        extend_size = input$track_extend_size,
        flank_only = input$track_repeats_flank_only,
        repeat_column = input$track_repeats_column,
        repeat.y = input$track_repeats_y,
        repeat.fill = input$track_repeats_fill
        )
  })
  
  output$helper<-renderText({
    #str(input$ciri_datatable_rows_all)
    #str(input$ciri_datatable_search_columns)
    #str(filtering[['criteria']])
    #input$add_filter
  })
})