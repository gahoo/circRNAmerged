library(shiny)
library(DT)
library(magrittr)
library(ggbio)
library(org.Hs.eg.db)
library(dplyr)
library(lazyeval)
library(d3heatmap)
library(pheatmap)
library(RColorBrewer)

source('functions.R')

c('circRNA_ID', 'sample', 'X.junction_reads', 'SM_MS_SMS', 'X.non_junction_reads', 
  'junction_reads_ratio', 'junction.Normal', 'junction.Tumor',
  'non_junction.Normal', 'non_junction.Tumor', 
  'depth.Normal', 'depth.Tumor', 'depth.Total',
  'ratio.Normal', 'ratio.Tumor', 'ratio.Diff', 
  'p.values', 'fdr') ->
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
      loadCIRI(nrow=input$preview_nrow)
  })
  
  ciri_rbind<-reactive({
    df<-do.call(rbind, ciri_list()) %>%
      mutate(length = circRNA_end - circRNA_start,
             ratio.Diff = ratio.Tumor-ratio.Normal,
             depth.Normal = junction.Normal + non_junction.Normal,
             depth.Tumor = junction.Tumor + non_junction.Tumor,
             depth.Total = depth.Normal + depth.Tumor)
  })
  
  ciri_rbind_gr<-reactive({
    ciri_rbind() %>% df2GRanges
  })
  
  ciri_rbind_anno<-reactive({
    progress <- shiny::Progress$new(session)
    on.exit(progress$close())
    progress$set(message = 'Annotate Genes',
                 detail = 'This may take a while...')
    
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
    progress <- shiny::Progress$new(session)
    on.exit(progress$close())
    progress$set(message = 'Calc Rank',
                 detail = 'This may take a while...')
    
    ciri_rbind() %>% rankCircRNA
  })
  
  ciri_merged<-reactive({
    progress <- shiny::Progress$new(session)
    on.exit(progress$close())
    progress$set(message = 'Building table',
                 detail = 'This may take a while...')
    
    ciri_rbind() %>%
      left_join(ciri_rbind_anno()) %>%
      left_join(ciri_rbind_rmsk()) %>%
      #left_join(ciri_rbind_db()) %>%
      left_join(ciri_rbind_rank())
  })
  
  ciri_merged_filter<-reactive({
    progress <- shiny::Progress$new(session)
    on.exit(progress$close())
    progress$set(message = 'Filtering table',
                 detail = 'This may take a while...')
    
    if(length(filtering[['criteria']])>0){
      ciri_merged() %>%
        filter_(.dots = filtering[['criteria']] )
    }else{
      ciri_merged()
    }
  })
  
  output$downloadTableData <- downloadHandler(
    filename = function() {
      dir(input$ciri_path) %>%
        gsub(pattern = '.CIRI.merged',
             replacement = '') %>%
        paste0(collapse = '_') ->
        ciri_basename
      paste(format(Sys.time(), "%Y-%b-%d_%X"), '.', ciri_basename, '.filtered.csv.gz', sep='')
    },
    content = function(con) {
      gzip <- gzfile(con, "w")
      ciri_merged_filter() %>%
        write.csv(file=gzip, row.names=F)
      close(gzip)
    }
  )
  
  output$ciri_datatable<-DT::renderDataTable({
    ciri_merged_filter() %>%
      datatable_template
  })
  
#   proxy = DT::dataTableProxy('ciri_datatable')
#   
#   observeEvent(input$clear_selection,{
#     selectColumns(proxy, NULL)
#   })
  
  ciri_column_class<-reactive({
    ciri_merged_filter() %>%
      sapply(class) %>%
      as.list
  })
  
  output$ciri_filtering_column<-renderUI({
    selectInput('ciri_filtering_column', 'Filter',
                choices = c('free', names(ciri_merged_filter())),
                selected = 'free')
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
    #update after reload
    input$preview_nrow
    filtering[['criteria']] %>% paste(collapse = ', ')
  })
  
  ciri_selected<-reactive({
    if(input$filter_only){
      ciri_conductor<-ciri_merged_filter
    }else{
      ciri_conductor<-ciri_merged
    }
    
    if(input$subsettingBy=='rows'){
      row_id<-input$ciri_datatable_rows_selected
    }else if(input$subsettingBy=='pages'){
      row_id<-input$ciri_datatable_rows_current
    }
    
    if(input$subsettingBy=='none'){
      ciri_conductor()
    }else{
      if(is.null(row_id)){
        row_id<-1
      }
      
      ciri_conductor() %>%
        selectedFilter(
          filterValues=row2filterValues(
            ciri_merged_filter(),
            row_id = row_id,
            columnName = input$showBy),
          columnName=input$showBy)
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
      datatable_template2(escape = F)
  })
  
  selected_rows_samples<-reactive({
    transformerFunc<-transformer(input$col2row)
    ciri_selected()[sample_columns] %>%
      transformerFunc
  })
  
  output$rows_sample_table<-DT::renderDataTable({
    selected_rows_samples() %>%
      datatable_template2
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
    progress <- shiny::Progress$new(session)
    on.exit(progress$close())
    progress$set(message = 'Ploting Relative Expression Ratio Pattern',
                 detail = 'This may take a while...')
    
    ciri_selected() %>%
      plotRelExpPattern(
        facet=input$ratio_pattern_facet,
        significant2alpha=input$ratio_pattern_map_significant,
        line=input$ratio_pattern_line,
        p.value=input$ratio_pattern_p)
  })
  
  output$ratio_d3heatmap<-renderD3heatmap({
    if(input$diff_ratio){
      colors_scheme = rev(brewer.pal(3,"RdYlGn"))
    }else{
      colors_scheme = 'Blues'
    }
    
    ciri_selected() %>%
      prepareHeatmapRatio(diff=input$diff_ratio, color_fix=T, rownames_fix=T) %>%
      d3heatmap(colors = colors_scheme,
                dendrogram = input$d3heatmap_dendrogram,
                scale = input$d3heatmap_scale,
                symm = input$d3heatmap_symm,
                xaxis_height=150, yaxis_width=350)
  })
  
  output$ratio_pheatmap<-renderPlot({
    if(input$diff_ratio){
      colors_scheme = colorRampPalette(c("green", "yellow", "red"))(41)
      breaks = seq(-1,1,by = 0.05)
      isSymbreaks <- T
    }else{
      colors_scheme = colorRampPalette(c("white", "blue"))(21)
      breaks = seq(0,1,by = 0.05)
      isSymbreaks <- F
    }
    
    ciri_selected() %>%
      prepareHeatmapRatio(
        diff=input$diff_ratio,
        color_fix=F,
        rownames_fix=F) %>%
      as.matrix %>%
      pheatmap(
        color = colors_scheme,
        breaks = breaks,
        scale = input$d3heatmap_scale,
        cluster_rows = input$pheatmap_cluster_rows,
        cluster_cols = input$pheatmap_cluster_cols
        )
  })
  
  output$arc_plot<-renderPlot({
    progress <- shiny::Progress$new(session)
    on.exit(progress$close())
    progress$set(message = 'Ploting Tracks',
                 detail = 'This may take a while...')
    
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
  
  output$downloadPlotData <- downloadHandler(
    filename = function() {
      dir(input$ciri_path) %>%
        gsub(pattern = '.CIRI.merged',
             replacement = '') %>%
        paste0(collapse = '_') ->
        ciri_basename
      paste(format(Sys.time(), "%Y-%b-%d_%X"), '.', ciri_basename, '.pdf', sep='')
    },
    content = function(con) {
      pdf(con)
      input$batch_ids %>%
        strsplit(split='\n') ->
        ids
      plotAllFig(ids, ciri_merged_filter())
      dev.off()
    })
  
  output$tbl_plot_ctrls<-renderUI({
    choices <- names(ciri_merged_filter())
    additionnal.choices<-c('NULL', choices)
    exclude.columns<-setdiff(additionnal.choices, c())
      
    list(
      column(6, style='background: rgba(255, 255, 255, 0.9);',
        selectInput('tbl_plot_x', 'X', choices = choices, selected='ratio.Diff'),
        selectInput('tbl_plot_y', 'Y', choices = choices, selected='p.values'),
        checkboxInput('tbl_plot_x_log', 'log(X)', F),
        checkboxInput('tbl_plot_y_log', 'log(Y)', F),
        checkboxInput('tbl_plot_flip', 'flip x y', F),
        selectInput('tbl_plot_func', 'function',
                    choices = c('geom_point', 'geom_bar','geom_boxplot'),
                    selected='geom_point'),
        textInput('tbl_plot_facet', 'facet by',
                  value = '. ~ sample')
        ),
      column(6, style='background: rgba(255, 255, 255, 0.9);',
        selectInput('tbl_plot_fill', 'fill', 
                    choices = additionnal.choices, selected='NULL'),
        selectInput('tbl_plot_color', 'color', 
                    choices = additionnal.choices, selected='NULL'),
        selectInput('tbl_plot_alpha', 'alpha', 
                    choices = additionnal.choices, selected='NULL'),
        selectInput('tbl_plot_group', 'group', 
                    choices = additionnal.choices, selected='NULL'),
        selectInput('tbl_plot_size', 'size',
                    choices = additionnal.choices, selected='NULL')
        )
    )
  })

  output$table_plot<-renderPlot({
    ciri_selected() %>%
      plotTable(
        x=input$tbl_plot_x,
        y=input$tbl_plot_y,
        log_x=input$tbl_plot_x_log,
        log_y=input$tbl_plot_y_log,
        flip=input$tbl_plot_flip,
        func=input$tbl_plot_func,
        fill=input$tbl_plot_fill,
        color=input$tbl_plot_color,
        alpha=input$tbl_plot_alpha,
        group=input$tbl_plot_group,
        size=input$tbl_plot_size,
        facet=input$tbl_plot_facet)
  })

  output$rows_summary_table<-DT::renderDataTable({
    ciri_selected()%>%
      summaryTblNumCols ->
      dt
    if(input$col2row){
      row.names(dt)<-dt$sample
      dt <- t(dt[,-1]) %>% as.data.frame
    }
    dt %>%
      datatable
  })

  output$helper<-renderText({
    #str(input$ciri_datatable_rows_all)
    #str(input$ciri_datatable_search_columns)
    #str(filtering[['criteria']])
    #input$add_filter
  })
})