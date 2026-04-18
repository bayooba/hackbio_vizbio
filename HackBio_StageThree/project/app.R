library(shiny)
library(ggplot2)
library(dplyr)
library(pheatmap)

meta_df <- read.csv('cell_metadata.csv')
expr_mat <- read.csv('expression_matrix.csv')

rownames(expr_mat) <- expr_mat[, 1]
expr_mat <- expr_mat[, -1]

umap_cord <- read.csv('umap_coordinates.csv')
merged_data <- left_join(meta_df, umap_cord, by = 'cell_id')

# helper/boiler functions
# 1. colors with adjustable transparency
t_col <- function(color,
                  percent = 50,
                  name = NULL) {
  # color: named color or hex code
  # percent: transparency percentage (0 = opaque, 100 = fully transparent)
  # name: optional name for the resulting color
  
  rgb.val <- col2rgb(color)
  
  t.col <- rgb(
    rgb.val[1],
    rgb.val[2],
    rgb.val[3],
    max = 255,
    alpha = (100 - percent) * 255 / 100,
    names = name
  )
  
  invisible(t.col)
}

# scale values: convert numeric values to a 0-100 range for consistent visualization
scale_0_100 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  
  if (rng[1] == rng[2]) {
    return(rep(50, length(x)))
  }
  
  (x - rng[1]) / (rng[2] - rng[1]) * 100
}

# compute per-gene statistics for a selected cell type
compute_gene_statistics <- function(expr_mat, meta_df, target_ct) {
  # identify cells inside and outside the selected cell type
  in_cells <- meta_df$cell_id[meta_df$cell_type == target_ct]
  out_cells <- meta_df$cell_id[meta_df$cell_type != target_ct]
  
  # subset expression matrix
  xin <- expr_mat[in_cells, , drop = FALSE]
  xout <- expr_mat[out_cells, , drop = FALSE]
  
  # detection rates
  det_in <- colMeans(xin > 0)
  det_out <- colMeans(xout > 0)
  
  # Mean expression
  mean_in <- colMeans(xin)
  mean_out <- colMeans(xout)
  
  #specificity difference
  diff <- mean_in - mean_out
  
  return(
    data.frame(
      gene = colnames(expr_mat),
      det_in = det_in,
      det_out = det_out,
      mean_in = mean_in,
      mean_out = mean_out,
      diff = diff,
      stringsAsFactors = FALSE
    )
  )
}

# this is the user interface (UI) component of the app
ui <- fluidPage(titlePanel("Simulated Cell Type Viewer"),
                
                sidebarLayout(
                  sidebarPanel(
                    selectInput(
                      inputId = 'cell_type',
                      label = 'Select Cell Type',
                      choices = c('All', unique(merged_data$cell_type)),
                      selected = 'All'
                    )
                  ),
                  
                  mainPanel(
                    plotOutput(outputId = 'scatter_plot'),
                    uiOutput(outputId = 'marker_text'),
                    conditionalPanel(condition = "input.cell_type != 'All'", tableOutput(outputId = 'per_gene_stat'))
                  )
                ))

# this is the server side (logic) component of the app; this part defines and controls what is displayed on the user interface component.
server <- function(input, output) {
  data_to_plot <- reactive({
    if (input$cell_type == 'All') {
      my_current_subset <- merged_data
    }
    
    else {
      my_current_subset <- subset(merged_data, cell_type == input$cell_type)
    }
    
    return(
      list(
        'x_value' = my_current_subset$UMAP_1,
        'y_value' = my_current_subset$UMAP_2,
        'col_value' = my_current_subset$cell_type
      )
    )
  })
  
  color_cells <- c(
    'T_cell' = '#6CB4EE',
    'B_cell' = '#1B1B1B',
    'NK_cell' = '#50C878',
    'Monocyte' = '#d56'
  )
  
  color_cells_gene_exp <- c(
    'T_cell' = '#0059CF',
    'B_cell' = '#100C08',
    'NK_cell' = '#00563B',
    'Monocyte' = '#660000'
  )
  
  output$scatter_plot <- renderPlot({
    if (input$cell_type == 'All') {
      plot_df <- as.data.frame(data_to_plot())
      ggplot(plot_df, aes(x = x_value, y = y_value, colour = col_value))+
        geom_point()+
        scale_color_manual(values = c(
          'T_cell' = '#6CB4EE',
          'B_cell' = '#1B1B1B',
          'NK_cell' = '#50C878',
          'Monocyte' = '#d56'
        ))+
        labs(
          x = 'UMAP_1',
          y = 'UMAP_2',
          color = 'Cell Type'
        )+
        theme_minimal()+
        theme(
          legend.title = element_text(size = 14, face = 'bold'),
          legend.text = element_text(size = 12),
          axis.text = element_text(size = 10),
          axis.title = element_text(size = 12)
        )
      # plot(
      #   x = data_to_plot()$x_value,
      #   y = data_to_plot()$y_value,
      #   col = NULL,
      #   bg = color_cells[data_to_plot()$col_value],
      #   pch = 21,
      #   cex = 1.2,
      #   xlab = 'UMAP_1',
      #   ylab = 'UMAP_2'
      # )
    } else {
      stat_result <- compute_gene_statistics(expr_mat, meta_df, input$cell_type)
      
      best_gene <- stat_result$gene[which.max(stat_result$diff)]
      
      expr_values <- expr_mat[meta_df$cell_id, best_gene]
      
      col_pallete <- sapply(scale_0_100(expr_values), function(x) {
        x <- (x / 100)^0.5 * 100
        t_col(color_cells_gene_exp[data_to_plot()$col_value], 100 - x)
      })
      
      plot_df <- as.data.frame(data_to_plot())
      rownames(plot_df) <- meta_df[meta_df$cell_type == input$cell_type, 'cell_id']
      plot_df$expression <- expr_mat[rownames(plot_df), best_gene]
      
      ggplot(plot_df, aes(x = x_value, y = y_value, color = expression)) +
        geom_point() +
        scale_color_gradient(
          low = 'white',
          high = color_cells_gene_exp[input$cell_type],
          name = 'Expression level'
        )+
        labs(x = 'UMAP_1', y = 'UMAP_2')+
        theme_minimal()+
        theme(
          #legend.title.align = 0.5, ::: alternative to hjust
          legend.title = element_text(size = 14, hjust = 0.5),
          
          #legend.text.align = 0.5,
          legend.text = element_text(size = 12),
          legend.key.width = unit(0.8, 'cm'),
          legend.key.height = unit(1.5, 'cm'),
          axis.text = element_text(size = 10)
        )
      
      # plot(
      #   x = data_to_plot()$x_value,
      #   y = data_to_plot()$y_value,
      #   col = col_pallete,
      #   #bg = expr_values,
      #   pch = 16,
      #   cex = 1.2,
      #   xlab = 'UMAP_1',
      #   ylab = 'UMAP_2'
      # )
      
    }
    
  })
  
  output$per_gene_stat <- renderTable({
    stat_result <- compute_gene_statistics(expr_mat, meta_df, input$cell_type)
    head(stat_result[order(-stat_result$diff), ], 7)
  })
  
  output$marker_text <- renderUI({
    if (input$cell_type == 'All') {
      return ("")
    }
    
    stat_result <- compute_gene_statistics(expr_mat, meta_df, input$cell_type)
    best_gene <- stat_result$gene[which.max(stat_result$diff)]
    
    HTML(
      paste0(
        'The marker gene for ',
        input$cell_type,
        ' is ',
        '<b style = "color: red">',
        best_gene,
        '</b>',
        ' with a diff value of ',
        '<b style = "color: red">',
        sprintf('%.2f', max(stat_result$diff)),
        '<b>',
        '<br>',
        '&nbsp;' # tool a while to get the "&nbsp;" code that added an empty line break to the page. R!!!!!
      )
    )
  })
}

# Run the application
shinyApp(ui = ui, server = server)
