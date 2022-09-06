#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("Please supply an input dir!", call.=FALSE)
} else if (length(args)==1) {
  require(plotly)
  require(tidyverse)
  files <- list.files(args[1], pattern = "_splits.coverage.tsv", full.names = T, recursive = T)
  
  cov <- files %>%
    map(read_tsv, col_names = F) %>%
    reduce(rbind) %>% 
    separate(X1, c("Isolate", "Reads"), sep = "_")
  
  colnames(cov)[3:5] <- c("CHR", "POS", "Cov")
  
  
  p1 <- cov %>% 
    filter(grepl("split", Reads)) %>% 
    plot_ly(x = ~POS, y = ~Cov, color = ~Isolate, type = 'scatter', mode = 'lines', legendgroup = ~Isolate) %>% 
    add_annotations(
      text = "SPLIT",
      x = 0.5,
      y = 1,
      yref = "paper",
      xref = "paper",
      xanchor = 0,
      yanchor = "bottom",
      showarrow = FALSE,
      font = list(size = 12)
    ) %>%
    layout(
      # showlegend = FALSE,
      shapes = list(
        type = "rect",
        x0 = 0,
        x1 = 1,
        xref = "paper",
        y0 = 0, 
        y1 = 16,
        yanchor = 1,
        yref = "paper",
        ysizemode = "pixel",
        fillcolor = toRGB("gray"),
        line = list(color = "transparent")
      )
    ) 
  
  p2 <- cov %>% 
    filter(grepl("all", Reads)) %>% 
    plot_ly(x = ~POS, y = ~Cov, color = ~Isolate, type = 'scatter', mode = 'lines', legendgroup = ~Isolate) %>% 
    add_annotations(
      text = "ALL",
      x = 0.5,
      y = 1,
      yref = "paper",
      xref = "paper",
      xanchor = 0,
      yanchor = "bottom",
      showarrow = FALSE,
      font = list(size = 12)
    ) %>%
    layout(
      # showlegend = FALSE,
      shapes = list(
        type = "rect",
        x0 = 0,
        x1 = 1,
        xref = "paper",
        y0 = 0, 
        y1 = 16,
        yanchor = 1,
        yref = "paper",
        ysizemode = "pixel",
        fillcolor = toRGB("gray"),
        line = list(color = "transparent")
      )
    ) 
  
  
  # plot 1, side by side of all and split
  p3 <- subplot(p1, p2, shareY = T)
  
  
  # plot 2, percentage of split vs all
  p4 <- cov %>% 
    pivot_wider(names_from = Reads, values_from = Cov, values_fill = 0) %>% 
    mutate(prop = split/all) %>% 
    plot_ly(x = ~POS, y = ~prop, color = ~Isolate, type = 'scatter', mode = 'lines', legendgroup = ~Isolate) %>% 
    add_annotations(
      text = "SPLIT:All",
      x = 0.5,
      y = 1,
      yref = "paper",
      xref = "paper",
      xanchor = 0,
      yanchor = "bottom",
      showarrow = FALSE,
      font = list(size = 12)
    ) %>%
    layout(
      # showlegend = FALSE,
      shapes = list(
        type = "rect",
        x0 = 0,
        x1 = 1,
        xref = "paper",
        y0 = 0, 
        y1 = 16,
        yanchor = 1,
        yref = "paper",
        ysizemode = "pixel",
        fillcolor = toRGB("gray"),
        line = list(color = "transparent")
      )
    ) 
  
  
  
  p_final <- subplot(p3, p4, nrows = 2) %>%
    layout(plot_bgcolor="#3A3A39",
           paper_bgcolor ="#3A3A39",
           xaxis = list(tickcolor = "#FEFBE8",
                        gridcolor = "rgba(254, 251, 232, 0.2)",
                        linecolor = "#FEFBE8"),
           yaxis = list(tickcolor = "#FEFBE8",
                        gridcolor = "rgba(254, 251, 232, 0.2)",
                        linecolor = "#FEFBE8"),
           legend = list(
             # x = 0.6, y = 0.4,
             font = list(
               family = "sans-serif",
               size = 12,
               color = "#FEFBE8"),
             bgcolor = "#3A3A39",
             bordercolor = "#3A3A39"
             # borderwidth = 2
           ),
           font = list(color = '#FEFBE8')
    )
  
  # orca(p_final, "splits.pdf")
  htmlwidgets::saveWidget(as_widget(p_final), "splits.html") 
}