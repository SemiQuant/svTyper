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
    plot_ly(x = ~POS, y = ~Cov, color = ~Isolate, type = 'scatter', mode = 'lines')
  
  p2 <- cov %>% 
    filter(grepl("all", Reads)) %>% 
    plot_ly(x = ~POS, y = ~Cov, color = ~Isolate, type = 'scatter', mode = 'lines')
  
  
  # plot 1, side by side of all and split
  p3 <- subplot(p1, p2, shareY = T)
  
  
  # plot 2, percentage of split vs all
  p4 <- cov %>% 
    pivot_wider(names_from = Reads, values_from = Cov, values_fill = 0) %>% 
    mutate(prop = split/all) %>% 
    plot_ly(x = ~POS, y = ~prop, color = ~Isolate, type = 'scatter', mode = 'lines')
  
  # orca(subplot(p3, p4, nrows = 2), "splits.pdf")
  htmlwidgets::saveWidget(as_widget(subplot(p3, p4, nrows = 2)), "splits.html")
}