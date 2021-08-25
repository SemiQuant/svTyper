#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
knitr::opts_chunk$set(dev="CairoPNG")

if (length(args)==0) {
  stop("Please supply an input dir!", call.=FALSE)
} else if (length(args)==5) {
  require(plotly, quietly = T)
  require(tidyverse, quietly = T)
  require(vcfR, quietly = T)
  require(DT, quietly = T)
  files <- list.files(args[1], pattern = "_splits.coverage.tsv", full.names = T, recursive = T)
  
  # cov <- files %>%
  #   map(read_tsv, col_names = F) %>%
  #   reduce(rbind) %>% 
  #   separate(X1, c("Isolate", "Reads"), sep = "_")
  # colnames(cov)[3:5] <- c("CHR", "POS", "Cov")
  
  # above was problematic with sample nanimgs so
  cov <- files %>%
    map(read_tsv, col_names = F) %>%
    reduce(rbind) %>% 
    mutate(Reads = ifelse(grepl("split", X1), "split", "all")) %>% 
    separate(X1, c("Isolate"), sep = "(?:.(?!_))+$")
  
  colnames(cov) <- c("Isolate", "CHR", "POS", "Cov", "Reads")
  
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
  
  
  
  p_final <- subplot(p3, p4, nrows = 2, margin = 0.05) %>%
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
  # htmlwidgets::saveWidget(as_widget(p_final), "splits.html") 
  
  
  

  
  gene_start <- args[2]
  gene_end <-args[3]
  is_files <- list.files(path = args[1], pattern = "_IS6110_table.txt$", full.names = T, recursive = T)
  ss_files <- list.files(path = args[1], pattern = "_gridss.vcf.gz$", full.names = T, recursive = T)
  lp_files <- list.files(path = args[1], pattern = "_lumpy_SVs.denovo.vcf.gz$", full.names = T, recursive = T)
  
  
  
  
  
  if (gene_start > gene_end){
    tmp <- gene_start
    gene_start <- gene_end
    gene_end <- tmp
  }
  
  
  is_tib <- data_frame(Indiv = gsub("_IS6110_table.txt", "", basename(is_files))) %>% 
    mutate(file_contents = map(is_files, ~ read_tsv(file.path(.)))) %>% 
    unnest %>% 
    filter((x >= gene_end & x <= gene_start) | (y <= gene_end & y >= gene_start)) %>% # this make sense?
    select(Indiv, x, y)
  
  try({
    cbind(is_tib, "IS element")
    colnames(is_tib) <- c("ID", "Start", "End", "IS elemnt")
    
  })
  
  
  ss_tib <- ss_files %>% 
    map(function(x) vcfR2tidy(read.vcfR(x))$gt) %>% 
    map_dfr(bind_rows) %>% 
    select(POS, Indiv, gt_AF) %>% 
    filter((POS >= gene_end & POS <= gene_start) | (POS <= gene_end & POS >= gene_start)) # this doesnt make sense
  
  
  lp_tib <- data_frame(Indiv = gsub("_lumpy_SVs.denovo.vcf.gz", "", basename(lp_files))) %>% 
    mutate(file_contents = map(lp_files, ~ vcfR2tidy(read.vcfR(file.path(.)))$gt)) %>% 
    unnest
  lp_tib <- lp_tib %>% 
    select(Indiv, POS, gt_AD, gt_GT_alleles) %>% 
    filter((POS >= gene_end & POS <= gene_start) | (POS <= gene_end & POS >= gene_start)) # this doesnt make sense
  
  
  save.image(file = paste0(args[1], "/", args[4], "_make_dash_dat.Rimg"))
  
  
  rmarkdown::render(
    input = paste0(args[5], "/make_dash.Rmd"),
    params = list(filename  = paste0(args[1], "/", args[4], "_make_dash_dat.Rimg")),
    output_file = paste0(args[4], "_dash.html"),
    output_dir = args[1],
    clean = T
    )
  
  # Twas a bitch to get this working in an old container 
  # apt-get install libgtk2.0-dev libcairo2-dev xvfb xauth xfonts-base libxt-dev r-cran-cairodevic
  # set knitr::opts_chunk$set(dev="CairoPNG")
  # install.packages("Cairo", dep = T)
  # install.packages("rmarkdown", dep = T)
  
}