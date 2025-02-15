#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("Please supply an input dir!", call.=FALSE)
} else if (length(args)==5) {
  require(plotly, quietly = T)
  require(tidyverse, quietly = T)
  require(vcfR, quietly = T)
  require(DT, quietly = T)
  
  # Read split coverage files
  files <- list.files(args[1], pattern = "_splits.coverage.tsv", full.names = T, recursive = T)
  
  # Process coverage data
  cov <- files %>%
    map(read_tsv, col_names = F) %>%
    reduce(rbind) %>% 
    mutate(Reads = ifelse(grepl("split", X1), "split", "all")) %>% 
    separate(X1, c("Isolate"), sep = "(?:.(?!_))+$")
  
  colnames(cov) <- c("Isolate", "CHR", "POS", "Cov", "Reads")
  
  # Create split reads plot
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
  
  # Create all reads plot
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
  
  # Combine plots side by side
  p3 <- subplot(p1, p2, shareY = T)
  
  # Create proportion plot
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
  
  # Create final combined plot
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
             font = list(
               family = "sans-serif",
               size = 12,
               color = "#FEFBE8"),
             bgcolor = "#3A3A39",
             bordercolor = "#3A3A39"
           ),
           font = list(color = '#FEFBE8')
    )
  
  # Process gene coordinates
  gene_start <- as.numeric(args[2])
  gene_end <- as.numeric(args[3])
  
  # Find relevant files
  is_files <- list.files(path = args[1], pattern = "_IS6110_table.txt$", full.names = T, recursive = T)
  ss_files <- list.files(path = args[1], pattern = "_gridss.vcf.gz$", full.names = T, recursive = T)
  lp_files <- list.files(path = args[1], pattern = "_sample.splitters.unsorted.bam$", full.names = T, recursive = T)
  
  if (gene_start > gene_end) {
    tmp <- gene_start
    gene_start <- gene_end
    gene_end <- tmp
  }
  
  # Process ISMap results
  if (length(is_files) > 0) {
    is_tib <- tibble(Indiv = gsub("_IS6110_table.txt", "", basename(is_files))) %>% 
      mutate(file_contents = map(is_files, ~ read_tsv(file.path(.)))) %>% 
      unnest %>% 
      filter((x >= gene_start & x <= gene_end) | (y <= gene_end & y >= gene_start)) %>%
      select(Indiv, x, y)
    
    colnames(is_tib) <- c("ID", "Start", "End", "IS element")
  } else {
    is_tib <- tibble(ID = character(), Start = numeric(), End = numeric())
  }
  
  # Process GRIDSS results
  if (length(ss_files) > 0) {
    ss_tib <- ss_files %>% 
      map(function(x) vcfR2tidy(read.vcfR(x))$gt) %>% 
      map_dfr(bind_rows) %>% 
      select(POS, Indiv, gt_AF) %>% 
      filter(POS >= gene_start & POS <= gene_end)
  } else {
    ss_tib <- tibble(POS = numeric(), Indiv = character(), gt_AF = numeric())
  }
  
  # Process split read alignment results
  if (length(lp_files) > 0) {
    lp_tib <- tibble(Indiv = gsub("_sample.splitters.unsorted.bam", "", basename(lp_files))) %>%
      select(Indiv)
  } else {
    lp_tib <- tibble(Indiv = character())
  }
  
  # Save workspace
  save.image(file = paste0(args[1], "/", args[4], "_make_dash_dat.Rimg"))
  
  # Render dashboard
  rmarkdown::render(
    input = paste0(args[5], "/make_dash.Rmd"),
    params = list(filename = paste0(args[1], "/", args[4], "_make_dash_dat.Rimg")),
    output_file = paste0(args[4], "_dash.html"),
    output_dir = args[1],
    clean = T
  )
}
