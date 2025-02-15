#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("Please supply an input dir!", call.=FALSE)
} else if (length(args)==5) {
  # Load required packages with error handling
  required_packages <- c("plotly", "tidyverse", "vcfR", "DT", "cli", "flexdashboard")
  for(pkg in required_packages) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      install.packages(pkg)
      library(pkg, character.only = TRUE)
    }
  }

  cli_alert_info("Processing coverage files...")
  
  # Read split coverage files
  files <- list.files(args[1], pattern = "_splits.coverage.tsv", full.names = T, recursive = T)
  
  if(length(files) == 0) {
    cli_alert_warning("No coverage files found!")
    return()
  }

  # Process coverage data with progress bar
  cli_progress_bar("Reading coverage files", total = length(files))
  cov <- tibble()
  for(f in files) {
    cli_alert_info(paste("Reading file:", f))
    # Read the raw data first
    tmp <- read_tsv(f, col_names = "raw_data", show_col_types = FALSE)
    
    # print("Initial raw data:")
    # print(head(tmp))
    
    # Split the single column into proper columns
    tmp <- tmp %>%
      separate(raw_data, 
               into = c("Name", "CHR", "POS", "Cov"),
               sep = "\\s+",    # Split on whitespace
               extra = "drop") %>%
      mutate(Reads = ifelse(grepl("split", Name), "split", "all"))
    
    # print("After splitting columns:")
    # print(head(tmp))
    
    cov <- bind_rows(cov, tmp)
    cli_progress_update()
  }
  cli_progress_done()
  
  # print("Combined data before cleanup:")
  # print(head(cov))
  # print(colnames(cov))
  
  # Clean up coverage data and ensure numeric columns
  cov <- cov %>%
    separate(Name, c("Isolate"), sep = "_(?=[^_]+$)", extra = "drop") %>%
    mutate(
      POS = as.numeric(as.character(POS)),
      Cov = as.numeric(as.character(Cov))
    )

  # print("Final data structure:")
  # print(head(cov))
  # print(str(cov))

  cli_alert_info("Creating plots...")

  # Create split reads plot with proper layout
  p1 <- cov %>% 
    filter(Reads == "split") %>% 
    plot_ly(data = ., x = ~POS, y = ~Cov, color = ~Isolate, 
            type = 'scatter', mode = 'lines', legendgroup = ~Isolate) %>% 
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
      xaxis = list(title = "Position"),
      yaxis = list(title = "Coverage"),
      showlegend = TRUE
    )

  # Create all reads plot
  p2 <- cov %>% 
    filter(Reads == "all") %>% 
    plot_ly(data = ., x = ~POS, y = ~Cov, color = ~Isolate, 
            type = 'scatter', mode = 'lines', legendgroup = ~Isolate) %>% 
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
      xaxis = list(title = "Position"),
      yaxis = list(title = "Coverage"),
      showlegend = FALSE
    )

  # Combine plots side by side
  p3 <- subplot(p1, p2, shareY = T)

  # Create proportion plot
  p4 <- cov %>% 
    pivot_wider(names_from = Reads, values_from = Cov, values_fill = 0) %>% 
    mutate(prop = split/all) %>% 
    plot_ly(data = ., x = ~POS, y = ~prop, color = ~Isolate, 
            type = 'scatter', mode = 'lines', legendgroup = ~Isolate) %>% 
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
      xaxis = list(title = "Position"),
      yaxis = list(title = "Proportion"),
      showlegend = TRUE
    )

  # Create final combined plot with proper styling
  p_final <- subplot(p3, p4, nrows = 2, margin = 0.05) %>%
    layout(
      plot_bgcolor = "#3A3A39",
      paper_bgcolor = "#3A3A39",
      xaxis = list(
        tickcolor = "#FEFBE8",
        gridcolor = "rgba(254, 251, 232, 0.2)",
        linecolor = "#FEFBE8"
      ),
      yaxis = list(
        tickcolor = "#FEFBE8",
        gridcolor = "rgba(254, 251, 232, 0.2)",
        linecolor = "#FEFBE8"
      ),
      legend = list(
        font = list(
          family = "sans-serif",
          size = 12,
          color = "#FEFBE8"
        ),
        bgcolor = "#3A3A39",
        bordercolor = "#3A3A39"
      ),
      font = list(color = '#FEFBE8')
    )

  # Process gene coordinates
  cli_alert_info("Processing gene coordinates...")
  gene_start <- as.numeric(args[2])
  gene_end <- as.numeric(args[3])

  if (gene_start > gene_end) {
    tmp <- gene_start
    gene_start <- gene_end
    gene_end <- tmp
  }

  # Find and process relevant files
  cli_alert_info("Processing additional data files...")
  
  # Find relevant files
  is_files <- list.files(path = args[1], pattern = "_IS6110_table.txt$", full.names = T, recursive = T)
  ss_files <- list.files(path = args[1], pattern = "_gridss.vcf.gz$", full.names = T, recursive = T)
  lp_files <- list.files(path = args[1], pattern = "_sample.splitters.unsorted.bam$", full.names = T, recursive = T)
  
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
  
  cli_alert_success("Dashboard data generation complete!")
}
