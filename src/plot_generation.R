# Existing code...

# Fix the data reading
df <- read_delim(input_file, 
                 delim="\t",
                 col_names=TRUE,
                 col_types=cols(.default = "c")) %>%
  separate_rows(everything(), sep="\t")

# Fix the plotly layout
p <- plot_ly(df, ...) %>%
  layout(
    shapes = list(
      list(
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
  )

# Add progress bar
with_progress({
  p <- ggplotly(p)
  message("Plot generated successfully")
})

# Existing code... 