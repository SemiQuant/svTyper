library(cli)
library(tidyverse)

create_progress_table <- function(steps) {
  steps_df <- tibble(
    Step = steps,
    Status = "Pending",
    Time = NA
  )
  
  cli::cli_table(steps_df)
}

update_progress <- function(table, step, status, time) {
  table %>%
    mutate(
      Status = if_else(Step == step, status, Status),
      Time = if_else(Step == step, time, Time)
    )
} 