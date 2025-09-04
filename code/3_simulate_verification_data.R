# Project: E(tive)Lice
# Tim Szewczyk
# tim.szewczyk@sams.ac.uk
# Simulate verification data


library(tidyverse)
library(glue)
dir("code/fn", ".R$", full.names=T) |> walk(source)

n_sims <- 20


for(i in 1:n_sims) {
  sim_number <- str_pad(length(dir("data/sim", "sim"))+1, 2, 'left', '0')
  render_qmd("code/3_simulate_verification_data.qmd",
             output_path=glue("data/sim/sim_{sim_number}/"),
             file_ext="html")
}

