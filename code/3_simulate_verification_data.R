# Project: E(tive)Lice
# Tim Szewczyk
# tim.szewczyk@sams.ac.uk
# Simulate verification data


library(tidyverse)
library(glue)
dir("code/fn", full.names=T) |> walk(source)

n_sims <- 30


for(i in 1:n_sims) {
  sim_number <- str_pad(length(dir("data/sim", "sim"))+1, 2, 'left', '0')
  render_qmd("code/3_simulate_verification_data_noGravid.qmd",
             output_path=glue("data/sim/sim_{sim_number}/"),
             file_ext="pdf")
}

