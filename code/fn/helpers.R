# Project: E(tive)Lice
# Tim Szewczyk
# tim.szewczyk@sams.ac.uk
# Helper functions


seq_range <- function(x, ...) {
  x_min <- min(x)
  x_max <- max(x)
  seq(x_min, x_max, ...)
}



make_compositional <- function(x, method="softmax") {
  switch(method,
         "proportional" = x/sum(x),
         "softmax" = exp(x)/sum(exp(x)),
         paste0("Error: Method must be 'proportional' or 'softmax', but was '", method, "'")
  )
}




render_qmd <- function(input_file, output_path, file_ext, ...) {
  # Extract just the input file name (without the file-extension)
  file_name <- xfun::sans_ext(input_file)

  # render the input document and output file will be in the
  # current working directory.
  quarto::quarto_render(input = input_file, output_format = file_ext, ...)

  # name of the rendered output file
  output_name <- paste0(file_name, ".", file_ext)

  # move the file to the output path
  fs::file_move(paste0(output_name), paste0(output_path, "sim_", str_sub(output_path, -3, -2), ".", file_ext))

  msg <- paste0(paste0(output_name, collapse = " and "), " moved to ", output_path)
  message(msg)
}
