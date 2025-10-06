
#save.image(file = "/tmp/Rdata")

# redirect output and messages/errors to the log
log <- file(snakemake@log[[1]], open="wt")
sink(log, type = "output")
sink(log, type = "message")


# get the run_dir,  and then create a microhaplot directory
# there with the Shiny elements in it.




bams <- snakemake@input$bams
bams_dir <- dirname(bams[1])

vcf <- snakemake@input$input_vcf

outfile <- snakemake@output$rds

num_threads <- snakemake@threads[1]



library(tidyverse)
library(microhaplotextract)

# get the short file names for the bams, in a tibble
bams_tibble <- tibble(
  bams = basename(bams),
  sample = str_replace(bams, "\\.bam", "")
)

# first create the label file by parsing the units.csv
# joining NMFS_DNA_ID on there and then writing it out
# into the bams directory
run_dir_with_species <- str_replace(bams_dir, "/bams.*$", "")
run_dir <- dirname(run_dir_with_species)
marker_set <- str_replace(basename(outfile), "--.*$", "")
label_file_path <- file.path(bams_dir, "labels.txt")
app_path <- file.path(run_dir_with_species, "microhaplot")
out_path <- tempfile()
dir.create(out_path, recursive = TRUE, showWarnings = FALSE)

units <- read_csv(file.path(run_dir, "units.csv")) %>%
  filter(Markers == marker_set)

samples <- read_csv(file.path(run_dir, "samples.csv")) %>%
  select(sample, NMFS_DNA_ID, Sample_Project)

labels_tibble <- left_join(units, samples, by = "sample") %>%
  select(sample, NMFS_DNA_ID, Sample_Project) %>%
  mutate(sample = paste0(sample, ".bam"))

write.table(
  labels_tibble,
  file = label_file_path,
  quote = FALSE,
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE
)


# set the run label to be the basename of the output file minus the rds
run_label = basename(outfile) %>%
  str_replace("\\.rds$", "")


haplo.read.tbl <- prepHaplotFiles(
  run.label = run_label,
  bam.path = normalizePath(bams_dir),
  out.path = out_path,
  label.path = label_file_path,
  vcf.path = vcf,
  app.path = app_path,
  n.jobs = num_threads
)
