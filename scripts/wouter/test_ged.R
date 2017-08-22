library(tidyverse)


net1 <- tribble(
  ~from, ~to,
  1, 2
)

net2 <- tribble(
  ~from, ~to,
  1, 2,
  2, 3,
  4, 5
)


calculate_ged <- function(net1, net2) {
  net1 <- net1 %>% mutate(dir="u") %>% select(from, dir, to)
  net2 <- net2 %>% mutate(dir="u") %>% select(from, dir, to)

  tempfolder <- tempdir()
  write.table(net1, file.path(tempfolder, "net1.sif"), row.names = FALSE, col.names = FALSE)
  write.table(net2, file.path(tempfolder, "net2.sif"), row.names = FALSE, col.names = FALSE)


  system(glue::glue("./inst/extra_code/GEDEVO/linux-x64/gedevo --groups a b --sif {tempfolder}/net1.sif a --sif {tempfolder}/net2.sif b --no-prematch --no-workfiles --save {tempfolder}/out"), ignore.stdout=TRUE)

  score <- read_file(paste0(tempfolder, "/out.txt")) %>% gsub("^.*GED score:[ ]*([0-9\\.]*).*", "\\1", .) %>% as.numeric()

  1-score
}


calculate_ged(net2, net1)

