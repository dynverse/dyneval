library(tidyverse)
library(googlesheets)

script_file <- "scripts/robrecht/scholar.py"
download.file("https://raw.githubusercontent.com/ckreibich/scholar.py/master/scholar.py", destfile = script_file)

num_citations_by_clusterid <- function(clusterid, scholar_file = script_file) {
  tryCatch({
    command <- paste0("python ", scholar_file, " -C ", clusterid, " --csv-header")
    output <- system(command, intern = T)
    tab <- readr::read_delim(paste(gsub("\n", " ", output), collapse = "\n"), delim = "|")
    sum(tab$num_citations)
  }, error = function(e) NA)
}

method_df <- gs_key("1Mug0yz8BebzWt8cmEW306ie645SBh_tDHwjVw4OFhlE") %>%
  gs_read(col_types = cols(GScholarClusterID = "c"))

method_df <- method_df %>%
  filter(Pseudotime) %>%
  mutate(Citations = pbapply::pbsapply(GScholarClusterID, num_citations_by_clusterid)) %>%
  arrange(desc(Citations))

method_df %>% select(Name, PubDate, Preprint, Journal, Citations) %>% as.data.frame
