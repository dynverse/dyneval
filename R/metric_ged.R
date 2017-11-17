
#' Compute the graph edit distance using gedevo
#'
#' @param net1 the first network to compare
#' @param net2 the second network to compare
#'
#' @importFrom readr read_file
#' @importFrom utils write.table
#' @importFrom GEDEVO run_GEDEVO
calculate_ged <- function(net1, net2) {
  net1 <- net1 %>% mutate(dir="u") %>% select(from, dir, to)
  net2 <- net2 %>% mutate(dir="u") %>% select(from, dir, to)

  tempfolder <- tempfile()
  dir.create(tempfolder)
  write.table(net1, file.path(tempfolder, "net1.sif"), row.names = FALSE, col.names = FALSE)
  write.table(net2, file.path(tempfolder, "net2.sif"), row.names = FALSE, col.names = FALSE)

  cmd <- pritt("gedevo --groups a b --sif {tempfolder}/net1.sif a --sif {tempfolder}/net2.sif b --no-prematch --no-workfiles --save {tempfolder}/out --maxiter 100 --maxsecs 1 --maxsame 100")
  GEDEVO::run_GEDEVO(cmd)

  score <- readr::read_file(paste0(tempfolder, "/out.matching")) %>% gsub("^.*GED score:[ ]*([0-9\\.]*).*", "\\1", .) %>% as.numeric()

  1-score
}


