plot_strip <- function(task1, task2) {
  task1$milestone_network <- task1$milestone_network %>% mutate(cumlength = c(0, cumsum(length)[-length(length)]))
  task2$milestone_network <- task2$milestone_network %>% mutate(cumlength = c(0, cumsum(length)[-length(length)]))

  prog1 <- task1$progression %>% left_join(task1$milestone_network, by=c("from", "to")) %>% mutate(cumpercentage=percentage*length + cumlength) %>% rename_at(vars(-cell_id), ~paste0(., 1))
  prog2 <- task2$progression %>% left_join(task2$milestone_network, by=c("from", "to")) %>% mutate(cumpercentage=percentage*length + cumlength) %>% rename_at(vars(-cell_id), ~paste0(., 2))

  prog <- full_join(prog1, prog2, by=c("cell_id"))

  ggplot(prog) +
    geom_point(aes(cumpercentage1, cumpercentage2)) +
    geom_vline(aes(xintercept=cumlength), data=task1$milestone_network) +
    geom_hline(aes(yintercept=cumlength), data=task2$milestone_network) +
    ggtitle(paste0(task1$id, " -> ", task2$id))
}
