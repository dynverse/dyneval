add_phantom_edges <- function(milestone_ids, milestone_network) {
  bind_rows(
    milestone_network,
    bind_rows(lapply(milestone_ids, function(x) {
      strx <- milestone_network %>%
        filter(from == x)

      if (nrow(strx) > 1) {
        strx <- strx %>%
          mutate(
            angle = seq(0, 120/360*pi*2, length.out = n()),
            x = length * cos(angle),
            y = length * sin(angle)
          )
        poss <- strx %>% select(x, y) %>% as.matrix
        rownames(poss) <- strx$to
        poss %>%
          dist %>%
          as.matrix %>%
          reshape2::melt(varnames = c("from", "to"), value.name = "length") %>%
          mutate(from = as.character(from), to = as.character(to)) %>%
          filter(from != to)
      } else {
        NULL
      }
    })),
    bind_rows(lapply(milestone_ids, function(x) {
      strx <- milestone_network %>%
        filter(to == x)

      if (nrow(strx) > 1) {
        strx <- strx %>%
          mutate(
            angle = seq(0, 120/360*pi*2, length.out = n()),
            x = length * cos(angle),
            y = length * sin(angle)
          )
        poss <- strx %>% select(x, y) %>% as.matrix
        rownames(poss) <- strx$from
        poss %>%
          dist %>%
          as.matrix %>%
          reshape2::melt(varnames = c("from", "to"), value.name = "length") %>%
          mutate(from = as.character(from), to = as.character(to)) %>%
          filter(from != to)
      } else {
        NULL
      }
    })))
}

