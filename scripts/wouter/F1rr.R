percs <- task$state_percentages %>% select(-id) %>% as.matrix()
task$state_network


branch_assignment <- apply((task$state_percentages[, unique(task$state_network$from)] > 0), 1, function(x) which(x)[[1]])

percentages = pred_output$state_percentages %>% select(-id)
network = pred_output$state_network

####

percentages = task$state_percentages %>% select(-id)
network = task$state_network
cal_branch_assignment = function(percentages, network) {
  branch_assignment <- apply(percentages, 1, function(x) {
    positives = names(which(x > 0))
    if(length(intersect(positives, network$from)) > 0) {
      intersect(positives, network$from)[[1]]
    } else {
      network %>% filter(to==positives[[1]]) %>% .$from %>% .[[1]]
    }
  }) %>% factor() %>% as.numeric()
}


branch_assignment_true = cal_branch_assignment(task$state_percentages %>% select(-id), task$state_network)
branch_assignment_observed = cal_branch_assignment(pred_output$state_percentages %>% select(-id), pred_output$state_network)

F1rr = function(labels1, labels2) {
  overlaps = map(unique(labels1), function(truelabel) {
    map_dbl(unique(labels2), function(observedlabel) {
      sum((labels1 == truelabel) & (labels2 == observedlabel)) / sum((labels1 == truelabel) | (labels2 == observedlabel))
    })
  }) %>% invoke(cbind, .)
  print(overlaps)
  2/(1/mean(apply(overlaps, 1, max)) + 1/mean(apply(overlaps, 2, max)))
}

F1rr(branch_assignment_true, branch_assignment_observed)

score_prediction_F1rr = function(task, pred_output) {
  branch_assignment_true = cal_branch_assignment(task$state_percentages %>% select(-id), task$state_network)
  branch_assignment_observed = cal_branch_assignment(pred_output$state_percentages %>% select(-id), pred_output$state_network)
  F1rr(branch_assignment_true, branch_assignment_observed)
}
map(pred_outputs, ~score_prediction_F1rr(task, .))
