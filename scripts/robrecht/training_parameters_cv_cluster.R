library(cowplot)
library(tidyverse)
library(dyneval)
library(mlrMBO)
library(parallelMap)
library(PRISM)

output_root_folder <- "results/output_dyngen_paramtraincv/"
dir.create(output_root_folder, recursive = T)

## load datasets
.datasets_location = "../dyngen/results/4/" # needs to be defined, to let dyngen know where the datasets are
# tasks <- load_datasets(8) # this function takes way too long due to the geodesic distances being calculated
# saveRDS(tasks, paste0(.datasets_location, "tasks.rds"))
tasks <- readRDS(paste0(.datasets_location, "tasks.rds")) %>%
  mutate(group = paste0(platform_id, "_", takesetting_type)) %>%
  group_by(group, ti_type) %>% mutate(subtask_ix = seq_len(n())) %>%
  ungroup()

## select datasets # limit for now
select_tasks <- tasks %>% filter(platform_id == "fluidigm_c1", takesetting_type == "snapshot") %>% arrange(ti_type, subtask_ix)

# limit even further
select_tasks <- select_tasks %>% filter(subtask_ix %in% c(1,2), ti_type == "consecutive_bifurcating")

dyneval:::benchmark_suite_submit(select_tasks, select_tasks$group, select_tasks$subtask_ix)




# load(paste0(output_root_folder, "temp.RData"))
#
# qsub_handle$stop_on_error <- F
# qsub_handle$verbose <- F
#
# post_fun <- function(rds_i, out_rds) {
#   grid_i <- rds_i
#   fold_i <- grid$fold_i[[rds_i]]
#   repeat_i <- grid$repeat_i[[rds_i]]
#   method_i <- grid$method_i[[rds_i]]
#   design <- out_rds$design
#   train_out <- out_rds$tune_train
#   test_out <- out_rds$tune_test
#
#   eval_summ_gath <- bind_rows(
#     data.frame(type = "train", train_out$opt.path$env$path) %>% mutate(
#       grid_i, repeat_i, fold_i, method_i,
#       param_i = seq_len(n()),
#       time = train_out$opt.path$env$exec.time
#     ),
#     data.frame(type = "test", test_out$opt.path$env$path) %>% mutate(
#       grid_i, repeat_i, fold_i, method_i,
#       param_i = seq_len(n()),
#       time = test_out$opt.path$env$exec.time
#     )
#   ) %>% filter(param_i <= nrow(train_out$opt.path$env$path)) %>%
#     dplyr::as_data_frame()
#
#   if (!all(eval_summ_gath$y_1 == -1)) {
#     eval_summ <- eval_summ_gath %>%
#       gather(eval_metric, score, y_1:y_3, time) %>%
#       mutate(comb = paste0(type, "_", eval_metric)) %>%
#       select(-type, -eval_metric) %>%
#       spread(comb, score) %>%
#       mutate(iteration_i = train_out$opt.path$env$dob[param_i]) %>%
#       arrange(param_i)
#
#     ## collect the scores per dataset individually
#     eval_ind <- bind_rows(lapply(seq_len(nrow(eval_summ)), function(param_i) {
#       iteration_i <- eval_summ$iteration_i[[param_i]]
#       bind_rows(
#         if (eval_summ$train_y_1[[param_i]] >= 0) { # did this execution finish correctly?
#           train_out$opt.path$env$extra[[param_i]]$.summary %>% mutate(grid_i, repeat_i, fold_i, method_i, param_i, iteration_i, fold_type = "train")
#         } else {
#           NULL
#         },
#         if (eval_summ$test_y_1[[param_i]] >= 0) { # did this execution finish correctly?
#           test_out$opt.path$env$extra[[param_i]]$.summary %>% mutate(grid_i, repeat_i, fold_i, method_i, param_i, iteration_i, fold_type = "test")
#         } else {
#           NULL
#         }
#       )
#     })) %>% left_join(tasks %>% dplyr::select(type, ti_type, name, experimentid, platformname, version, dataset_i, subdataset_i), by = c("task_name" = "name")) %>%
#       as_data_frame
#
#     ## group them together per ti_type
#     eval_grp <- eval_ind %>% group_by(ti_type, grid_i, repeat_i, fold_i, method_i, iteration_i, param_i, fold_type) %>% summarise_at(metrics, mean) %>% ungroup()
#   } else {
#     eval_summ <- NULL
#     eval_summ_gath <- NULL
#     eval_ind <- NULL
#     eval_grp <- NULL
#   }
#
#   dplyr::lst(eval_summ, eval_summ_gath, eval_ind, eval_grp)
# }
#
# output <- qsub_retrieve(qsub_handle, post_fun = post_fun, wait = F)
#
# # save(output, file = paste0(output_root_folder, "temp_output.RData"))
# load(paste0(output_root_folder, "temp_output.RData"))
#
# method_names <- methods %>% map_chr(~ .$name)
# grid %>% as_data_frame %>% mutate(method_str = method_names[method_i])
#
# output <- output[!sapply(output, function(x) length(x) == 1 && is.na(x))]
#
# scores <- output %>% map_df(~ .$eval_summ) %>% mutate(grid_str = paste0("Repeat ", repeat_i, ", fold ", fold_i), method_str = method_names[method_i])
# individual_scores <- output %>% map_df(~ .$eval_ind) %>% mutate(grid_str = paste0("Repeat ", repeat_i, ", fold ", fold_i), method_str = method_names[method_i])
# grouped_scores <- output %>% map_df(~ .$eval_grp) %>% mutate(grid_str = paste0("Repeat ", repeat_i, ", fold ", fold_i), method_str = method_names[method_i])
#
# scores %>% group_by(repeat_i, fold_i, method_i) %>% slice(1) %>% ungroup() %>% group_by(repeat_i, fold_i) %>% summarise(n = n()) %>% ungroup()
# scores %>% group_by(repeat_i, fold_i, method_str) %>% slice(1) %>% ungroup() %>% group_by(method_str) %>% summarise(n = n()) %>% ungroup()
# scores %>% group_by(repeat_i, fold_i, method_str) %>% slice(1) %>% ungroup() %>% group_by(method_str, fold_i) %>% summarise(n = n()) %>% ungroup() %>% spread(fold_i, n, fill = 0)
#
# best_scores <- scores %>% group_by(repeat_i, fold_i, method_i) %>% arrange(desc(train_y_1)) %>% slice(1) %>% ungroup() %>%
#   arrange(desc(train_y_1)) %>% mutate(method_rank = paste0("#", seq_len(n()), ": ", method_str), method_fac = factor(method_rank, levels = method_rank))
# best_param <- best_scores %>% select(repeat_i, fold_i, method_i, param_i, method_rank, method_fac)
# best_individual_scores <- individual_scores %>% right_join(best_param, by = c("repeat_i", "fold_i", "method_i", "param_i"))
# best_grouped_scores <- grouped_scores %>% right_join(best_param, by = c("repeat_i", "fold_i", "method_i", "param_i"))
#
# ggplot(scores %>% filter(train_y_1 >= 0, test_y_1 >= 0)) +
#   geom_vline(aes(xintercept = train_y_1), scores %>% filter(param_i == 1, train_y_1 >= 0, test_y_1 >= 0)) +
#   geom_hline(aes(yintercept = test_y_1), scores %>% filter(param_i == 1, train_y_1 >= 0, test_y_1 >= 0)) +
#   geom_vline(aes(xintercept = train_y_1), best_scores %>% filter(train_y_1 >= 0, test_y_1 >= 0), colour = "red") +
#   geom_hline(aes(yintercept = test_y_1), best_scores %>% filter(train_y_1 >= 0, test_y_1 >= 0), colour = "red") +
#   geom_point(aes(train_y_1, test_y_1, colour = iteration_i)) +
#   scale_colour_distiller(palette = "RdBu") +
#   facet_grid(method_str~grid_str) +
#   coord_equal()
#
# aggregated_scores <- best_scores %>% group_by(method_i) %>% mutate_at(c(paste0("train_y_", 1:3), paste0("test_y_", 1:3)), mean) %>% summarise_all(head1)
# aggregated_scores_spr <- aggregated_scores %>% gather(eval_metric, score, contains("_y_"))
# #ggplot(aggregated_scores) + geom_point(aes(train_y_1, test_y_1, colour = method_str))
# ggplot(aggregated_scores_spr %>% filter(eval_metric %in% c("train_y_1", "test_y_1")) %>% mutate(eval_metric = factor(eval_metric, levels = c("train_y_1", "test_y_1")))) +
#   geom_bar(aes(method_str, score, fill = eval_metric), stat = "identity", position = "dodge")
#
# grspr <- grouped_scores %>% select(-Q_local, -correlation) %>% spread(fold_type, Q_global)
# bgrspr <- best_grouped_scores %>% select(-Q_local, -correlation) %>% spread(fold_type, Q_global)
# ggplot(grspr) +
#   geom_vline(aes(xintercept = train), grspr %>% filter(param_i == 1)) +
#   geom_hline(aes(yintercept = test), grspr %>% filter(param_i == 1)) +
#   geom_vline(aes(xintercept = train), bgrspr, colour = "red") +
#   geom_hline(aes(yintercept = test), bgrspr, colour = "red") +
#   geom_point(aes(train, test, colour = iteration_i)) +
#   scale_colour_distiller(palette = "RdBu") +
#   facet_grid(ti_type~method_str) +
#   coord_equal()
#
# ggplot(best_grouped_scores) + geom_bar(aes(method_fac, Q_global, fill = method_fac), stat = "identity", position = "dodge") + facet_grid(fold_type~ti_type)

#
# pdf(paste0(output_root_folder, "1_celltree_overallscore.pdf"), 10, 4)
# # png(paste0(output_root_folder, "1_celltree_overallscore.png"), 1000, 400)
# ggplot(scores %>% filter(y_train >= 0, y_test >= 0)) +
#   geom_vline(aes(xintercept = y_train), scores %>% filter(param_i == 1)) +
#   geom_hline(aes(yintercept = y_test), scores %>% filter(param_i == 1)) +
#   geom_vline(aes(xintercept = y_train), best_scores, colour = "red") +
#   geom_hline(aes(yintercept = y_test), best_scores, colour = "red") +
#   geom_point(aes(y_train, y_test, colour = iteration_i)) +
#   scale_colour_distiller(palette = "RdBu") +
#   facet_wrap(~grid_str, nrow = 1) +
#   coord_equal()
# dev.off()
#
# pdf(paste0(output_root_folder, "2_celltree_groupedscore.pdf"), 16, 12)
# # png(paste0(output_root_folder, "2_celltree_groupedscore.png"), 1600, 1200)
# ggplot(grouped_scores %>% select(-auc_lcmc) %>% spread(fold_type, max_lcmc)) +
#   geom_vline(aes(xintercept = train), grouped_scores %>% select(-auc_lcmc) %>% spread(fold_type, max_lcmc) %>% filter(param_i == 1)) +
#   geom_hline(aes(yintercept = test), grouped_scores %>% select(-auc_lcmc) %>% spread(fold_type, max_lcmc) %>% filter(param_i == 1)) +
#   geom_vline(aes(xintercept = train), best_grouped_scores %>% select(-auc_lcmc) %>% spread(fold_type, max_lcmc), colour = "red") +
#   geom_hline(aes(yintercept = test), best_grouped_scores %>% select(-auc_lcmc) %>% spread(fold_type, max_lcmc), colour = "red") +
#   geom_point(aes(train, test, colour = iteration_i)) +
#   scale_colour_distiller(palette = "RdBu") +
#   facet_grid(grid_str~ti_type) +
#   coord_equal()
# dev.off()
#
#
#
# # ggplot(best_individual_scores) + geom_bar(aes(factor(i), max_lcmc, fill = ti_type), stat = "identity") + facet_wrap(~ti_type)
# # png("rplot3.png", 1000, 300)
# # ggplot(best_grouped_scores) + geom_bar(aes(ti_type, max_lcmc, fill = factor(fold_type, levels = c("train", "test"))), stat = "identity", position = "dodge") + labs(fill = "Fold type")
# # dev.off()
# ggplot() +
#   geom_bar(aes(ti_type, max_lcmc, group = factor(fold_type, levels = c("train", "test")), colour = "best"), fill = NA, stat = "identity", position = "dodge", best_grouped_scores) +
#   geom_boxplot(aes(ti_type, max_lcmc, fill = factor(fold_type, levels = c("train", "test"))), position = "dodge", grouped_scores) +
#   labs(fill = "Fold type") +
#   scale_fill_brewer(palette = "Set2") +
#   scale_colour_brewer(palette = "Set1")
#
#
# g <- ggplot() +
#   geom_bar(aes(ti_type, max_lcmc, group = factor(fold_type, levels = c("train", "test")), colour = "best"), fill = NA, stat = "identity", position = "dodge", best_grouped_scores) +
#   geom_boxplot(aes(ti_type, max_lcmc, fill = factor(fold_type, levels = c("train", "test"))), position = "dodge", grouped_scores) +
#   labs(fill = "Fold type") +
#   scale_fill_brewer(palette = "Set2") +
#   scale_colour_brewer(palette = "Set1")
# g
#
# png("rplot4.png", 1000, 600)
# g
# dev.off()
#
# ggplot(individual_scores) + geom_point(aes(run_i, max_lcmc, colour = fold_type)) + facet_grid(i ~ ti_type)
#
# # ggplot(grouped_scores) + geom_point(aes(point, max_lcmc, colour = ti_type)) + scale_colour_brewer(palette = "Dark2")
# # ggplot(grouped_scores) + geom_point(aes(point, max_lcmc, colour = ti_type)) + scale_colour_brewer(palette = "Dark2") + facet_wrap(~ti_type, nrow = 1) + geom_vline(xintercept = which(diff(scores$iteration) != 0)+.5)
# # ggplot(grouped_scores) + geom_point(aes(point, max_lcmc, colour = ti_type)) + scale_colour_brewer(palette = "Dark2") + facet_wrap(~ti_type, ncol = 1) + geom_vline(xintercept = which(diff(scores$iteration) != 0)+.5)
# #
# ggplot(scores) + geom_point(aes(i, y_train)) + geom_vline(xintercept = which(diff(scores$iteration) != 0)+.5)
# ggplot(scores) + geom_point(aes(i, y_test)) + geom_vline(xintercept = which(diff(scores$iteration) != 0)+.5)
# ggplot(scores) + geom_point(aes(y_train, y_test))
# ggplot(scores %>% filter(y_train > 0)) + geom_point(aes(y_train, y_test, colour = iteration)) + scale_colour_distiller(palette = "RdBu")
# ggplot(grouped_scores %>% select(-auc_lcmc) %>% spread(fold_type, max_lcmc)) + geom_point(aes(train, test, colour = iteration)) + scale_colour_distiller(palette = "RdBu") + facet_wrap(~ti_type)
# # ggplot(scores) + geom_point(aes(i, y, colour = num_topics_lower)) + geom_vline(xintercept = which(diff(scores$iteration) != 0)+.5) + scale_colour_distiller(palette="RdBu")
# # ggplot(scores) + geom_point(aes(i, y, colour = num_topics_upper)) + geom_vline(xintercept = which(diff(scores$iteration) != 0)+.5) + scale_colour_distiller(palette="RdBu")
# # ggplot(scores) + geom_point(aes(i, y, colour = sd_filter)) + geom_vline(xintercept = which(diff(scores$iteration) != 0)+.5) + scale_colour_distiller(palette="RdBu")
# # ggplot(scores) + geom_point(aes(i, y, colour = tot_iter)) + geom_vline(xintercept = which(diff(scores$iteration) != 0)+.5) + scale_colour_distiller(palette="RdBu")
# # ggplot(scores) + geom_point(aes(i, y, colour = tolerance)) + geom_vline(xintercept = which(diff(scores$iteration) != 0)+.5) + scale_colour_distiller(palette="RdBu")
# # ggplot(scores) + geom_point(aes(i, y, colour = width_scale_factor)) + geom_vline(xintercept = which(diff(scores$iteration) != 0)+.5) + scale_colour_distiller(palette="RdBu")
# # ggplot(scores) + geom_point(aes(i, y, colour = outlier_tolerance_factor)) + geom_vline(xintercept = which(diff(scores$iteration) != 0)+.5) + scale_colour_distiller(palette="RdBu")
#
# ggplot(scores) + geom_point(aes(i, y, colour = "train")) + geom_point(aes(i, y_test, colour = "test")) + geom_vline(xintercept = which(diff(scores$iteration) != 0)+.5)
# ggplot(grouped_scores_test) + geom_point(aes(point, max_lcmc, colour = ti_type)) + scale_colour_brewer(palette = "Dark2") + facet_wrap(~ti_type, ncol = 1) + geom_vline(xintercept = which(diff(scores$iteration) != 0)+.5)
#
# ggplot(grouped_scores) + geom_point(aes(point, max_lcmc, colour = fold_type)) + scale_colour_brewer(palette = "Dark2") + facet_wrap(~ti_type, ncol = 1) + geom_vline(xintercept = which(diff(scores$iteration) != 0)+.5)
# ggplot(grouped_scores) + geom_point(aes(point, max_lcmc, colour = fold_type)) + scale_colour_brewer(palette = "Dark2") + facet_grid(ti_type~fold_type) + geom_vline(xintercept = which(diff(scores$iteration) != 0)+.5)
