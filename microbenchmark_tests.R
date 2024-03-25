
func_input <- gen_func_input()
design_matrices <- func_input$design_matrices
semi_def_matrices <- func_input$semi_def_matrices
outcome_list <- func_input$outcome_list
summary_model <- func_input$summary_model


find_mle_parameters(c(1,1,1), design_matrices, semi_def_matrices,outcome_list)
summary_model

microbenchmark(find_mle_parameters(c(1,1,1), design_matrices, semi_def_matrices,outcome_list), lme4::lmer(y ~ 1 + (1|klasse), data=DF, REML=FALSE))

profvis(find_mle_parameters(c(1,1,1), design_matrices, semi_def_matrices,outcome_list))


mean_mbm_slow <- mean(c(7.185012,8.090789,7.88007,8.555246,6.854877))
