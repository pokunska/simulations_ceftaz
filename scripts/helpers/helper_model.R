##Full Parameter Table
#Purpose: Derive parameter table values
#Arguments:
#tab - data frame
#keep fixed - logical - should fixed parameters be retained in the table (default = FALSE)
full_parameter_table <- function(tab, keepfixed=FALSE, ...){
  if(keepfixed == FALSE) tab <- tab %>% filter(fixed == FALSE)
  tab <- tab %>%
    mutate(rse = ifelse(fixed == TRUE, NA_real_, abs(stderr/estimate*100)),
           ci95_lb = signif(estimate - stderr*1.96, digits = 3),
           ci95_ub = signif(estimate + stderr*1.96, digits = 3),
           trans_estimate = case_when(is.na(param_type) ~ NA_real_,
                                      param_type == "OMEGA" & fixed == TRUE ~ NA_real_,
                                      param_type == "THETA" & grepl("LOGD", label, ignore.case = TRUE) ~ exp(estimate),
                                      param_type == "THETA" & grepl("LOGITD", label, ignore.case = TRUE) ~ exp(estimate)/(exp(estimate)+1),
                                      param_type == "THETA" & unit == "TS" & label == "LOG ADD" ~ estimate*100,
                                      param_type == "THETA" & unit == "TS" & label == "ADD" ~ estimate,
                                      param_type == "THETA" & unit == "TS" & label == "PROP" ~ estimate*100,
                                      param_type == "OMEGA" & diag == FALSE ~ random_effect_sd,
                                      param_type == "OMEGA" ~ sqrt(exp(estimate)-1)*100,
                                      param_type == "SIGMA" ~ sqrt(estimate)*100,
                                      TRUE ~ NA_real_),
           trans_ci95_lb = case_when(is.na(param_type) ~ NA_real_,
                                     fixed == TRUE ~ NA_real_,
                                     param_type == "THETA" & grepl("LOGD", label, ignore.case = TRUE) ~ exp(ci95_lb),
                                     param_type == "THETA" & grepl("LOGITD", label, ignore.case = TRUE) ~ exp(ci95_lb)/(exp(ci95_lb)+1),
                                     param_type == "THETA" & unit == "TS" & label == "LOG ADD" ~ ci95_lb*100,
                                     param_type == "THETA" & unit == "TS" & label == "ADD" ~ ci95_lb,
                                     param_type == "THETA" & unit == "TS" & label == "PROP" ~ ci95_lb*100,
                                     param_type == "OMEGA" & diag == FALSE ~  NA_real_,
                                     param_type == "OMEGA" & ci95_lb < 0 ~ 0,
                                     param_type == "OMEGA" ~ sqrt(exp(ci95_lb)-1)*100,
                                     param_type == "SIGMA" ~ sqrt(ci95_lb)*100,
                                     TRUE ~ NA_real_),
           trans_ci95_ub = case_when(is.na(param_type) ~ NA_real_,
                                     fixed == TRUE ~ NA_real_,
                                     param_type == "THETA" & grepl("LOGD", label, ignore.case = TRUE) ~ exp(ci95_ub),
                                     param_type == "THETA" & grepl("LOGITD", label, ignore.case = TRUE) ~ exp(ci95_ub)/(exp(ci95_ub)+1),
                                     param_type == "THETA" & unit == "TS" & label == "LOG ADD" ~ ci95_ub*100,
                                     param_type == "THETA" & unit == "TS" & label == "ADD" ~ ci95_ub,
                                     param_type == "THETA" & unit == "TS" & label == "PROP" ~ ci95_ub*100,
                                     param_type == "OMEGA" & diag == FALSE ~  NA_real_,
                                     param_type == "OMEGA"~ sqrt(exp(ci95_ub)-1)*100,
                                     param_type == "SIGMA" ~ sqrt(ci95_ub)*100,
                                     TRUE ~ NA_real_))
}

##Formatted Parameter Table
#Purpose: Derive a formatted parameter table for output and reporting
#Arguments:
#tab - data frame - full_parameter_table output data frame
#sum - list - bbi model_summary() output list
format_parameter_table <- function(tab, sum, model_dir, run, ...){
  ofv <- summary_log(model_dir)$ofv[summary_log(model_dir)$run == run]
  cn <- sum$condition_number[[1]]$condition_number
  nsubj <- sum$run_details$number_of_subjects
  nobs <- sum$run_details$number_of_obs
  nrecords <- sum$run_details$number_of_data_records
  min_term <- sum$run_heuristics$minimization_terminated
  cov_abort <- sum$run_heuristics$covariance_step_aborted
  bound <- sum$run_heuristics$parameter_near_boundary
  zero_grad <- sum$run_heuristics$has_final_zero_gradient
  
  run_details <- data.frame(
    "parameter_names" = c("ofv", "condition_number", "number_of_subjects", "number_of_obs", "number_of_records",
                          "minimization_terminated", "covariance_step_aborted", "parameter_near_boundary", "has_final_zero_gradient"),
    "estimate" = c(ofv, cn, nsubj, nobs, nrecords,
                   min_term, cov_abort, bound, zero_grad)
  )
  
  tab <- tab %>%
    mutate(trans_estimate = signif(trans_estimate, digits = 3),
           ci95 = paste0("(", signif(ci95_lb, digits=3), "; " ,signif(ci95_ub, digits=3), ")"),
           trans_ci95 = ifelse(is.na(trans_estimate), NA_character_,
                               paste0("(", signif(trans_ci95_lb, digits=3), "; " ,signif(trans_ci95_ub, digits=3), ")")),
           unit = case_when(type == "[P]"~ "% CV",
                            type == "[A]"~ "SD",
                            TRUE ~ unit)) %>%
    select(parameter_names, param_type, label, unit, type, estimate, stderr, rse,ci95, trans_estimate, trans_ci95, shrinkage) %>%
    bind_rows(run_details)
  
  return(tab)
}