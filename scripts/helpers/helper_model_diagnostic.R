
### Functions

### For parameter tables ####  -------------------------------------------------------
# Most CV equations https://ascpt.onlinelibrary.wiley.com/doi/full/10.1002/psp4.12404
getCV_lognormO  <- function(v) sqrt(exp(v) - 1) * 100
#VAR 2 method produces SD values for prop error in the theta block
getCV_propS  <- function(v) (v) * 100
### shouldn't report CV% for add etas unless it's constrained to be positive
# getCV_addO  <- function(v, t) (sqrt(v) / abs(t)) * 100

#' SD and %CV for logit-normal distributions
#' Decided CV% not appropriate for logit transforms so report SD
#' %CV for random variable Y,
#'   Y = a + b * (1 / (1 + exp(-X))),
#' where X ~ N(mu, sigma) is the normally-distributed logit term
#' e.g. for PARAM = 1 / (1 + EXP(-(THETA1 + ETA1)))
#'
#' @param .mean mean of the logit term (THETA1 in the example)
#' @param .var  variance of the logit term (OMEGA(1,1) in the example)
#' @param .a  additive term
#' @param .b  proportional term
getSD_logitO <- function(.mean, .var, .a = 0, .b = 1) {
  # cat("Mean", .mean, " Var: ", .var, "\n")
  sdList = NA
  for (i in 1:length(.mean)) {
    m = .mean[i]
    v = .var[i]
    if (is.na(m) | is.na(v)) {
      sd = NA_real_ } else {
        moments <- logitnorm::momentsLogitnorm(mu = m, sigma = sqrt(v))
        sd <- sqrt(moments[["var"]])
        # cv <- .b * sqrt(moments[["var"]]) / (.b * abs(moments[["mean"]]) + .a) * 100
      }
    sdList = c(sdList, sd)
  }
  return(sdList[-1])
}

# Confidence intervals
lowerCI <- function(est, se) est - 1.96*se
upperCI <- function(est, se) est + 1.96*se
parensSQ <- function(x) paste0('[',x,']')

parensSQ <- function(x) paste0('[',x,']')
parensSQ_CV <- function(.x) glue::glue("[CV\\%=<<.x>>]", .open = "<<", .close  = ">>")
parensSQ_corr <- function(.x) glue::glue("[Corr=<<.x>>]", .open = "<<", .close  = ">>")
parensSQ_se <- function(.x) glue::glue("[SD=<<.x>>]", .open = "<<", .close  = ">>")
getEvenNo = function(x) x[which(x %% 2 == 0)]


# Greek number helper functions
mathMode <- function(.x) glue::glue("$<<.x>>$", .open = "<<", .close  = ">>")
gtGreek <- function(.x) glue::glue("\\<<.x>>", .open = "<<", .close  = ">>")
greekNum <- function(.x, .y) glue::glue("<<.x>>_{<<.y>>}", .open = "<<", .close  = ">>")
expGreek  <- function(.x, .y) glue::glue("$\\exp(\\<<.x>>_{<<.y>>})$", .open = "<<", .close  = ">>")
logitGreek  <- function(.x, .y) glue::glue("$\\exp(\\<<.x>>_{<<.y>>}) / \\newline(1 + \\exp(\\<<.x>>_{<<.y>>}))$", .open = "<<", .close  = ">>")


## check whether "~" is used to signify the associated THETA
checkTransforms <- function(df){
  df$transTHETA = NA
  if(any(str_detect(df$trans, "~"))){
    # if there is a '~' in df$trans, replace NA with the value
    # e.g. logitOmSD ~ THETA1 puts "THETA1" in transTHETA
    df$transTHETA[which(str_detect(df$trans, "~"))] =
      stringr::str_split(df$trans, fixed("~")) %>% map(trimws) %>% map(2) %>% unlist
    
    # Then, remove everything after the "~" in the trans column and
    # replace THETAx in transTHETA with corresponding estimate
    df = df %>%
      mutate(trans = case_when(str_detect(trans, "~") ~
                                 stringr::str_split(trans, fixed("~")) %>% map(trimws) %>% map(1) %>% unlist,
                               TRUE ~ trans),
             ## second, replace THETA with corresponding estimate
             transTHETA = estimate[match(transTHETA, parameter_names)]
      )
  }
  return(df)
}

## Define a series of true/false columns to make filter easier later
defineRows <- function(df){
  df %>%
    mutate(
      TH = stringr::str_detect(name, "TH"),
      OM = stringr::str_detect(name, "OM"),
      S = stringr::str_detect(name, "S"),
      LOG = (trans=="logTrans"),
      LOGIT = (trans=="logitTrans"),
      lognormO = (trans=="lognormalOm"),
      Osd = (trans=="OmSD"),
      logitOsd = (trans=="logitOmSD"),
      propErr = (trans=="propErr"),
      addErr = (trans=="addErr")
    )
}

## calculate 95% confidence intervals
get95CI <- function(df){
  df %>%
    mutate(lower = lowerCI(value, se),
           upper = upperCI(value, se))
}


## calculate % RSE - not used but included if needed
# Note, this is appropriate when parameters are estimated untransformed or in the log
# it may not be appropriate if any other transformations (such as logit) were used
getpRSE <- function(df){
  df %>%
    mutate(pRSE = case_when(fixed ~ "-",
                            # pRSE of a log-trans TH is equivalent to the CV% of a log-trans TH
                            TH & LOG ~ sig  (sqrt(exp(se^2)-1)*100),
                            TH & !LOG & !LOGIT ~ sig ((se/abs(value)) * 100),
                            diag & !LOG & !LOGIT ~ sig ((se/abs(value)) * 100),
                            TRUE ~ "-"))
}

## Back transform parameters estimated in the log domain
# make sure any other calculations, such as CI (and pRSE) are
# done before back-calculating these values
backTrans_log <- function(df){
  df %>%
    mutate(value = case_when(LOG ~ exp(value), TRUE ~ value),
           lower = case_when(LOG ~ exp(lower), TRUE ~ lower),
           upper = case_when(LOG ~ exp(upper), TRUE ~ upper))
}
backTrans_logit <- function(df){
  df %>%
    mutate(value = case_when(LOGIT ~ exp(value)/(1+exp(value)), TRUE ~ value),
           lower = case_when(LOGIT ~ exp(lower)/(1+exp(lower)), TRUE ~ lower),
           upper = case_when(LOGIT ~ exp(upper)/(1+exp(upper)), TRUE ~ upper))
}
## Calculate CV%
getpCV <- function(df){
  df %>%
    mutate(cv = case_when(diag & OM & lognormO ~ sig(getCV_lognormO(value)),
                          #diag & OM & logitOsd ~ sig(getSD_logitO(.mean=transTHETA, .var = value)),
                          S & propErr ~ sig(getCV_propS(value)),
                          TRUE ~ "-"))
}

# value should have estimate [something]
#   theta = estimate only                           # use estimate column
#   omega diagonals = variance [%CV]                # estimate [CV from estimate, stderr]
#   omega off-diagonals = covariance [corr coeff]   # estimate [random_effect_sd]
#   sigma diagonal proportional = variance [%CV]    # estimate [CV from estimate, stderr]
#   sigma diagonal additive = variance [SD]         # estimate [random_effect_sd]
getValueSE <- function(df){
  df %>%
    mutate(value = estimate,
           se = stderr,
           corr_SD = case_when(OM & !diag |
                                 S & diag & addErr ~ sig(random_effect_sd),
                               TRUE ~ "-")
    )
}


# 95% CI should show lower, upper or FIXED
# rounding for display in report
# define what is in estimate column and what is in square brackets
formatValues <- function(df){
  df %>%
    # back transform any parameters here
    backTrans_log() %>%   # back transform from log domain
    backTrans_logit() %>%
    getpCV() %>%          # get % CV
    # format the values for the final table
    mutate(ci = paste0(sig(lower), ', ', sig(upper)),
           ci = if_else(fixed, "FIXED", ci),
           # get sd if needed
           sd = case_when(diag & OM & Osd ~ sig(random_effect_sd),
                          diag & OM & logitOsd ~ sig(getSD_logitO(.mean=transTHETA, .var = value)),
                          TRUE ~ "-"
           ),
           # round values for report table
           value = sig(value),
           # define which values appear where
           value = case_when(diag & OM & Osd |
                               diag & OM & logitOsd ~ glue::glue("{value} {parensSQ_se(sd)}"),
                             
                             diag & OM |
                               diag & S & propErr ~
                               glue::glue("{value} {parensSQ_CV(cv)}"),
                             !diag & OM ~ glue::glue("{value} {parensSQ_corr(corr_SD)}"),
                             diag & S & addErr ~ glue::glue("{value} {parensSQ_se(corr_SD)}"),
                             !diag & S ~ glue::glue("{value} {parensSQ_corr(corr_SD)}"),
                             TRUE ~ value),
           # round shrinkage values for report table
           shrinkage = case_when(is.na(shrinkage) ~ "-",
                                 TRUE ~ sig(shrinkage)))
  
}

## Format the THETA/OMEGA/SIGMA values to display as greek letters with
# subscript numbers
formatGreekNames <- function(df){
  df %>%
    mutate(greekName = name) %>%
    # make column with greek letters and parameter numbers
    separate(greekName,
             into = c("text", "num"),
             sep = "(?<=[A-Za-z])(?=[0-9])"
    ) %>%
    separate(parameter_names,
             into = c("text2", "num2"),
             sep = "A"
    ) %>%
    select(-num, -text2) %>%
    mutate(text = case_when(OM ~ "Omega",
                            S ~ "Sigma",
                            TRUE ~ tolower(text)),
           greek = case_when(TH & LOG ~ expGreek(text, num2),
                             TH & LOGIT ~ logitGreek(text, num2),
                             TRUE ~ mathMode(greekNum(gtGreek(text), num2))
           )
    )
}


## Define which parameters should appear under which panel name in the final table
getPanelName = function(df){
  df %>%
    mutate(type = case_when(panel=="RV" ~ "Residual SD",
                            OM & !diag ~ "Interindividual covariance parameters",
                            OM & diag & panel=="IIV" ~ "Interindividual variance parameters",
                            # IOV not used here but included for convenience
                            OM & diag & panel=="IOV"  ~ "Interoccasion variance parameters",
                            panel=="cov" ~ "Covariate effect parameters",
                            panel=="struct" ~ "Structural model parameters"),
           # Make type a factor and use to sort, this ensures all parameters
           # of the same type are together - needed to make sure panels pull out
           # correct rows
           type_f = case_when(panel=="RV" ~ 6,
                              OM & !diag ~ 5,
                              OM & diag & panel=="IIV" ~ 3,
                              # IOV not used here but included for convenience
                              OM & diag & panel=="IOV"  ~ 4,
                              panel=="cov" ~ 2,
                              panel=="struct" ~ 1)
    ) %>%
    arrange(type_f)
}


##Pooled GOF Plots
#Purpose: post-process nlmixr2 output to derive pooled overlay GOF plots
#Arguments
#fit = nlmixr2 fit object
#strat_var = character, variable to stratify plots by (e.g., "DOSE")
#units_dv = character, units of the dependent variable (e.g., "mg/L")
#units_time = character, units of the time variable (e.g., "hours")
#page_no = integer, page number to print passed to facet_wrap_paginate
#nrow = integer, number of rows of plots per page (Default = 2)
#ncol = integer, number of columns of plots per page (Default = 2)
#log_y = logical, log-scale y-axis (default = FALSE)
create_gof <- function(fit,
                       strat_var,
                       strat_var_val,
                       units_time = "hours",
                       units_dv = "mg/L",
                       cfb = FALSE,
                       obs = TRUE,
                       log_y=FALSE,
                       scales="fixed"){
  
  fit <- fit %>%
    filter(!!as.symbol(strat_var) %in% strat_var_val)
  
  dvu <- ifelse(cfb == FALSE, units_dv, "% Change")
  
  xlab <- paste0("Time (", units_time, ")")
  ylab <- paste0("Concentration (", units_dv, ")")
  
  plot <- ggplot(fit, aes(x = TIME, y=DV, grp = !!as.symbol(strat_var))) +
    stat_summary(fun = mean, geom= "line", linewidth=1.5, color = "navy")+
    stat_summary(aes(x=NTIME, y=PRED),fun = mean, geom= "line",  linewidth=1.5, color = "darkred", inherit.aes = FALSE)+
    stat_summary(aes(x=NTIME, y=IPRED),fun = mean, geom= "line",  linewidth=1.5, color = "darkgreen", inherit.aes = FALSE)+
    facet_wrap(as.formula(paste("~",strat_var)), scales = scales)+
    labs(x=xlab, y=ylab)
  
  if(obs == TRUE) plot <- plot +  geom_point(shape=1, size=3, alpha = 0.5)
  
  if(log_y==TRUE) plot <- plot + scale_y_log10()
  
  if(cfb == TRUE) plot <- plot + geom_hline(yintercept = 0, linewidth = 1, linetype = "dashed")
  
  return(plot)
}





fig_model_fits_fun_cef <- function(df) {
  ggplot(data = subset(df, EVID == 0 & CMT == 1), aes(x = TIME, y = DV)) + 
    geom_point(size = 0.9, color = "#4a3a27") + 
    geom_line(data = subset(df, EVID != 1 & CMT == 1), aes(x = TIME, y = IPRED), linetype = 1, color = "bisque4") +
    geom_line(data = subset(df, EVID != 1 & CMT == 1), aes(x = TIME, y = PRED), linetype = 2, color = "#9d7b54") + 
    scale_y_continuous(limits = c(0.1, 100),
                       trans = 'log10',
                       breaks = trans_breaks('log10', function(x) 10^x, n=4),
                       labels = trans_format('log10', math_format(10^.x))) +
    ylab(" Ceftolozane concentrations, mg/ml") +
    xlab("Time, h") + 
    my_theme +
    facet_wrap(.~ID, nrow = 3)
}

fig_model_fits_fun_taz <- function(df) {
  ggplot(data = subset(df, EVID == 0 & CMT == 3), aes(x = TIME, y = DV)) + 
    geom_point(size = 0.9, color = "#00253d") + 
    geom_line(data = subset(df, EVID != 1 & CMT == 3), aes(x = TIME, y = IPRED), linetype = 1, color = "deepskyblue4") +
    geom_line(data = subset(df, EVID != 1 & CMT == 3), aes(x = TIME, y = PRED), linetype = 2, color = "#00426f") + 
    scale_y_continuous(limits = c(0.1, 100),
                       trans = 'log10',
                       breaks = trans_breaks('log10', function(x) 10^x, n=4),
                       labels = trans_format('log10', math_format(10^.x))) +
    ylab("Tazobactam concentrations, mg/ml") +
    xlab("Time, h") + 
    my_theme +
    facet_wrap(.~ID, nrow = 3)
}

fig_gof_fun <- function(df) {
  
  p1 <- ggplot(subset(df, MDV == 0), aes(x = IPRED, y = DV)) +
    geom_point(data = subset(df, EVID == 0 & CMT == 1), color = "bisque4", size = 0.95) +
    geom_point(data = subset(df, EVID == 0 & CMT == 3), color = "deepskyblue4", size = 0.95) +
    geom_abline(intercept = 0, slope = 1, linetype = 2) +
    labs(x = "IPRED", y = "DV") + 
    scale_y_continuous(limits = c(0.1, 100),
                       trans ='log10',
                       breaks = trans_breaks('log10', function(x) 10^x, n=4),
                       labels = trans_format('log10', math_format(10^.x))) + 
    scale_x_continuous(limits = c(0.1, 100),
                       trans = 'log10',
                       breaks = trans_breaks('log10', function(x) 10^x, n=4),
                       labels = trans_format('log10', math_format(10^.x)))
  
  p2 <- ggplot(subset(df, MDV == 0), aes(x = PRED, y = DV)) +
    geom_point(data = subset(df, EVID == 0 & CMT == 1), color = "bisque4", size = 0.95) +
    geom_point(data = subset(df, EVID == 0 & CMT == 3), color = "deepskyblue4", size = 0.95) +
    geom_abline(intercept = 0, slope = 1, linetype = 2) +
    labs(x = "PRED", y = "DV") + 
    scale_y_continuous(limits = c(0.1, 100),
                       trans = 'log10',
                       breaks = trans_breaks('log10', function(x) 10^x, n=4),
                       labels = trans_format('log10', math_format(10^.x))) + 
    scale_x_continuous(limits = c(0.1, 100),
                       trans = 'log10',
                       breaks = trans_breaks('log10', function(x) 10^x, n=4),
                       labels = trans_format('log10', math_format(10^.x)))
  
  p3 <- ggplot(subset(df, MDV == 0), aes(x = TIME, y = CWRES)) +
    geom_point(data = subset(df, EVID == 0 & CMT == 1), color = "bisque4", size = 0.95) +
    geom_point(data = subset(df, EVID == 0 & CMT == 3), color = "deepskyblue4", size = 0.95) +
    geom_hline(yintercept = 0, linetype = 2) + 
    labs(x ="Time (h)", y = "CWRES")
  
  p4 <- ggplot(subset(df, MDV == 0), aes(x = PRED, y = CWRES)) +
    geom_point(data = subset(df, EVID == 0 & CMT == 1), color = "bisque4", size = 0.95) +
    geom_point(data = subset(df, EVID == 0 & CMT == 3), color = "deepskyblue4", size = 0.95) +
    geom_hline(yintercept = 0, linetype = 2) +
    labs(x = "PRED", y = "CWRES") + 
    scale_x_continuous(limits = c(0.1, 100),
                       trans = 'log10',
                       breaks = trans_breaks('log10', function(x) 10^x, n=4),
                       labels = trans_format('log10', math_format(10^.x)))
  
  fig_gof <- (p1+p2)/(p3+p4)
  
  return(fig_gof)
}