fig_model_fits_fun_cef <- function(df) {
  ggplot(data = subset(df, EVID == 0 & CMT == 1), aes(x = TIME, y = DV)) + 
    geom_point(size = 0.9, color = "#4a3a27") + 
    geom_line(data = subset(df, EVID != 1 & CMT == 1), aes(x = TIME, y = IPRED), linetype = 1, color = "bisque4") +
    geom_line(data = subset(df, EVID != 1 & CMT == 1), aes(x = TIME, y = PRED), linetype = 2, color = "#9d7b54") + 
    scale_y_continuous(limits = c(0.1, 100),
                       trans = 'log10',
                       breaks = trans_breaks('log10', function(x) 10^x, n=4),
                       labels = trans_format('log10', math_format(10^.x))) +
    ylab(" Ceftolozane concentrations, ng/ml") +
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
    ylab("Tazobactam concentrations, ng/ml") +
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