#------------------------------------------------------------------------------
fun_cat_cov<-function(cov,etaM){
  p1<-ggplot(etaM, aes_string(x=gsub("//.*","",cov), y="medianInd")) + 
    geom_boxplot()+
    geom_jitter()+
    geom_hline(yintercept = 0, color="gray", linetype="dotted", size=1) + 
    labs(y="Etas",x=gsub(".*//","",cov))+
    facet_wrap(~type, scales = "free_y")
  return(p1)
}

fun_cont_cov<-function(cov,etaM){
  p1<-ggplot(etaM, aes_string(x=gsub("//.*","",cov), y="medianInd")) + 
    geom_point()+
    geom_errorbar(aes(ymin=lbInd, ymax=ubInd), width=.2, position=position_dodge(0.05))+
    geom_hline(yintercept = 0, color="gray", linetype="dotted", size=1) + 
    labs(y="Etas",x=gsub(".*//","",cov))+
    facet_wrap(~type, scales = "free_y")
  return(p1)
}

#------------------------------------------------------------------------------
extract_etas_fun <- function(mod_to_plot,xdata){
  idata = xdata %>%
    distinct(ID,.keep_all = T)%>%
    select(ID,AGE,BW, SEX) %>% 
    mutate(i = 1:n()) %>%
    mutate(SEX = if_else(SEX==1,"Male", "Female"))
  
  etaM <-  posterior::as_draws_df(mod_to_plot$draws("etaM")) %>% 
    tidybayes::spread_draws(etaM[i,type])%>%
    mutate(etaM = log(etaM)) %>%
    group_by(i,type) %>%
    summarize(lbInd = quantile(etaM, probs = 0.05, na.rm = TRUE),
              medianInd = quantile(etaM, probs = 0.5, na.rm = TRUE),
              meanInd = mean(etaM, na.rm = TRUE),
              ubInd = quantile(etaM, probs = 0.95, na.rm = TRUE)) %>%
    mutate(type = case_when(
      type == "1" ~ "CLC",
      type == "2" ~ "QC",
      type == "3" ~ "V1C",
      type == "4" ~ "V2C",
      type == "5" ~ "CLT",
      type == "6" ~ "QT",
      type == "7" ~ "V1T", 
      type == "8" ~ "V2T"))%>%
    left_join(idata)
}