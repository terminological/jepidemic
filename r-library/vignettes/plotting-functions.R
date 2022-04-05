## Plotting ----

#' logit scale
#'
#' @description it perform logit scaling with right axis formatting. To not be used directly but with ggplot (e.g. scale_y_continuous(trans = "logit") )
#'
#' @return A scales object
#'
#' @examples
#'
#' library(ggplot2)
#' library(tibble)
#'
#' tibble(pvalue = c(0.001, 0.05, 0.1), fold_change = 1:3) %>%
#'  ggplot(aes(fold_change , pvalue)) +
#'  geom_point() +
#'  scale_y_continuous(trans = "logit")
#'
# TODO: figure out where this best lives #' @export
logit_trans <- function(){
  
  
  if (find.package("functional", quiet = TRUE) %>% length %>% equals(0)) {
    message("Installing functional needed for analyses")
    install.packages("functional", repos = "https://cloud.r-project.org")
  }
  
  trans <- stats::qlogis
  inv <- stats::plogis
  
  trans_new("logit",
            transform = trans,
            inverse = inv,
            breaks = functional::Compose(trans, scales::extended_breaks(), inv),
            format = scales::label_scientific(digits = 2)
  )
}

outputter = function(directory = here::here("output"), datedFile=TRUE, datedSubdirectory=FALSE) {
  directory = fs::path_norm(directory)
  fs::dir_create(directory)
  message("directing output to: ",directory)
  return(function(filename) {
    if(datedSubdirectory) {
      directory = fs::path(directory,Sys.Date())
    }
    ext = fs::path_ext(filename)
    if(datedFile) filename = paste0(fs::path_ext_remove(filename),"-",Sys.Date()) %>% fs::path_ext_set(ext)
    fs::path(directory,filename)
  })
}


## Plots ----
# 
# rtPlot = function(out, ylim=c(0,NA)) {
#   ggplot(out %>% filter(!is.nan(Rt.Mean)),aes(x=Rt.EndDate, y=Rt.Mean, colour=as.factor(Rt.Window)))+
#     geom_point(size=0.5) +
#     geom_errorbarh(aes(x=Rt.EndDate, xmin=Rt.StartDate, xmax=Rt.EndDate), alpha=0.3) +
#     geom_errorbar(aes(y=Rt.Quantile.0.5, ymin = Rt.Quantile.0.025, ymax = Rt.Quantile.0.975), alpha=0.3)+
#     guides(colour="none")+coord_cartesian(ylim=ylim)+
#     geom_hline(yintercept=1,colour="grey50",inherit.aes=FALSE)+
#     xlab("date")+
#     ylab(latex2exp::TeX("$R_t$"))
# }
# 
# rtPanel = function(out, ylim=c(0,5), colourExpr=expr(as.factor(Rt.Window))) {
#   ggplot(out %>% filter(!is.nan(Rt.Mean)), aes(x=Rt.EndDate, y=Rt.Mean, colour=!!colourExpr))+
#     geom_point(size=1) +
#     geom_errorbarh(aes(x=Rt.EndDate, xmin=Rt.StartDate, xmax=Rt.EndDate), alpha=0.3,colour = "grey") +
#     geom_errorbar(aes(y=Rt.Quantile.0.5, ymin = Rt.Quantile.0.025, ymax = Rt.Quantile.0.975), alpha=0.3,colour = "grey")+
#     guides(colour="none")+coord_cartesian(ylim=ylim)
# }
# 
# rtRibbon = function(out, ylim=c(0,NA), mapping = aes(), ...) {
#   defaultAes = aes(x=Rt.EndDate, y=Rt.Mean)
#   dots = rlang::list2(...)
#   combinedAes = modifyList(defaultAes,mapping)
#   ggplot(out %>% filter(!is.nan(Rt.Mean)), mapping = combinedAes)+
#     geom_ribbon(aes(ymin = Rt.Quantile.0.025, ymax = Rt.Quantile.0.975), alpha=0.2,colour = NA,fill="grey20")+
#     geom_ribbon(aes(ymin = Rt.Quantile.0.25, ymax = Rt.Quantile.0.75), alpha=0.3,colour = NA,fill="grey20")+
#     geom_line(...) +
#     geom_hline(yintercept=1,colour="grey50",inherit.aes=FALSE)+
#     guides(colour="none")+coord_cartesian(ylim=ylim)
# }
# 
# incidencePlot = function(out, ilim = c(NA,NA)) {
#   ggplot(out, aes(x=date,y=value,group=subgroup))+
#     geom_point(size=0.5)+
#     geom_line(alpha=0.1)+
#     ylab("cases")
# }
