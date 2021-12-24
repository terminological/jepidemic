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
#' @export
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