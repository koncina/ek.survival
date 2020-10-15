NULL

#' Fit Proportional Hazards Regression Model for the survival data
#'
#' @details
#' Define a standalone function to pretty print the model in the kaplan_scan object
#'
#' @importFrom survival coxph Surv
#'
#' @param time Follow up time
#' @param event The status indicator, normally 0=alive, 1=dead.
#' @param expression The expression stratus (normally 'low' vs 'high')
#'
#' @return an object of class coxph representing the fit. See coxph.object for details.
ks_coxph <- function(time, event, expression) {
  coxph(formula = Surv(time = time, event = event) ~ expression)
}

#' @importFrom crayon make_style finish
#' @importFrom cli cat_line
#' @importFrom stats start
#' @export
print.kaplan_scan <- function(x, ...) {
  `%+%` <- crayon::`%+%`
  darkgey <- make_style("darkgrey")
  orange <- make_style("orange")
  percent_low <- round(100 * x[["cutoff_rank"]] / nrow(x[["data"]]), 2)
  percent_high <- 100 - percent_low
  cat_line("Best cutoff is: " %+% orange(percent_low) %+% orange("% low") %+% " and " %+% orange(percent_high)  %+% orange("% high"))
  cat_line(paste0(x[["p_value_adj_method"]], " adjusted p-value: ", round(x[["p_value_adj"]], 3), " (unadjusted: ", round(x[["p_value"]], 3), ")"))
  cat(start(darkgey))
  cat_line("Model:")
  print(x[["model"]])
  cat(finish(darkgey))
}

#' Augment data with information from a kaplan_scan object
#'
#' @importFrom broom augment
#' @importFrom dplyr mutate
#'
#' @param x A kaplan_scan object
#' @param ... Extra arguments, not used
#'
#' @return A tibble with 4 columns: time, event, expression and expression_level.
#'
#' @method augment kaplan_scan
#' @export
augment.kaplan_scan <- function(x) {
  mutate(x[["data"]], expression_level = x[["expression_level"]])
}

#' Stratify participants using the KaplanScan algorithm.
#'
#' Use the KaplanScan algorithm from the \href{http://r2.amc.nl/}{R2 Genomics Analysis and Visualisation Platform}
#'
#'
#' @details
#' Implemented as described in the \href{http://hgserver1.amc.nl/r2/help/manual/R2-Manual.html#section_5.2.5}{R2 manual}: \emph{To this end, we start by ordering the samples within a dataset from the lowest expressor to the highest expressor for a given gene. In the next step we start at the lower end of expression and determine the logrank test for the first 8 samples against all the higher expressed samples. We then increase the group by the addition of the 9th sample to the first group and repeat this procedure until the second group contains 8 samples. In the final step the best logrank test is chosen to represent the optimal cut-off point for the current gene. As we are performing a whole series of tests to find the most optimal cut-off, all the p-values of the logrank tests are corrected for multiple testing by bonferoni.}
#' The procedure describesd above can be performed iteratively to all genes on the array, which would thereby identify the most significant logrank tests. This is exactly what the kaplan scan for a group of genes does.
#'
#' @importFrom broom tidy
#' @importFrom stats p.adjust
#' @importFrom dplyr select arrange
#' @importFrom purrr set_names map map_dbl
#' @importFrom rlang abort enquos
#' @importFrom tibble rowid_to_column
#'
#' @param data A data frame.
#' @param time,event,expression Names of the columns containing the follow up time, status indicator and gene expression.
#' @param p.adjust.method p value correction method (see \link[stats]{p.adjust}).
#'
#' @return \code{kaplan_scan} returns an object of \link[base]{class} "kaplan_scan"
#'
#' @export
kaplan_scan <- function(data, time, event, expression, p.adjust.method = "bonferroni") {

  n_rows <- nrow(data)
  if (n_rows < 16) rlang::abort("nrow(data) must be greater than 16")

  # select variables of interest
  ks_vars <- enquos(time = time, event = event, expression = expression)
  ks_data <- select(data, !!!ks_vars)
  if (!identical(names(ks_data), c("time", "event", "expression")))
    rlang::abort("Provide the time, event and expression columns")
  # sort by increasing expression
  ks_data <- rowid_to_column(ks_data, "row")
  ks_data <- arrange(ks_data, expression)

  # Split expression in low vs high starting from rank 8 up to n - 8:
  high_low_split <- map(set_names(seq(from = 8, length.out = n_rows - 15)), seq, from = 1)
  expression_level <- map(high_low_split, replace, x = rep("high", n_rows), values = "low")

  ks_models <- map(expression_level, ks_coxph, time = ks_data[["time"]], event = ks_data[["event"]])

  p_values <- map(ks_models, broom::tidy, exponentiate = TRUE)
  hr_values <- map_dbl(p_values, "estimate")
  p_values <- map_dbl(p_values, "p.value")


  cutoff_rank <- which(rank(p_values, ties.method = "first") == 1)

  structure(
    list(p_value = p_values[[cutoff_rank]],
         p_value_adj = p.adjust(p_values, method = p.adjust.method)[[cutoff_rank]],
         p_value_adj_method = p.adjust.method,
         data = select(data, !!!unname(ks_vars)),
         model = ks_models[[cutoff_rank]],
         expression_level = expression_level[[cutoff_rank]][order(ks_data[["row"]])],
         p_values = p_values,
         hr_values = hr_values,
         cutoff_rank = unname(cutoff_rank + 7)), class = "kaplan_scan")
}

#' Extract the threshold expression value able to achieve the biggest survival curve separation upon dichotomization.
#'
#' Uses the KaplanScan algorithm to determine the threshold.
#'
#' @param time,event,expression Numerical vectors containing the follow up time, status indicator and gene expression.
#'
#' @export
best_threshold <- function(time, event, expression) {
  x <- data.frame(time = time,
                  event = event,
                  expression = expression)
  ks <- kaplan_scan(x, time, event, expression)
  sort(expression)[ks[["cutoff_rank"]]]
}
