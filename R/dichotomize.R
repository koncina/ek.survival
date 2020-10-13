NULL
#' @export
dichotomize_survival <- function(data, time, event, expression, method, ...) {
  method_name <- quo_name(enquo(method))
  if (method_name == "kaplan_scan") {
    x <- kaplan_scan({{data}}, {{time}}, {{event}}, {{expression}}, ...)
    expression_level <- x[["expression_level"]]
    cox_model <- x[["model"]]
  } else {
    expr <- pull(data, {{expression}})
    expression_level <- if_else( expr < method(expr), "low", "high")
    cox_model <- ks_coxph(pull(data, {{time}}),pull(data, {{event}}), expression_level)
  }

  tibble(method = method_name,
         expression_level = list(expression_level),
         model = list(cox_model))
}
