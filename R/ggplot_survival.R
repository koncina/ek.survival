#' Wrapper around `ggplot` and `ggkm` to produce repeatedly used Kaplan-Meier curves
#'
#' The function add statistical annotations (HR and p-value). Grouping the data will show facets with the annotations.
#'
#' @importFrom rlang quo_name enquo
#' @importFrom ggplot2 ggplot expand_limits theme_classic
#' @importFrom ggkm geom_km geom_kmticks
#'
#' @param data A tibble
#' @param time,event,predictor the variable names to produce the Kaplan-Meier plot.
#' @param hr Boolean telling whether to show the hazard ratio and p-value.
#' @param count_npcx,count_npcy position of the counts. If count_npcx is NULL, will try to spread them horizontally.
#'
#' @export
#'
ggplot_survival <- function(data, time, event, predictor = NULL, hr = TRUE, count_npcx = NULL, count_npcy = 0.05) {
  p <- ggplot(data, aes(time = {{time}}, status = {{event}}, colour = {{predictor}})) +
    ggkm::geom_km() +
    ggkm::geom_kmticks(size = 1) +
    expand_limits(y = 0) +
    theme_classic(14)

  # TODO: find a better test...
  if (quo_name(enquo(predictor)) == "NULL") {
    n_levels <- 1
  } else {
    n_levels <- pull(data, {{predictor}}) %>% as_factor() %>% fct_unique() %>% length()
  }

  if (!is_empty(groups(data))) p <- p + facet_wrap(groups(data))

  if (isTRUE(hr) && n_levels == 2) {
    get_stats <- tidy_coxph(data, time = {{time}}, event = {{event}}, predictor = {{predictor}}) %>%
      mutate(across(c(p_val), round_nsmall, 3),
             across(starts_with("hr"), round_nsmall, 2))
    p <- p +
      ggpmisc::geom_text_npc(data = get_stats ,
                             aes(label = paste("p-val =", p_val)),
                             npcx = 0.1, npcy = 0.15, hjust = 0) +
      ggpmisc::geom_text_npc(data = get_stats,
                             aes(label = paste("HR =", hr)), npcx = 0.1, npcy = 0.25, hjust = 0)

  } else if (n_levels > 2 && isTRUE(hr)) {
    warning("too many levels for automated annotations", call. = FALSE)
  }

  if (n_levels > 1) {
    count_names <- pull(data, {{predictor}}) %>%
      unique()
    if (is.null(count_npcx)) {
      count_npcx <- seq(n_levels) / (n_levels + 1)
      count_npcx <- scale(count_npcx) / 10 + count_npcx
    }

    if (length(count_npcy) == 1) count_npcy <- rep(count_npcy, n_levels)

    count_npcx <- set_names(count_npcx, count_names)
    count_npcy <- set_names(count_npcy, count_names)

    p <- p + ggpmisc::geom_text_npc(aes(time = {{time}}, status = {{event}}, colour = {{predictor}},
                                        npcx = count_npcx[{{predictor}}], npcy = count_npcy[{{predictor}}]),
                                    stat = "km_summary", hjust = 0.5, fontface = "italic", size = 3)
  } else {
    if (is.null(count_npcx)) count_npcx <- 0.5
    p <- p + ggpmisc::geom_text_npc(aes(time = {{time}}, status = {{event}},
                                        npcx = count_npcx[1], npcy = count_npcy[1]),
                                    stat = "km_summary", hjust = 0.5, fontface = "italic", size = 3)
  }
  p
}
