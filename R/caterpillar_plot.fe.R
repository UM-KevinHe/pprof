#' Caterpillar Plot from a Fitted `logis_fe` object
#'
#' This function creates a caterpillar plot of the provider effects or standardization ratios/rates from a fitted `logis_fe` model.
#'
#' @param fit model obtained from `logis_fe`.
#'
#' @param parm specify a subset of which providers are to be plotted. All providers are included by default.
#'
#' @param level the level of confidence intervals. The default value is `0.95`.
#'
#' @param option the confidence interval for the function's output, whether it is for gamma or standardization ratios/rates.
#'   \itemize{
#'   \item "gamma": provider effect
#'   \item "SR": standardization ratios/rates
#'   }
#'
#' @param test a character string specifying the type of testing method. Defaulting to "exact".
#'   \itemize{
#'   \item "exact": exact test
#'   \item "wald": wald test
#'   \item "score": score test
#'   }
#'
#' @param stdz if option = 'SR', a character string specifying the standardization method. Defaulting to "indirect".
#'   \itemize{
#'   \item "indirect": using indirect standardized method
#'   \item "direct": using direct standardized method
#'   }
#'
#' @param measure if option = 'SR', a boolean indicating whether the output measure is "ratio" or "rate". Both "rate" and "ratio" will be provided by default.
#'   \itemize{
#'   \item "rate": output the standardized rate. The "rate" has been restricted to 0% - 100%.
#'   \item "ratio":  output the standardized ratio
#'   }
#'
#' @param dot_size Size of the points in the plot.
#'
#' @param jitter_width Width of the jitter applied to the points to reduce overlap.
#'
#' @param median_color Color of the median line.
#'
#' @param median_linetype Type of line used for the median line, defaulting to 'dashed'.
#'
#' @param median_size Thickness of the median line.
#'
#' @param errorbar_width Width of the end caps on the error bars.
#'
#' @param errorbar_size Thickness of the error bars.
#'
#' @param ...
#'
#' @return A plot is produced, and nothing is returned.
#'
#' @importFrom ggplot2 ggplot aes theme element_text element_blank element_line margin labs ggtitle geom_hline geom_point scale_x_continuous scale_y_continuous scale_linetype_manual scale_fill_manual position_jitter geom_errorbar coord_flip theme_minimal
#' @importFrom ggpubr annotate_figure ggarrange text_grob
#'
#' @export
#'
#' @examples
#' data(data_FE)
#' fit_fe <- logis_fe(data_FE$Y, data_FE$Z, data_FE$ID, message = FALSE)
#' caterpillar_plot.fe(fit_fe, option = "gamma")

caterpillar_plot.fe <- function(fit, parm, level = 0.95, test = "exact", option = "SR",
                          stdz = "indirect", measure = c("rate", "ratio"),
                          point_size = 2, jitter_width = 0,
                          median_linecolor = "red", median_linetype = "dashed", median_linesize = 1,
                          errorbar_width = 0.2, errorbar_size = 0.5) {

  CI = confint(fit, parm, level, test, option, stdz, measure)

  if (option == "gamma") {
    CI$prov <- rownames(CI)

    # CI$inf_lower <- is.infinite(CI$gamma.lower)
    # CI$inf_upper <- is.infinite(CI$gamma.upper)
    #
    # CI$gamma.lower <- ifelse(is.infinite(CI$gamma.lower), min(CI$gamma)-10, CI$gamma.lower)
    # CI$gamma.upper <- ifelse(is.infinite(CI$gamma.upper), max(CI$gamma)+10, CI$gamma.upper)

    # p <- ggplot(CI, aes(x = reorder(prov, gamma), y = gamma)) +
    #   geom_point(position = position_jitter(width = jitter_width), size = dot_size) +
    #   geom_errorbar(aes(ymin = gamma.lower, ymax = gamma.upper,
    #                     linetype = ifelse(CI$inf_lower | CI$inf_upper, "dotted", "solid")), width = errorbar_width) +
    #   geom_segment(data = CI[CI$inf_lower, ], aes(x = reorder(prov, gamma),
    #                                               xend = reorder(prov, gamma),
    #                                               y = gamma, yend = gamma.lower),
    #                arrow = arrow(length = unit(0.02, "npc")), linetype = "dotted") +
    #   geom_segment(data = CI[CI$inf_upper, ], aes(x = reorder(prov, gamma),
    #                                               xend = reorder(prov, gamma),
    #                                               y = gamma, yend = gamma.upper),
    #                arrow = arrow(length = unit(0.02, "npc")), linetype = "dotted") +
    #   geom_hline(aes(yintercept = median(CI$gamma, na.rm = TRUE)), color = median_color, linetype = median_linetype) +
    #   scale_linetype_manual(values = c("dotted", "solid")) +
    #   coord_flip() +
    #   labs(x = xlab, y = ylab, title = title, subtitle = subtitle) +
    #   theme

    CI <- CI[!is.infinite(CI$gamma.lower) & !is.infinite(CI$gamma.upper), ]

    p <- ggplot(CI, aes(x = reorder(prov, gamma), y = gamma)) +
      geom_point(position = position_jitter(width = jitter_width), size = point_size) +
      geom_errorbar(aes(ymin = gamma.lower, ymax = gamma.upper), width = errorbar_width, linewidth = errorbar_size) +
      geom_hline(aes(yintercept = median(gamma, na.rm = TRUE)), color = median_linecolor, linetype = median_linetype, size = median_linesize) +
      coord_flip() +
      labs(x = "Provider ID", y = "Estimated Value", title = "Provider Effects Confidence Intervals") +
      theme_minimal()

    return(p)
  }

  else if (option == "SR") {
    return_ls <- list()
    if ("indirect" %in% stdz) {
      if ("ratio" %in% measure){
        CI$CI.indirect_ratio$prov <- rownames(CI$CI.indirect_ratio)

        CI$CI.indirect_ratio <- CI$CI.indirect_ratio[!is.infinite(CI$CI.indirect_ratio$CI_ratio.lower) & !is.infinite(CI$CI.indirect_ratio$CI_ratio.upper), ]

        return_ls$indirect_ratio <- ggplot(CI$CI.indirect_ratio, aes(x = reorder(prov, indirect_ratio), y = indirect_ratio)) +
          geom_point(position = position_jitter(width = jitter_width), size = dot_size) +
          geom_errorbar(aes(ymin = CI_ratio.lower, ymax = CI_ratio.upper), width = errorbar_width, linewidth = errorbar_size) +
          geom_hline(aes(yintercept = median(CI$CI.indirect_ratio$indirect_ratio, na.rm = TRUE)),
                     color = median_linecolor, linetype = median_linetype, size = median_linesize) +
          coord_flip() +
          labs(x = "Provider ID", y = "Estimated Value", title = "Indirect Ratio Confidence Intervals") +
          theme_minimal()
      }

      if ("rate" %in% measure) {
        CI$CI.indirect_rate$prov <- rownames(CI$CI.indirect_rate)

        CI$CI.indirect_rate <- CI$CI.indirect_rate[!is.infinite(CI$CI.indirect_rate$CI_rate.lower) & !is.infinite(CI$CI.indirect_rate$CI_rate.upper), ]

        return_ls$indirect_rate <- ggplot(CI$CI.indirect_rate, aes(x = reorder(prov, indirect_rate), y = indirect_rate)) +
          geom_point(position = position_jitter(width = jitter_width), size = dot_size) +
          geom_errorbar(aes(ymin = CI_rate.lower, ymax = CI_rate.upper), width = errorbar_width, linewidth = errorbar_size) +
          geom_hline(aes(yintercept = median(CI$CI.indirect_rate$indirect_rate, na.rm = TRUE)),
                     color = median_linecolor, linetype = median_linetype, size = median_linesize) +
          coord_flip() +
          labs(x = "Provider ID", y = "Estimated Value", title = "Indirect Rate Confidence Intervals") +
          theme_minimal()
      }
    }

    if ("direct" %in% stdz) {
      if ("ratio" %in% measure){
        CI$CI.direct_ratio$prov <- rownames(CI$CI.direct_ratio)

        CI$CI.direct_ratio <- CI$CI.direct_ratio[!is.infinite(CI$CI.direct_ratio$CI_ratio.lower) & !is.infinite(CI$CI.direct_ratio$CI_ratio.upper), ]

        return_ls$direct_ratio <- ggplot(CI$CI.direct_ratio, aes(x = reorder(prov, direct_ratio), y = direct_ratio)) +
          geom_point(position = position_jitter(width = jitter_width), size = dot_size) +
          geom_errorbar(aes(ymin = CI_ratio.lower, ymax = CI_ratio.upper), width = errorbar_width, linewidth = errorbar_size) +
          geom_hline(aes(yintercept = median(CI$CI.direct_ratio$direct_ratio, na.rm = TRUE)),
                     color = median_linecolor, linetype = median_linetype, size = median_linesize) +
          coord_flip() +
          labs(x = "Provider ID", y = "Estimated Value", title = "Direct Ratio Confidence Intervals") +
          theme_minimal()
      }

      if ("rate" %in% measure) {
        CI$CI.direct_rate$prov <- rownames(CI$CI.direct_rate)

        CI$CI.direct_rate <- CI$CI.direct_rate[!is.infinite(CI$CI.direct_rate$CI_rate.lower) & !is.infinite(CI$CI.direct_rate$CI_rate.upper), ]

        return_ls$direct_rate <- ggplot(CI$CI.direct_rate, aes(x = reorder(prov, direct_rate), y = direct_rate)) +
          geom_point(position = position_jitter(width = jitter_width), size = dot_size) +
          geom_errorbar(aes(ymin = CI_rate.lower, ymax = CI_rate.upper), width = errorbar_width, linewidth = errorbar_size) +
          geom_hline(aes(yintercept = median(CI$CI.direct_rate$direct_rate, na.rm = TRUE)),
                     color = median_linecolor, linetype = median_linetype, size = median_linesize) +
          coord_flip() +
          labs(x = "Provider ID", y = "Estimated Value", title = "Direct Rate Confidence Intervals") +
          theme_minimal()
      }
    }

    return(return_ls)
  }

}
