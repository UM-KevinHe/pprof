#' Get Caterpillar Plot for Standardized Measures
#'
#' Provide caterpillar plot for standardized measures from a fixed effect logistic model.
#'
#' @param fit a model fitted from \code{logis_fe}.
#' @param theme theme for the plot. The default is \code{theme_bw()}.
#' @param point_size size of the points in the caterpillar plot. The default value is 2.
#' @param point_color color of the points in the plot. The default value is "#475569".
#' @param medianline_value value of the horizontal median line. The default value is NULL.
#' @param medianline_color color of the median line. The default value is "#64748b".
#' @param medianline_size size of the median line. The default value is 1.
#' @param medianline_type line type for the median line. The default value is "dashed".
#' @param errorbar_width the width of the error bars (horizontal ends of the CI bars). The default value is 0.
#' @param errorbar_size the thickness of the error bars. The default value is 0.5.
#' @param errorbar_alpha transparency level for the error bars. A value between 0 and 1, where 0 is completely transparent and 1 is fully opaque. The default value is 0.5.
#' @param errorbar_color color of the error bars. The default value is "#94a3b8".
#' @param use_flag logical; if \code{TRUE}, the error bars are colored to show providers' flags based on their performance. The default is \code{FALSE}.
#' @param flag_color vector of colors used for flagging providers when \code{use_flag = TRUE}. The default value is \code{c("#E69F00", "#56B4E9", "#009E73")}.
#' @param \dots Additional arguments to pass to the internal \code{confint.linear_fe} function for a \code{linear_fe} object or
#' the internal \code{confint.linear_re} function for a \code{linear_re} object. These control the type of standardized measure and
#' confidence intervals to be displayed. The \code{option} argument is fixed to \code{"SM"}, meaning that only standardized measures are supported.
#'
#' @details
#' This function creates caterpillar plots to visualize the standardized ratio/rate (indirect or direct) from a fitted fixed effect logistic model.
#' Each provider's standardized measure value is represented as a point, and a reference line is shown at the value specified by `medianline_value` (default is NULL).
#' If `medianline_value` is not specified, for standardized ratios, the reference line will be set at 1;
#' for standardized rates, the reference line will be set at the population rate, which represents the average outcome across all data points.
#' The confidence intervals (CI) are displayed as error bars: for \code{alternative = "two.sided"}, the full CI is shown;
#' for \code{alternative = "greater"}, the error bar extends from the lower bound to the standardized measure value;
#' and for \code{alternative = "less"}, it extends from the standardized measure value to the upper bound.
#'
#' When \code{use_flag = TRUE}, the plot will use colors specified by `flag_color` to show the flags of providers.
#' Each error bar will be colored to reflect the flag, making it easy to identify providers with different performance levels.
#' When \code{use_flag = FALSE}, all error bars will have the same color, specified by `errorbar_color`.
#' This provides a simpler visualization without flagging individual providers.
#'
#' @return A list of ggplot objects containing caterpillar plots for indirect and/or direct standardized ratio/rate based on the value of `stdz` and `measure`.
#'
#' @examples
#' data(data_FE)
#' fit_fe <- logis_fe(Y = data_FE$Y, Z = data_FE$Z, ID = data_FE$ID, message = FALSE)
#' plot(fit_fe)
#'
#' @seealso \code{\link{logis_fe}}, \code{\link{confint.logis_fe}}
#'
#' @importFrom ggplot2 ggplot aes theme element_text element_blank element_line margin labs ggtitle geom_hline geom_point scale_x_continuous scale_y_continuous scale_linetype_manual scale_fill_manual position_jitter geom_errorbar coord_flip theme_minimal
#' @importFrom ggpubr annotate_figure ggarrange text_grob
#'
#' @exportS3Method plot logis_fe


plot.logis_fe <- function(fit, theme = theme_bw(), point_size = 2, point_color = "#475569",
                          medianline_value = NULL, medianline_color = "#64748b", medianline_size = 1, medianline_type = "dashed",
                          errorbar_width = 0, errorbar_size = 0.5, errorbar_alpha = 0.5, errorbar_color = "#94a3b8",
                          use_flag = FALSE, flag_color = c("#E69F00", "#56B4E9", "#009E73"),
                          ...) {
  if (missing(fit)) stop ("Argument 'fit' is required!",call.=F)
  if (!class(fit) %in% c("logis_fe")) stop("Object fit is not of the classes 'logis_fe'!",call.=F)

  args <- list(...)
  stdz <- if ("stdz" %in% names(args)) args$stdz else "indirect"
  measure <- if ("measure" %in% names(args)) args$measure else c("rate", "ratio")
  alternative <- if("alternative" %in% names(args)) args$alternative else "two.sided"
  if ("option" %in% names(args)) {
    if (args$option == "gamma") {
      stop("Caterpillar plot only supports standardized measure and the 'option' argument must be 'SM'.", call. = FALSE)
    }
  }
  CI <- do.call(confint, c(list(fit), args))

  return_ls <- list()



  if ("indirect" %in% stdz) {
    if ("ratio" %in% measure){
      CI$CI.indirect_ratio$prov <- rownames(CI$CI.indirect_ratio)
      medianline_ratio <- if (is.null(medianline_value)) 1 else medianline_value
      #CI$CI.indirect_ratio <- CI$CI.indirect_ratio[!is.infinite(CI$CI.indirect_ratio$CI_ratio.lower) & !is.infinite(CI$CI.indirect_ratio$CI_ratio.upper), ]
      if (alternative == "two.sided") {
        CI$CI.indirect_ratio$flag <- ifelse(CI$CI.indirect_ratio$CI_ratio.upper < medianline_ratio, "Lower",
                                      ifelse(CI$CI.indirect_ratio$CI_ratio.lower > medianline_ratio, "Higher", "Normal"))
      } else if (alternative == "greater") {
        CI$CI.indirect_ratio$flag <- ifelse(CI$CI.indirect_ratio$CI_ratio.lower > medianline_ratio, "Higher", "Normal")
      } else if (alternative == "less") {
        CI$CI.indirect_ratio$flag <- ifelse(CI$CI.indirect_ratio$CI_ratio.upper < medianline_ratio, "Lower", "Normal")
      }

      p_indirect_ratio <- ggplot(CI$CI.indirect_ratio, aes(x = reorder(prov, indirect_ratio), y = indirect_ratio))
      if (use_flag == TRUE) {
        p_indirect_ratio <- p_indirect_ratio +
          geom_errorbar(aes(ymin = if (alternative == "less") indirect_ratio else CI_ratio.lower,
                            ymax = if (alternative == "greater") indirect_ratio else CI_ratio.upper,
                            color = flag),
                        width = errorbar_width, linewidth = errorbar_size, alpha = errorbar_alpha) +
          scale_color_manual(values = flag_color, guide = guide_legend(title = NULL, box.linetype = "solid",
                                                                       override.aes = list(linewidth = 1.5) ))
      } else {
        p_indirect_ratio <- p_indirect_ratio +
          geom_errorbar(aes(ymin = if (alternative == "less") indirect_ratio else CI_ratio.lower,
                            ymax = if (alternative == "greater") indirect_ratio else CI_ratio.upper),
                        width = errorbar_width, linewidth = errorbar_size, alpha = errorbar_alpha, color = errorbar_color)
      }
      p_indirect_ratio <- p_indirect_ratio +
        geom_point(size = point_size, color = point_color) +
        geom_hline(aes(yintercept = medianline_ratio),
                   color = medianline_color, linetype = medianline_type, linewidth = medianline_size) +
        labs(x = "Provider", y = "Indirect Standardized Ratio", title = "Indirect Standardized Ratio Caterpillar Plot") +
        theme +
        theme(
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title = element_text(size = 18),
          axis.text = element_text(size = 15),
          plot.title = element_text(size = 20),
          legend.position = c(0.95, 0.05),
          legend.justification = c("right", "bottom"),
          legend.box.background = element_rect(color = "black", linewidth = 0.5),
          legend.box.margin = margin(5, 5, 5, 5),
          legend.text = element_text(size = 15)
        )
      return_ls$indirect_ratio <- p_indirect_ratio
    }

    if ("rate" %in% measure) {
      CI$CI.indirect_rate$prov <- rownames(CI$CI.indirect_rate)
      medianline_rate <- if (is.null(medianline_value)) attr(CI$CI.indirect_rate, "population_rate") else medianline_value

      if (alternative == "two.sided") {
        CI$CI.indirect_rate$flag <- ifelse(CI$CI.indirect_rate$CI_rate.upper < medianline_rate, "Lower",
                                           ifelse(CI$CI.indirect_rate$CI_rate.lower > medianline_rate, "Higher", "Normal"))
      } else if (alternative == "greater") {
        CI$CI.indirect_rate$flag <- ifelse(CI$CI.indirect_rate$CI_rate.lower > medianline_rate, "Higher", "Normal")
      } else if (alternative == "less") {
        CI$CI.indirect_rate$flag <- ifelse(CI$CI.indirect_rate$CI_rate.upper < medianline_rate, "Lower", "Normal")
      }

      p_indirect_rate <- ggplot(CI$CI.indirect_rate, aes(x = reorder(prov, indirect_rate), y = indirect_rate))
      if (use_flag == TRUE) {
        p_indirect_rate <- p_indirect_rate +
          geom_errorbar(aes(ymin = if (alternative == "less") indirect_rate else CI_rate.lower,
                            ymax = if (alternative == "greater") indirect_rate else CI_rate.upper,
                            color = flag),
                        width = errorbar_width, linewidth = errorbar_size, alpha = errorbar_alpha) +
          scale_color_manual(values = flag_color, guide = guide_legend(title = NULL, box.linetype = "solid",
                                                                       override.aes = list(linewidth = 1.5) ))
      } else {
        p_indirect_rate <- p_indirect_rate +
          geom_errorbar(aes(ymin = if (alternative == "less") indirect_rate else CI_rate.lower,
                            ymax = if (alternative == "greater") indirect_rate else CI_rate.upper),
                        width = errorbar_width, linewidth = errorbar_size, alpha = errorbar_alpha, color = errorbar_color)
      }
      p_indirect_rate <- p_indirect_rate +
        geom_point(size = point_size, color = point_color) +
        geom_hline(aes(yintercept = medianline_rate),
                   color = medianline_color, linetype = medianline_type, linewidth = medianline_size) +
        labs(x = "Provider", y = "Indirect Standardized Rate", title = "Indirect Standardized Rate Caterpillar Plot") +
        theme +
        theme(
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title = element_text(size = 18),
          axis.text = element_text(size = 15),
          plot.title = element_text(size = 20),
          legend.position = c(0.95, 0.05),
          legend.justification = c("right", "bottom"),
          legend.box.background = element_rect(color = "black", linewidth = 0.5),
          legend.box.margin = margin(5, 5, 5, 5),
          legend.text = element_text(size = 15)
        )
      return_ls$indirect_rate <- p_indirect_rate
    }
  }

  if ("direct" %in% stdz) {
    if ("ratio" %in% measure){
      CI$CI.direct_ratio$prov <- rownames(CI$CI.direct_ratio)
      medianline_ratio <- if (is.null(medianline_value)) 1 else medianline_value

      if (alternative == "two.sided") {
        CI$CI.direct_ratio$flag <- ifelse(CI$CI.direct_ratio$CI_ratio.upper < medianline_ratio, "Lower",
                                          ifelse(CI$CI.direct_ratio$CI_ratio.lower > medianline_ratio, "Higher", "Normal"))
      } else if (alternative == "greater") {
        CI$CI.direct_ratio$flag <- ifelse(CI$CI.direct_ratio$CI_ratio.lower > medianline_ratio, "Higher", "Normal")
      } else if (alternative == "less") {
        CI$CI.direct_ratio$flag <- ifelse(CI$CI.direct_ratio$CI_ratio.upper < medianline_ratio, "Lower", "Normal")
      }

      p_direct_ratio <- ggplot(CI$CI.direct_ratio, aes(x = reorder(prov, direct_ratio), y = direct_ratio))
      if (use_flag == TRUE) {
        p_direct_ratio <- p_direct_ratio +
          geom_errorbar(aes(ymin = if (alternative == "less") direct_ratio else CI_ratio.lower,
                            ymax = if (alternative == "greater") direct_ratio else CI_ratio.upper,
                            color = flag),
                        width = errorbar_width, linewidth = errorbar_size, alpha = errorbar_alpha) +
          scale_color_manual(values = flag_color, guide = guide_legend(title = NULL, box.linetype = "solid",
                                                                       override.aes = list(linewidth = 1.5) ))
      } else {
        p_direct_ratio <- p_direct_ratio +
          geom_errorbar(aes(ymin = if (alternative == "less") direct_ratio else CI_ratio.lower,
                            ymax = if (alternative == "greater") direct_ratio else CI_ratio.upper),
                        width = errorbar_width, linewidth = errorbar_size, alpha = errorbar_alpha, color = errorbar_color)
      }
      p_direct_ratio <- p_direct_ratio +
        geom_point(size = point_size, color = point_color) +
        geom_hline(aes(yintercept = medianline_ratio),
                   color = medianline_color, linetype = medianline_type, linewidth = medianline_size) +
        labs(x = "Provider", y = "Direct Standardized Ratio", title = "Direct Standardized Ratio Caterpillar Plot") +
        theme +
        theme(
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title = element_text(size = 18),
          axis.text = element_text(size = 15),
          plot.title = element_text(size = 20),
          legend.position = c(0.95, 0.05),
          legend.justification = c("right", "bottom"),
          legend.box.background = element_rect(color = "black", linewidth = 0.5),
          legend.box.margin = margin(5, 5, 5, 5),
          legend.text = element_text(size = 15)
        )
      return_ls$direct_ratio <- p_direct_ratio
    }

    if ("rate" %in% measure) {
      CI$CI.direct_rate$prov <- rownames(CI$CI.direct_rate)
      medianline_rate <- if (is.null(medianline_value)) attr(CI$CI.direct_rate, "population_rate") else medianline_value

      if (alternative == "two.sided") {
        CI$CI.direct_rate$flag <- ifelse(CI$CI.direct_rate$CI_rate.upper < medianline_rate, "Lower",
                                           ifelse(CI$CI.direct_rate$CI_rate.lower > medianline_rate, "Higher", "Normal"))
      } else if (alternative == "greater") {
        CI$CI.direct_rate$flag <- ifelse(CI$CI.direct_rate$CI_rate.lower > medianline_rate, "Higher", "Normal")
      } else if (alternative == "less") {
        CI$CI.direct_rate$flag <- ifelse(CI$CI.direct_rate$CI_rate.upper < medianline_rate, "Lower", "Normal")
      }

      p_direct_rate <- ggplot(CI$CI.direct_rate, aes(x = reorder(prov, direct_rate), y = direct_rate))
      if (use_flag == TRUE) {
        p_direct_rate <- p_direct_rate +
          geom_errorbar(aes(ymin = if (alternative == "less") direct_rate else CI_rate.lower,
                            ymax = if (alternative == "greater") direct_rate else CI_rate.upper,
                            color = flag),
                        width = errorbar_width, linewidth = errorbar_size, alpha = errorbar_alpha) +
          scale_color_manual(values = flag_color, guide = guide_legend(title = NULL, box.linetype = "solid",
                                                                       override.aes = list(linewidth = 1.5) ))
      } else {
        p_direct_rate <- p_direct_rate +
          geom_errorbar(aes(ymin = if (alternative == "less") direct_rate else CI_rate.lower,
                            ymax = if (alternative == "greater") direct_rate else CI_rate.upper),
                        width = errorbar_width, linewidth = errorbar_size, alpha = errorbar_alpha, color = errorbar_color)
      }
      p_direct_rate <- p_direct_rate +
        geom_point(size = point_size, color = point_color) +
        geom_hline(aes(yintercept = medianline_rate),
                   color = medianline_color, linetype = medianline_type, linewidth = medianline_size) +
        labs(x = "Provider", y = "Direct Standardized Rate", title = "Direct Standardized Rate Caterpillar Plot") +
        theme +
        theme(
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title = element_text(size = 18),
          axis.text = element_text(size = 15),
          plot.title = element_text(size = 20),
          legend.position = c(0.95, 0.05),
          legend.justification = c("right", "bottom"),
          legend.box.background = element_rect(color = "black", linewidth = 0.5),
          legend.box.margin = margin(5, 5, 5, 5),
          legend.text = element_text(size = 15)
        )
      return_ls$direct_rate <- p_direct_rate
    }
  }

  return(return_ls)
}
