#' Get Caterpillar Plot for Standardized Measures
#'
#' Provide caterpillar plot for standardized measures from a fixed/random effect linear model.
#'
#' @param fit a model fitted from \code{linear_fe} or \code{linear_re}.
#' @param theme theme for the plot. The default is \code{theme_bw()}.
#' @param point_size size of the points in the caterpillar plot. The default value is 2.
#' @param point_color color of the points in the plot. The default value is "#475569".
#' @param medianline_value value of the horizontal median line. The default value is 0.
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
#' This function creates caterpillar plots to visualize the standardized measures (indirect or direct) from a fitted fixed effect linear model.
#' Each provider's standardized measure value is represented as a point, and a reference line is shown at the value specified by `medianline_value` (default is 0).
#' The confidence intervals (CI) are displayed as error bars: for \code{alternative = "two.sided"}, the full CI is shown;
#' for \code{alternative = "greater"}, the error bar extends from the lower bound to the standardized measure value;
#' and for \code{alternative = "less"}, it extends from the standardized measure value to the upper bound.
#'
#' When \code{use_flag = TRUE}, the plot will use colors specified by `flag_color` to show the flags of providers.
#' Each error bar will be colored to reflect the flag, making it easy to identify providers with different performance levels.
#' When \code{use_flag = FALSE}, all error bars will have the same color, specified by `errorbar_color`.
#' This provides a simpler visualization without flagging individual providers.
#'
#' @return A list of ggplot objects containing caterpillar plots for indirect and/or direct standardized measures.
#' \item{indirect}{a ggplot object for the indirect standardized measure, if `stdz` includes \code{"indirect"}.}
#' \item{direct}{a ggplot object for the direct standardized measure, if `stdz` includes \code{"direct"}.}
#'
#' @examples
#' data(ExampleDataLinear)
#' Y <- ExampleDataLinear$Y
#' ID <- ExampleDataLinear$ID
#' Z <- ExampleDataLinear$Z
#'
#' fit_fe <- linear_fe(Y = Y, Z = Z, ID = ID)
#' plot(fit_fe)
#'
#' @seealso \code{\link{linear_fe}}, \code{\link{confint.linear_fe}}
#'
#' @importFrom ggplot2 ggplot aes theme element_blank labs ggtitle geom_hline geom_point scale_x_continuous scale_y_continuous scale_linetype_manual scale_fill_manual position_jitter geom_errorbar coord_flip theme_minimal theme_bw
#' @importFrom ggpubr annotate_figure ggarrange text_grob
#'
#' @exportS3Method plot linear_fe

plot.linear_fe <- function(fit, theme = theme_bw(), point_size = 2, point_color = "#475569",
                           medianline_value = 0, medianline_color = "#64748b", medianline_size = 1, medianline_type = "dashed",
                           errorbar_width = 0, errorbar_size = 0.5, errorbar_alpha = 0.5, errorbar_color = "#94a3b8",
                           use_flag = FALSE, flag_color = c("#E69F00", "#56B4E9", "#009E73"),
                           ...) {
  if (missing(fit)) stop ("Argument 'fit' is required!",call.=F)
  if (!class(fit) %in% c("linear_fe")) stop("Object fit is not of the classes 'linear_fe'!",call.=F)

  args <- list(...)
  stdz <- if ("stdz" %in% names(args)) args$stdz else "indirect"
  alternative <- if("alternative" %in% names(args)) args$alternative else "two.sided"
  if ("option" %in% names(args)) {
    if (args$option == "gamma") {
      stop("Caterpillar plot only supports standardized measure and the 'option' argument must be 'SM'.", call. = FALSE)
    }
  }
  CI <- do.call(confint, c(list(fit), args))

  return_ls <- list()

  if ("indirect" %in% stdz) {
    CI$CI.indirect$prov <- rownames(CI$CI.indirect)
    if (alternative == "two.sided") {
      CI$CI.indirect$flag <- ifelse(CI$CI.indirect$indirect.Upper < medianline_value, "Lower",
                                    ifelse(CI$CI.indirect$indirect.Lower > medianline_value, "Higher", "Normal"))
    } else if (alternative == "greater") {
      CI$CI.indirect$flag <- ifelse(CI$CI.indirect$indirect.Lower > medianline_value, "Higher", "Normal")
    } else if (alternative == "less") {
      CI$CI.indirect$flag <- ifelse(CI$CI.indirect$indirect.Upper < medianline_value, "Lower", "Normal")
    }

    p_indirect <- ggplot(CI$CI.indirect, aes(x = reorder(prov, Indirect.Difference), y = Indirect.Difference))
    if (use_flag == TRUE) {
      p_indirect <- p_indirect +
        geom_errorbar(aes(ymin = if (alternative == "less") Indirect.Difference else indirect.Lower,
                          ymax = if (alternative == "greater") Indirect.Difference else indirect.Upper,
                          color = flag),
                      width = errorbar_width, linewidth = errorbar_size, alpha = errorbar_alpha) +
        scale_color_manual(values = flag_color, guide = guide_legend(title = NULL, box.linetype = "solid",
                                                                     override.aes = list(linewidth = 1.5) ))
    } else {
      p_indirect <- p_indirect +
        geom_errorbar(aes(ymin = if (alternative == "less") Indirect.Difference else indirect.Lower,
                          ymax = if (alternative == "greater") Indirect.Difference else indirect.Upper),
                      width = errorbar_width, linewidth = errorbar_size, alpha = errorbar_alpha, color = errorbar_color)
    }

    p_indirect <- p_indirect +
      geom_point(size = point_size, color = point_color) +
      geom_hline(aes(yintercept = medianline_value),
                 color = medianline_color, linetype = medianline_type, linewidth = medianline_size) +
      scale_x_discrete(expand = expansion(add = 5)) +
      labs(x = "Provider", y = "Indirect Standardized Difference", title = "Indirect Standardized Difference Caterpillar Plot") +
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
    return_ls$indirect <- p_indirect
  }

  if ("direct" %in% stdz) {
    CI$CI.direct$prov <- rownames(CI$CI.direct)
    if (alternative == "two.sided") {
      CI$CI.direct$flag <- ifelse(CI$CI.direct$direct.Upper < medianline_value, "Lower",
                                    ifelse(CI$CI.direct$direct.Lower > medianline_value, "Higher", "Normal"))
    } else if (alternative == "greater") {
      CI$CI.direct$flag <- ifelse(CI$CI.direct$direct.Lower > medianline_value, "Higher", "Normal")
    } else if (alternative == "less") {
      CI$CI.direct$flag <- ifelse(CI$CI.direct$direct.Upper < medianline_value, "Lower", "Normal")
    }

    p_direct <- ggplot(CI$CI.direct, aes(x = reorder(prov, Direct.Difference), y = Direct.Difference))
    if (use_flag == TRUE) {
      p_direct <- p_direct +
        geom_errorbar(aes(ymin = if (alternative == "less") Direct.Difference else direct.Lower,
                          ymax = if (alternative == "greater") Direct.Difference else direct.Upper,
                          color = flag),
                      width = errorbar_width, linewidth = errorbar_size, alpha = errorbar_alpha) +
        scale_color_manual(values = flag_color, guide = guide_legend(title = NULL, box.linetype = "solid",
                                                                     override.aes = list(linewidth = 1.5) ))
    } else {
      p_direct <- p_direct +
        geom_errorbar(aes(ymin = if (alternative == "less") Direct.Difference else direct.Lower,
                          ymax = if (alternative == "greater") Direct.Difference else direct.Upper),
                      width = errorbar_width, linewidth = errorbar_size, alpha = errorbar_alpha, color = errorbar_color)
    }

    p_direct <- p_direct +
      geom_point(size = point_size, color = point_color) +
      geom_hline(aes(yintercept = medianline_value),
                 color = medianline_color, linetype = medianline_type, linewidth = medianline_size) +
      labs(x = "Provider", y = "Direct Standardized Difference", title = "Direct Standardized Measure Caterpillar Plot") +
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
    return_ls$direct <- p_direct
  }

  return(return_ls)
}
