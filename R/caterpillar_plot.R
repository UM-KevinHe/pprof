#' Get Caterpillar Plot for Standardized Measures
#'
#' Generates a caterpillar plot for standardized measures from different models using a provided CI dataframe.
#'
#' @param CI a dataframe from `confint` containing the standardized measure values, along with their confidence interval lower and upper bounds.
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
#'
#' @details
#' This function creates caterpillar plots to visualize the standardized measures (indirect or direct).
#' The input `CI` must be a dataframe output from package `pprof`'s `confint` function.
#' Each provider's standardized measure value is represented as a point, and a reference line is shown at the value specified by `medianline_value` (default is NULL).
#' If `medianline_value` is not specified, for linear FE or RE models with indirect or direct standardized measures, it will be set to 0;
#' for logistic FE models with indirect or direct ratios, it will be set to 1;
#' and for logistic FE with indirect or direct rates, it will be set to the population rate, which represents the average outcome across all data points.
#' The confidence intervals (CI) are displayed as error bars: for \code{alternative = "two.sided"}, the full CI is shown;
#' for \code{alternative = "greater"}, the error bar extends from the lower bound to the standardized measure value;
#' and for \code{alternative = "less"}, it extends from the standardized measure value to the upper bound.
#'
#' When \code{use_flag = TRUE}, the plot will use colors specified by `flag_color` to show the flags of providers.
#' Each error bar will be colored to reflect the flag, making it easy to identify providers with different performance levels.
#' When \code{use_flag = FALSE}, all error bars will have the same color, specified by `errorbar_color`.
#' This provides a simpler visualization without flagging individual providers.
#'
#' @return A ggplot object which is a caterpillar plot for the standardized measures.
#'
#' @examples
#' data(ExampleDataLinear)
#' fit_linear <- linear_fe(Y = ExampleDataLinear$Y, Z = ExampleDataLinear$Z, ID = ExampleDataLinear$ID)
#' CI_linear <- confint(fit_linear)
#' caterpillar_plot(CI_linear$CI.indirect, use_flag =T)
#'
#' data(data_FE)
#' fit_logis <- logis_fe(Y = data_FE$Y, Z = data_FE$Z, ID = data_FE$ID, message = FALSE)
#' CI_logis <- confint(fit_linear)
#' caterpillar_plot(CI$CI.indirect_rate)
#'
#' @seealso \code{\link{confint.linear_fe}}, \code{\link{confint.linear_re}}, \code{\link{confint.logis_fe}}
#'
#' @importFrom ggplot2 ggplot aes theme element_blank labs ggtitle geom_hline geom_point scale_x_continuous scale_y_continuous scale_linetype_manual scale_fill_manual position_jitter geom_errorbar coord_flip theme_minimal theme_bw
#' @importFrom ggpubr annotate_figure ggarrange text_grob
#' @importFrom dplyr arrange
#'
#' @export

caterpillar_plot <- function(CI, theme = theme_bw(), point_size = 2, point_color = "#475569",
                             medianline_value = NULL, medianline_color = "#64748b", medianline_size = 1, medianline_type = "dashed",
                             errorbar_width = 0, errorbar_size = 0.5, errorbar_alpha = 0.5, errorbar_color = "#94a3b8",
                             use_flag = FALSE, flag_color = c("#E69F00", "#56B4E9", "#009E73")) {
  if (missing(CI)) stop ("Argument 'CI' is required!",call.=F)
  if (!class(CI) %in% c("data.frame")) stop("Object CI should be a data frame!",call.=F)

  colnames(CI) <- c("SM", "Lower", "Upper")
  CI$prov <- rownames(CI)

  if (attr(CI, "model") == "FE linear" | attr(CI, "model") == "RE linear") {
    medianline_value <- if (is.null(medianline_value)) 0 else medianline_value
  }
  else if (attr(CI, "model") == "FE logis") {
    if (grepl("Ratio", attr(CI, "description"))) {
      medianline_value <- if (is.null(medianline_value)) 1 else medianline_value
    }
    else if (grepl("Rate", attr(CI, "description"))) {
      medianline_value <- if (is.null(medianline_value)) attr(CI, "population_rate") else medianline_value
    }
  }

  if (attr(CI, "type") == "two-sided") {
    CI$flag <- ifelse(CI$Upper < medianline_value, "Lower",
                      ifelse(CI$Lower > medianline_value, "Higher", "Normal"))
    # CI$flag <- factor(CI$flag, levels = c("Normal", "Lower", "Higher"), ordered = T)
  } else if (attr(CI, "type") == "upper one-sided") {
    CI$flag <- ifelse(CI$Lower > medianline_value, "Higher", "Normal")
  } else if (attr(CI, "type") == "lower one-sided") {
    CI$flag <- ifelse(CI$Upper < medianline_value, "Lower", "Normal")
  }

  # CI$flag <- factor(CI$flag, levels = c("Normal", "Lower", "Higher"), ordered = T)
  caterpillar_p <- ggplot(CI, aes(x = reorder(prov, SM), y = SM))
  if (use_flag == TRUE) {
    caterpillar_p <- caterpillar_p +
      geom_errorbar(aes(ymin = if (attr(CI, "type") == "lower one-sided") SM else Lower,
                        ymax = if (attr(CI, "type") == "upper one-sided") SM else Upper,
                        color = flag),
                    width = errorbar_width, linewidth = errorbar_size, alpha = errorbar_alpha) +
      scale_color_manual(values = flag_color, guide = guide_legend(title = NULL, box.linetype = "solid",
                                                                   override.aes = list(linewidth = 1.5) ))

  } else {
    caterpillar_p <- caterpillar_p +
      geom_errorbar(aes(ymin = if (attr(CI, "type") == "lower one-sided") SM else Lower,
                        ymax = if (attr(CI, "type") == "upper one-sided") SM else Upper),
                    width = errorbar_width, linewidth = errorbar_size, alpha = errorbar_alpha, color = errorbar_color)
  }

  caterpillar_p <- caterpillar_p +
    geom_point(size = point_size, color = point_color) +
    geom_hline(aes(yintercept = medianline_value),
               color = medianline_color, linetype = medianline_type, linewidth = medianline_size) +
    scale_x_discrete(expand = expansion(add = 5)) +
    labs(x = "Provider", y = attr(CI, "description"), title = paste(attr(CI, "description"), "Caterpillar Plot")) +
    theme +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title = element_text(size = 18, face = "bold"),
      axis.text = element_text(size = 15, face = "bold"),
      plot.title = element_text(size = 20, face = "bold"),
      legend.position = c(0.95, 0.05),
      legend.justification = c("right", "bottom"),
      legend.box.background = element_rect(color = "black", linewidth = 0.5),
      legend.box.margin = margin(5, 5, 5, 5),
      legend.text = element_text(size = 15, face = "bold")
    )

  return(caterpillar_p)
}
