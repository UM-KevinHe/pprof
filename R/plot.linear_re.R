#' @rdname plot.linear_fe
#'
#' @examples
#' data(ExampleDataLinear)
#' Y <- ExampleDataLinear$Y
#' ID <- ExampleDataLinear$ID
#' Z <- ExampleDataLinear$Z
#'
#' fit_re <- linear_re(Y = Y, Z = Z, ID = ID)
#' plot(fit_fe)
#'
#' @seealso \code{\link{linear_re}}, \code{\link{confint.linear_re}}
#'
#' @importFrom ggplot2 ggplot aes theme element_blank labs ggtitle geom_hline geom_point scale_x_continuous scale_y_continuous scale_linetype_manual scale_fill_manual position_jitter geom_errorbar coord_flip theme_minimal theme_bw
#' @importFrom ggpubr annotate_figure ggarrange text_grob
#'
#' @exportS3Method plot linear_re

plot.linear_re <- function(fit, theme = theme_bw(), point_size = 2, point_color = "#475569",
                           medianline_value = 0, medianline_color = "#64748b", medianline_size = 1, medianline_type = "dashed",
                           errorbar_width = 0, errorbar_size = 0.5, errorbar_alpha = 0.5, errorbar_color = "#94a3b8",
                           use_flag = FALSE, flag_color = c("#E69F00", "#56B4E9", "#009E73"),
                           ...) {
  if (missing(fit)) stop ("Argument 'fit' is required!",call.=F)
  if (!class(fit) %in% c("linear_re")) stop("Object fit is not of the classes 'linear_re'!",call.=F)

  args <- list(...)
  stdz <- if ("stdz" %in% names(args)) args$stdz else "indirect"
  alternative <- if("alternative" %in% names(args)) args$alternative else "two.sided"
  if ("option" %in% names(args)) {
    if (args$option == "alpha") {
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
