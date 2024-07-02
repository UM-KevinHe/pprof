#' @import dplyr rlang
calculate_scores <- function(input_data, obs, exp, ratio, test, alpha) {

  # Check if input_data is a data frame
  if (!is.data.frame(input_data)) stop("input_data must be a data frame")

  # Check if obs, exp, ratio, flag are present in data
  if (!all(c(as.character(rlang::ensym(obs)), as.character(rlang::ensym(exp)),
    as.character(rlang::ensym(ratio))) %in% names(input_data))) {
    stop(paste(as.character(rlang::ensym(obs)), as.character(rlang::ensym(exp)),
          as.character(rlang::ensym(ratio)), "must be column names in data"))
  }

  if (!(test %in% c("exact", "score"))) stop("Argument 'test' NOT as required!", call.=F)

  # Check if alpha is numeric
  if (!is.numeric(alpha)) stop("alpha must be a numeric value")

  if (test == "exact") {
    data <- input_data %>%
      mutate(
        prob = 1 - ppois({{ obs }}, {{ exp }}) + 0.5 * dpois({{ obs }}, {{ exp }}),
        flag = factor(
          ifelse(
            prob < alpha / 2, 1,
            ifelse(prob <= 1 - alpha/2, 0, -1)
          )
        ),
        p = 2 * pmin(prob, 1 - prob)
      )
  }
  else if (test == "score") {
    data <- input_data %>%
      mutate(
        z = ({{ obs }} - {{ exp }}) / sqrt({{ exp }}),
        p = 2 * pmin(pnorm(z), 1 - pnorm(z)),
        flag = factor(
          ifelse(p < alpha, ifelse({{ ratio }} > 1, 1, -1), 0)
        )
      )
  }

  return(data)
}


#' @import dplyr tidyr
calculate_control_limits <- function(data, exp, ratio, flag, test, alpha, target) {
  if (!is.data.frame(data)) stop("data must be a data frame")

  # Check if obs, exp, ratio, flag are present in data
  if (!all(c(as.character(rlang::ensym(exp)), as.character(rlang::ensym(ratio)),
             as.character(rlang::ensym(flag))) %in% names(data))) {
    stop(paste(as.character(rlang::ensym(exp)), as.character(rlang::ensym(ratio)),
               as.character(rlang::ensym(flag)), "must be column names in data"))
  }

  if (!(test %in% c("exact", "score"))) stop("Argument 'test' NOT as required!", call.=F)

  if (!is.numeric(alpha)) stop("alpha must be a numeric vector")

  # Check if method_type is a character
  if (!is.numeric(target)) stop("target must be a numeric value")

  cl_lower <- function(E, alpha) {
    # lower CL for obs
    o_lower <- qpois(alpha / 2, E)
    o_lower <- ifelse(ppois(o_lower - 1, E) + 0.5 * dpois(o_lower, E) >= alpha / 2, o_lower, o_lower + 1)
    lambda_lower <- (dpois(o_lower, E) + 2 * ppois(o_lower - 1, E) - alpha) / (dpois(o_lower, E) + dpois(o_lower - 1, E))
    lower <- pmax(o_lower - lambda_lower, 0)
    return(lower)
  }

  cl_upper <- function(E, alpha) {
    # upper CL for obs
    o_upper <- qpois(1 - alpha / 2, E)
    o_upper <- ifelse(ppois(o_upper - 1, E) + 0.5 * dpois(o_upper, E) >= 1 - alpha / 2, o_upper, o_upper + 1)
    lambda_upper <- (dpois(o_upper, E) + 2 * ppois(o_upper - 1, E) - 2 + alpha) / (dpois(o_upper - 1, E) + dpois(o_upper, E))
    upper <- o_upper - lambda_upper
    return(upper)
  }

  ctrl_limits <- data %>%
    mutate(se = sqrt({{ exp }}) / {{ exp }}) %>%
    arrange(desc(se)) %>%
    mutate(
      flag = {{ flag }},
      precision = 1 / se^2,
      indicator = {{ ratio }},
      exp = {{ exp }}) %>%
    select(se, precision, indicator, exp, flag) %>%
    cross_join(tibble(alpha))

  if (test == "score") {
    ctrl_limits <- ctrl_limits %>%
      mutate(
        lower = target - qnorm(1 - alpha / 2) * sqrt(1 / precision),
        upper = target + qnorm(1 - alpha / 2) * sqrt(1 / precision)
      )
  } else if (test == "exact") {
    ctrl_limits <- ctrl_limits %>%
      mutate(lower = cl_lower(exp, alpha) / exp,
             upper = cl_upper(exp, alpha) / exp)
  }

  ctrl_limits <- ctrl_limits %>%
    mutate(
      alpha = factor(alpha),
      lower = pmax(lower, 0)
    )

  return(ctrl_limits)
}


#' @import ggplot2 RColorBrewer
create_funnel_plot <- function(processed_data,
                               target,
                               alpha,
                               color_palette,
                               labels,
                               shapes,
                               point_size,
                               point_alpha,
                               line_size,
                               target_line_type,
                               xlab = "Precision",
                               ylab = "Outcome",
                               legend_justification = c(1, 1),
                               legend_position = c(1, 1),
                               point_legend_title = "Flagging",
                               linetype_legend_title = "Significance",
                               legend_title_size = 14,
                               legend_size = 14,
                               legend_box = "horizontal",
                               axis_title_size = 14,
                               axis_text_size = 14,
                               plot_title = "Funnel Plot",
                               plot_title_size = 18
) {

  # Check if processed_data is a data frame
  if (!is.data.frame(processed_data)) {
    stop("processed_data must be a data frame")
  }

  data <- processed_data %>% filter(alpha == alpha[1])

  # Ensure that data$flag is a factor
  data$flag <- factor(data$flag, levels = c(-1, 0, 1))

  # Create labels for the legend
  labs_color <- paste0(labels, " (", table(data$flag), ")")

  # Add dummy rows for missing levels with NA values
  missing_levels <- setdiff(c(-1, 0, 1), unique(data$flag))
  if (length(missing_levels) > 0) {
    dummy_data <- data.frame(flag = factor(missing_levels, levels = c(-1, 0, 1)),
                             precision = NA,
                             indicator = NA)
    data <- bind_rows(data, dummy_data)
  }

  num_levels <- length(levels(data$flag))

  # # Check if the length of shapes and color_palette is the same as the number of levels
  # if (length(color_palette) != num_levels) {
  #   stop("The length of color_palette must be the same as the number of levels in data$flag")
  # }
  #
  # # Check if the length of labels is the same as the number of levels
  # if (length(labels) != num_levels) {
  #   stop("The length of labels must be the same as the number of levels in data$flag")
  # }


  # Assign each level of data$flag to a color from the palette
  color_mapping <- setNames(color_palette, levels(data$flag))

  # Assign each level of data$flag to a shape
  shapes_mapping <- setNames(shapes, levels(data$flag))

  # Create a named vector of lables for the legend
  labels <- setNames(labels, levels(data$flag))


  xmax <- max(1 / data$se^2)
  ymax <- max(processed_data$upper)

  labs_linetype <- paste0((1 - alpha) * 100, "%")

  values_linetype <- c('solid', 'dashed', 'dotted', 'dotdash', 'longdash', 'twodash')[1:length(alpha)]

  values_linetype <- values_linetype[order(alpha)]
  labs_linetype <- labs_linetype[order(alpha)]

  plot <-
    ggplot() +
    scale_x_continuous(limits = c(0, xmax),
                       expand = c(1, 1)/50) +
    scale_y_continuous(breaks = round(seq(0, ymax, by=1), 1),
                       limits = c(0, ymax),
                       expand = c(1, 1)/50) +
    geom_point(data = data, aes(x = precision, y = indicator, shape = flag, color = flag), size = point_size, alpha = point_alpha) +
    scale_shape_manual(
      name = bquote(.(point_legend_title) ~ "(" * alpha == .(alpha[1]) * ")"),
      labels = labs_color,
      values = shapes_mapping
    ) +
    scale_color_manual(
      name = bquote(.(point_legend_title) ~ "(" * alpha == .(alpha[1]) * ")"),
      labels = labs_color,
      values = color_mapping
    ) +
    geom_line(data = processed_data, aes(x = precision, y = lower, group = alpha, linetype = alpha), linewidth = line_size) +
    geom_line(data = processed_data, aes(x = precision, y = upper, group = alpha, linetype = alpha), linewidth = line_size) +
    scale_linetype_manual(
      name =  linetype_legend_title,
      values = values_linetype,
      labels = labs_linetype
    ) +
    guides(shape = guide_legend(order = 1), color = guide_legend(order = 1), linetype = guide_legend(reverse = TRUE, order = 2)) +
    geom_hline(yintercept = target, linewidth = line_size, linetype = target_line_type) +
    theme_classic() +
    theme(
      legend.justification = legend_justification,
      legend.position = legend_position,
      legend.box = legend_box,
      legend.title = element_text(size = legend_title_size),
      legend.text = element_text(size = legend_size),
      axis.title = element_text(size = axis_title_size, margin = margin(t = 0, r = 0, b = 0, l = 0)),
      axis.text = element_text(size = axis_text_size),
      plot.title = element_text(hjust = 0.5, size = plot_title_size),
      text = element_text(size = 13),
      legend.background = element_rect(fill = "transparent", colour = NULL, linewidth = 0, linetype = "solid"),
    ) +
    labs(
      x = xlab,
      y = ylab,
      title = plot_title
    )


  return(plot)
}


#' Main Function for Creating a Funnel Plot
#'
#' This function orchestrates the process of creating a funnel plot.
#'
#' @param data A data frame or tibble containing study information.
#' @param target A numeric value representing the target outcome.
#' @param alphas A vector of significance levels.
#' @param color_palette A vector of colors for the plot.
#' @param labels A vector of labels for the plot.
#' @param shapes A vector of shapes for the plot.
#'
#' @import dplyr rlang
#' @return A ggplot object representing the funnel plot.
#'
#' @export
funnel_plot.fe <- function(fit, test = "exact", target = 1, alpha = c(0.05, 0.01),
                        color_palette = c("#E69F00", "#56B4E9", "#009E73"),
                        labels = c("lower", "expected", "higher"),
                        shapes = c(15, 17, 19),
                        point_size = 2,
                        point_alpha = 0.8,
                        line_size = 0.8,
                        target_line_type = "longdash"
) {
  if (missing(fit)) stop ("Argument 'fit' is required!", call.=F)
  if (!class(fit) %in% c("logis_fe")) stop("Object fit is not of the classes 'logis_fe'!", call.=F)

  SR = SR_output(fit, stdz = "indirect", measure = "ratio")
  data = cbind(SR$indirect.ratio, SR$OE$OE_indirect)
  colnames(data) = c("ratio", "Obs", "Exp")

  processed_data <- data %>%
    calculate_scores(Obs, Exp, ratio, test, alpha = alpha[1]) %>%
    calculate_control_limits(Exp, ratio, flag, test, alpha, target)

  plot <- create_funnel_plot(processed_data, target,
                             alpha,
                             color_palette,
                             labels,
                             shapes,
                             point_size,
                             point_alpha,
                             line_size,
                             target_line_type
                             )

  return(plot)
}
