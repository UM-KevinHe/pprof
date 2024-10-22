#' Get Bar Plot for flags of each provider
#'
#' Generates a bar plot for flagging percentage based on grouped sample sizes.
#'
#' @param flag_df a data frame from `test` function containing the flag of each provider.
#' @param group_num number of groups to divide providers into based on their sample size. The default is 4.
#' @param bar_colors a vector of colors used to fill the bars representing the categories. The default is c("#66c2a5", "#fc8d62", "#8da0cb").
#' @param bar_width Width of the bars in the bar chart. The default is 0.7.
#' @param theme theme for the plot. The default is theme_minimal().
#' @param label_color color of the text labels inside the bars. The default is "black".
#' @param label_size size of the text labels inside the bars. The default is 4.
#'
#' @details
#' This function generates a bar chart to visualize the percentage of flagging results based on provider size.
#' The input data frame `test_df` must be the output from package `pprof`'s `test` function.
#' The providers are grouped into a specified number of groups (`group_num`) based on their sample sizes,
#' and an additional "overall" group is included to show the flagging results across all providers.
#'
#' @return A ggplot object representing the bar chart of flagging results.
#'
#' @examples
#' data(ExampleDataLinear)
#' fit_linear <- linear_fe(Y = ExampleDataLinear$Y, Z = ExampleDataLinear$Z, ID = ExampleDataLinear$ID)
#' test_linear <- test(fit_linear)
#' bar_plot(CI_linear$CI.indirect, use_flag =T)
#'
#' data(data_FE)
#' fit_logis <- logis_fe(Y = data_FE$Y, Z = data_FE$Z, ID = data_FE$ID, message = FALSE)
#' test_logis <- test(fit_linear)
#' bar_plot(CI$CI.indirect_rate)
#'
#' @seealso \code{\link{test.linear_fe}}, \code{\link{test.linear_re}}, \code{\link{test.logis_fe}}
#'
#' @importFrom ggplot2 ggplot geom_bar geom_text labs aes theme element_blank ggtitle scale_x_continuous scale_y_continuous scale_linetype_manual scale_fill_manual theme_minimal theme_bw
#' @importFrom tidyverse group_by summerise mutate
#'
#' @export

bar_plot <- function(flag_df, group_num = 4,
                     bar_colors = c("#66c2a5", "#fc8d62", "#8da0cb"), bar_width = 0.7,
                     theme = theme_minimal(), label_color = "black", label_size = 4) {
  if (missing(flag_df)) stop ("Argument 'flag_df' is required!",call.=F)
  if (!class(flag_df) %in% c("data.frame")) stop("Object flag_df should be a data frame!",call.=F)
  if (!"flag" %in% colnames(flag_df) || is.null(attr(flag_df, "provider size"))) {
    stop("Dataframe must contain a 'flag' column and an attribute 'provider size'.")
  }

  flag_df$category <- factor(flag_df$flag, levels = c(1, 0, -1), labels = c("higher", "as expected", "lower"))

  provider_size <- attr(flag_df, "provider size")
  flag_df$size <- cut(provider_size,
                      breaks = quantile(provider_size, probs = (0:group_num)/group_num, na.rm = TRUE),
                      include.lowest = TRUE,
                      labels = paste0("Q", 1:group_num))
  flag_df$size <- factor(flag_df$size, levels = c(paste0("Q", 1:group_num), "Overall"))

  flag_df_overall <- flag_df
  flag_df_overall$size <- "Overall"
  flag_df <- rbind(flag_df, flag_df_overall)

  # Calculate percentage of each flag in each group
  df_long <- flag_df %>%
    group_by(size, category) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(size) %>%
    mutate(value = count / sum(count))

  # Plot the bar chart
  p <- ggplot(df_long, aes(x = size, y = value, fill = category)) +
    geom_bar(stat = "identity", color = "black", width = 0.7) +
    geom_text(aes(label = scales::percent(value, accuracy = 0.1)),
              position = position_stack(vjust = 0.5),
              color = label_color, size = label_size) +
    labs(x = "Provider Size",
         y = "Flagging Percentage",
         title = "Flagging Results Based on Provider Size",
         fill = "Category") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
    theme +
    scale_fill_manual(values = bar_colors) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      axis.title.x = element_text(face = "bold", size = 14),
      axis.title.y = element_text(face = "bold", size = 14),
      axis.text = element_text(size = 12),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 12),
      panel.grid.major = element_line(size = 0.5, linetype = 'solid', color = 'grey80'),
      panel.grid.minor = element_blank()
    )

  return(p)
}

# flag_df %>% filter(size == "Q1") %>% mutate(cat_low = category == "as expected") %>% pull(cat_low) %>% mean()
