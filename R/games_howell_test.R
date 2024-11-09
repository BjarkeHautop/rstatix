#'Games Howell Post-hoc Tests with Adjustable P-value Correction
#'
#'@description Performs Games-Howell test with customizable p-value adjustment.
#'@param data A data frame containing the variables in the formula.
#'@param formula A formula specifying the outcome and group variables.
#'@param conf.level Confidence level for confidence intervals.
#'@param detailed If TRUE, return additional columns such as standard error, method, and statistic.
#'@param p.adjust.method The p-value adjustment method (default is "Tukey").
#'@return A data frame with the results of the Games-Howell test.
#'@export
games_howell_test <- function(data, formula, conf.level = 0.95, detailed = FALSE, p.adjust.method = "Tukey") {
  args <- as.list(environment()) %>%
    .add_item(method = "games_howell_test")

  results <- data %>%
    doo(.games_howell_test, formula, conf.level = conf.level, p.adjust.method = p.adjust.method)

  if (!detailed) {
    results <- results %>%
      select(-.data$se, -.data$method, -.data$statistic, -.data$df, -.data$n1, -.data$n2)
  }

  results %>%
    set_attrs(args = args) %>%
    add_class(c("rstatix_test", "games_howell_test"))
}

.games_howell_test <- function(data, formula, conf.level = 0.95, p.adjust.method = "Tukey") {
  outcome <- get_formula_left_hand_side(formula)
  group <- get_formula_right_hand_side(formula)
  number.of.groups <- guess_number_of_groups(data, group)

  if (number.of.groups == 1) {
    stop("All observations are in the same group")
  }

  data <- data %>%
    select(!!!syms(c(outcome, group))) %>%
    get_complete_cases() %>%
    .as_factor(group)

  x <- data %>% pull(!!outcome)
  g <- data %>% pull(!!group)

  if (!all(is.finite(g))) stop("All group levels must be finite")

  grp.sizes <- tapply(x, g, length)
  nb.groups <- length(grp.sizes)
  grp.means <- tapply(x, g, mean)
  grp.vars <- tapply(x, g, stats::var)

  get_mean_diff <- function(i, j) grp.means[i] - grp.means[j]
  get_weltch_sd <- function(i, j) sqrt((grp.vars[i] / grp.sizes[i]) + (grp.vars[j] / grp.sizes[j]))
  get_degree_of_freedom <- function(i, j) {
    A <- ((grp.vars[i] / grp.sizes[i]) + (grp.vars[j] / grp.sizes[j]))^2
    B <- ((grp.vars[i] / grp.sizes[i])^2) / (grp.sizes[i] - 1)
    C <- ((grp.vars[j] / grp.sizes[j])^2) / (grp.sizes[j] - 1)
    A / (B + C)
  }

  mean.diff <- stats::pairwise.table(get_mean_diff, levels(g), p.adjust.method = "none") %>% tidy_squared_matrix()
  weltch.sd <- stats::pairwise.table(get_weltch_sd, levels(g), p.adjust.method = "none") %>% tidy_squared_matrix()
  df <- stats::pairwise.table(get_degree_of_freedom, levels(g), p.adjust.method = "none") %>% tidy_squared_matrix()

  t <- abs(mean.diff$value) / weltch.sd$value
  p.raw <- 2 * stats::pt(-t, df = df$value)  # Calculate the raw p-values (two-tailed)
  p.adj <- stats::p.adjust(p.raw, method = p.adjust.method)  # Adjust p-values using specified method

  se <- weltch.sd$value * sqrt(0.5)
  q <- stats::qtukey(p = conf.level, nb.groups, df = df$value)
  conf.high <- mean.diff$value + q * se
  conf.low <- mean.diff$value - q * se

  n1 <- grp.sizes[mean.diff$group1]
  n2 <- grp.sizes[mean.diff$group2]

  results <- mean.diff %>%
    rename(estimate = .data$value) %>%
    mutate(
      conf.low = conf.low, conf.high = conf.high,
      se = se, statistic = t, df = df$value,
      p.value = p.raw,  # Add raw p-values
      p.adj = p.adj  # Add adjusted p-values
    ) %>%
    add_column(n1 = n1, n2 = n2, .after = "group2") %>%
    add_column(.y. = outcome, .before = "group1") %>%
    add_significance("p.adj") %>%
    mutate(method = "Games-Howell")

  results
}
