library(tidyverse)

set.seed(2021)
N_PLOT <- 120000 # Sample size for plotting ideal distribution

# 1. Sum of normal variables -> Normal Distribution

## Generate Data

n <- 20000
x1 <- rnorm(n = n, mean = 3, sd = 0.6)
x2 <- rnorm(n = n, mean = 1, sd = 1.3)
x3 <- rnorm(n = n, mean = -2, sd = 0.9)

# Show histograms of original distribution

tibble(
    variable = c(rep("x1", n), rep("x2", n), rep("x3", n)),
    value = c(x1, x2, x3)
) %>% 
    ggplot(mapping = aes(x = value, fill = variable)) +
    geom_histogram(color="#e9ecef", alpha=0.6, position = 'identity')

x_sum <- 2 * x1 + 1 * x2 + (-3) * x3

tibble(value = x_sum) %>% 
    ggplot(mapping = aes(x = value)) +
    geom_histogram(color="#e9ecef", alpha=0.6, position = 'identity')

tibble(value = x_sum) %>% 
    ggplot(aes(sample = value)) +
    stat_qq() +
    stat_qq_line(colour = "blue")

# 2. Sum of squared normal variables -> Chi-square Distribution

## Generate Data

n <- 20000
x1 <- rnorm(n = n, mean = 0, sd = 1)
x2 <- rnorm(n = n, mean = 0, sd = 1)
x3 <- rnorm(n = n, mean = 0, sd = 1)

x_sum_sqr <- x1^2 + x2^2 + x3^2

tibble(value = x_sum_sqr) %>% 
    ggplot() +
    geom_histogram(mapping = aes(x = value, y = stat(count / sum(count))), 
                   binwidth = 0.5, color="#e9ecef", alpha=0.6, position = 'identity') +
    geom_freqpoly(data = tibble(value = rchisq(N_PLOT, df = 3)),
                  mapping = aes(x = value, y = stat(count / sum(count))),
                  binwidth = 0.5, size = 1.3, alpha = 0.6, colour = "red") +
    ylab("proportion") +
    coord_cartesian(xlim =c(0, 20))

# 3. Division of two sum of squared normal variables 

## Generate Data

n <- 20000
x1 <- rnorm(n = n, mean = 0, sd = 1)
x2 <- rnorm(n = n, mean = 0, sd = 1)
x3 <- rnorm(n = n, mean = 0, sd = 1)
x4 <- rnorm(n = n, mean = 0, sd = 1)
x5 <- rnorm(n = n, mean = 0, sd = 1)
x6 <- rnorm(n = n, mean = 0, sd = 1)
x7 <- rnorm(n = n, mean = 0, sd = 1)

## Check the two Chi-square distribution first

x_sum_sqr_1 <- x1^2 + x2^2 + x3^2        # degree of freedom = 3
x_sum_sqr_2 <- x4^2 + x5^2 + x6^2 + x7^2 # degree of freedom = 4

tibble(variable = c(rep("u1", n), rep("u2", n)), 
       value = c(x_sum_sqr_1, x_sum_sqr_2)) %>% 
    ggplot() +
    geom_histogram(mapping = aes(x = value, y = stat(count / sum(count) * 2), fill = variable), 
                   binwidth = 0.5, color="#e9ecef", alpha=0.6, position = 'identity') +
    geom_freqpoly(data = tibble(value = rchisq(N_PLOT, df = 3), variable = "u1"),
                  mapping = aes(x = value, y = stat(count / sum(count))),
                  binwidth = 0.5, size = 1.3, alpha = 0.8, colour = "grey70") +
    geom_freqpoly(data = tibble(value = rchisq(N_PLOT, df = 4), variable = "u2"),
                  mapping = aes(x = value, y = stat(count / sum(count))),
                  binwidth = 0.5, size = 1.3, alpha = 0.8, colour = "grey70") +
    facet_grid(~variable) +
    ylab("proportion") +
    coord_cartesian(xlim =c(0, 20))

## After division, they became F distributed

x_sum_sqr_div <- (x_sum_sqr_1 / 3) / (x_sum_sqr_2 / 4)

tibble(value = x_sum_sqr_div) %>% 
    ggplot() +
    geom_histogram(mapping = aes(x = value, y = stat(count / sum(count))), 
                   binwidth = 0.5, color="#e9ecef", alpha=0.6, position = 'identity') +
    geom_freqpoly(data = tibble(value = rf(N_PLOT, df1 = 3, df2 = 4)),
                  mapping = aes(x = value, y = stat(count / sum(count))), 
                  binwidth = 0.5, size = 1.3, alpha = 0.6, colour = "red") +
    ylab("proportion") +
    coord_cartesian(xlim =c(0, 20))

# Proportion testing: Z test vs Chi-square test

## Generate Data

n <- 40000 # sample size for each group
ab_data <- tibble(
    group = c(rep("A", n), rep("B", n)),
    response = c(rbernoulli(n = n, p = 0.4), rbernoulli(n = n, p = 0.405))
)
ab_data_contig <- with(ab_data, table(group, response))
ab_data_contig
ab_data_contig %>% prop.table(margin = 1) %>% round(4)

ab_data_list <- with(
    ab_data_contig %>% 
        as_tibble %>% 
        mutate(n_true = ifelse(response == "TRUE", n, 0)) %>%
        group_by(group) %>%
        summarise(n = sum(n),
                  response = sum(n_true)),
    list(response = response, n = n)
)
ab_data_list

## Built-in testing

### Proportion Z test
with(
    ab_data_list,
    prop.test(x = response, n = n, correct = FALSE)
)
### Chi-square test
with(ab_data, chisq.test(group, response, correct = FALSE))

## Calculate Z test manually

p_pooled <- with(ab_data_list, sum(response) / sum(n))
p_diff <- with(
    ab_data_list, 
    ((response / n) * c(-1, 1)) %>% sum
)
z_stat <- p_diff / sqrt(p_pooled * (1 - p_pooled) * sum(1 / ab_data_list$n))
p_val_z_test <- 2 * pnorm(z_stat, lower.tail = FALSE)
p_val_z_test

## Visualize Z test result

tibble(x = c(-3.5, 3.5)) %>% 
    ggplot(aes(x)) +
        stat_function(fun = dnorm, n = N_PLOT, 
                      args = list(mean = 0,  sd = 1), size = 1) +
        geom_area(stat = "function", fun = dnorm,
                  fill = "red", xlim = c(z_stat, 3.5), alpha = 0.3) +
        geom_area(stat = "function", fun = dnorm,
                  fill = "red", xlim = c(-3.5, -z_stat), alpha = 0.3) +
        labs(y = "density") +
        scale_x_continuous(breaks = seq(-3, 3))

## Calculate Chi-square test manually

chi_stat <- ab_data_contig %>% 
    as_tibble %>%
    mutate(
        multiplier = ifelse(response == "TRUE", p_pooled, 1 - p_pooled),
        expected = get("n", envir = .GlobalEnv) * multiplier,
        chi = (n - expected) ^ 2 / expected
    ) %>%
    with(., sum(chi))
p_val_chi_test <- pchisq(chi_stat, df = 1, lower.tail = FALSE)
p_val_chi_test

## Visualize Chi-square test result

tibble(x = c(0, 5)) %>% 
    ggplot(aes(x)) +
        stat_function(fun = dchisq, n = N_PLOT, 
                      args = list(df = 1), size = 1) +
        geom_area(stat = "function", fun = dchisq, args = list(df = 1),
                  fill = "red", xlim = c(z_stat, 5), alpha = 0.3) +
        labs(y = "density") +
        scale_x_continuous(breaks = seq(0, 5)) +
        coord_cartesian(ylim = c(0, 2))

