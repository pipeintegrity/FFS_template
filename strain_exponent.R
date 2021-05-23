
# CorLas Equations 05202021 -----------------------------------------------

## These equations were taken from a published paper "Applicability of Existing
## Fracture Initiation Models to Modern Line Pipe Steels"

# All units are ksi, inches and in^2
library(readxl)
library(tidyverse)
library(tidymodels)
# library(ggmisc)
theme_set(theme_minimal(14,"serif"))


tensile <-
  read_excel(path = "C:\\Users\\Joel\\OneDrive - RSI Pipeline Solutions\\PGE\\charpy\\MasterDB-SQL-2020-09-24.xlsx",
             sheet = "Tensile") %>%
  janitor::clean_names() %>%
  # select(x0_5_percent_eul_ys_ksi, uts_ksi, n) %>%
  rename(ys = x0_5_percent_eul_ys_ksi, uts = uts_ksi) %>%
  group_by(group, feature) %>%
  summarise(
    ys = mean(ys, na.rm = T),
    uts = mean(uts, na.rm = T),
    n = mean(n, na.rm = T)
  ) %>%
  select(group, feature, ys, uts, n) %>%
  drop_na() %>%
  mutate(ys_uts = ys/uts)

comp <- read_excel(path = "C:\\Users\\Joel\\OneDrive - RSI Pipeline Solutions\\PGE\\charpy\\MasterDB-SQL-2020-09-24.xlsx",
                   sheet = "Composition") %>%
  janitor::clean_names() %>%
  mutate(across(fe:pb, as.numeric)) %>%
  group_by(group, feature) %>%
  summarise(c = mean(c, na.rm = T), mn = mean(mn, na.rm=T))

join <- left_join(tensile, comp, by=c("group", "feature"))

tensile %>%
  ggplot(aes(ys / uts, n)) +
  geom_point(col = 'orangered', alpha = 0.5) +
  geom_smooth(method = "lm", se = F) +
  theme_minimal(14, "serif") +
  labs(title = "Strain Hardening Exponent to YS/UTS",
       x = "YS/UTS")

n_mod <- lm(n ~ ys/uts, data = tensile)

n_mod2 <- lm(n ~ ys+uts, data = tensile)

n_mod3 <- lm(n ~ ys_uts, data = tensile)

n_mod4 <- lm(n ~ ys*uts, data = tensile)

mod4_fit <- predict(n_mod4)



g1 <- glance(n_mod)
g2 <- glance(n_mod2)
g3 <- glance(n_mod3)
g4 <- glance(n_mod4)

bind_rows(g1, g2, g3, g4)

aug1 <- augment(n_mod)
aug2 <- augment(n_mod2)
aug3 <- augment(n_mod3)
aug4 <- augment(n_mod4)



ssq1 <- sum(aug1$.resid^2)
ssq2 <- sum(aug3$.resid^2)

ssq1
ssq2

lm_eqn <- function(df) {
  n_mod3 <- lm(n ~ ys_uts, data = df)
  eq <-
    substitute(
      italic(n) == a ~ b %.% italic(YS/UTS) * "," ~  ~ italic(r) ^2 ~ "=" ~ r2,
      list(
        a = format(unname(coef(n_mod3)[1]), digits = 3),
        b = format(unname(coef(n_mod3)[2]), digits = 3),
        c = format(unname(coef(n_mod3)[3]), digits = 3),
        r2 = format(summary(n_mod3)$r.squared, digits = 2)
      )
    )
  as.character(as.expression(eq))

}


n_mod3 %>%
  ggplot(aes(.fitted, n)) +
  geom_point(alpha = 0.6, col = 'sienna1') +
  # geom_smooth(method = "lm", se = F) +
  labs(title = "Observed to Predicted Strain Hardening Exponent (n)",
       x = "Predicted n",
       y = "Observed n") +
  annotate(
    "text",
    x = 0.1,
    y = 0.25,
    label = lm_eqn(tensile),
    parse = TRUE,
    size = 5
  )+
  geom_abline(lty=2, col='grey50')+
  coord_obs_pred()



n_mod %>%
  ggplot(aes(.fitted, n))+
  geom_point(alpha=0.5, col='red')+
  geom_smooth(method = "lm", se=F)+
  labs(title = "Mod1")


# tidymodels  -------------------------------------------------------------

lm_model <-
  linear_reg() %>%
  set_mode("regression") %>%
  set_engine("lm")

## More consistent results in double dipping when center & scale is used
## Doesn't change the R^2 but more consistent grade predictions later on
n_rec <- recipe(n ~ ys_uts ,
                      data = tensile) %>%
  step_center(all_predictors()) %>%
  step_scale(all_predictors()) %>%
  # step_interact(terms = ~ mn:c) %>%
  prep() #retain data =TRUE is default

n_bake <- bake(n_rec, new_data = NULL)

fitted <-  fit(lm_model,
                formula = n ~ ys_uts ,
                data = n_bake) # fitted model object)

n_workflow <- workflow() %>%
  add_model(lm_model) %>%
  add_recipe(n_rec)

fit_n <-  n_workflow %>%
  fit(n_bake)

glance(fit_n)


# quadratic ---------------------------------------------------------------

mod_p <- lm(n~poly(ys_uts,2), data = tensile)
tidy(mod_p)
aug_p <- augment(mod_p, interval = "prediction")
fit_p <- predict(mod_p)

aug <- augment(fitted$fit, interval = "prediction") %>%
  mutate(fit_p = (-0.00546+0.556*tensile$ys_uts-0.547*tensile$ys_uts^2)) %>%
  relocate(fit_p,.after = .fitted) %>%
  rename(lm_model=.fitted, paper=fit_p)

aug_long <- aug %>%
  pivot_longer(.lower:.upper)

aug %>%
  pivot_longer(lm_model:paper) %>%
  ggplot(aes(value, n)) +
  geom_point(aes(col=name), alpha = 0.6) +
  geom_smooth(aes(col=name),
              method = "lm",
              lwd = 0.8,
              se=F)+
  coord_obs_pred()+
  geom_abline(col='grey50', lty=2)



tidy(fit_n)

aug %>%
  ggplot(aes(.pred, n)) +
  geom_point(alpha = 0.6, col = 'sienna1') +
  # geom_smooth(method = "lm", se = F) +
  labs(title = "Observed to Predicted Strain Hardening Exponent (n)",
       x = "Predicted n",
       y = "Observed n") +
  # annotate(
  #   "text",
  #   x = 0.1,
  #   y = 0.25,
  #   label = lm_eqn(tensile),
  #   parse = TRUE,
  #   size = 5
  # )+
  geom_abline(lty=2, col='grey50')+
  coord_obs_pred()


# Corlas Functions --------------------------------------------------------

sig_y <- 43.4 #ksi - YS
sig_u <- 65 #ksi - UTS

D <- 24

R <- D/2

n <- 0.167 # from testing data

a <- 0.10

L <- 2.5

c <- L/2
t <- 0.313

c^2/(R*t)
at <- a/t

F_sf1 <- a / c + 2 * t / (pi * a) * (1 - a / c) * tan(pi * a / (2 * t))

F_sf2 <- a / c + (1 - a / c) * (8.515 + (162 / t) * (a / t - 0.95))

Q_f <-
  1.2581 - 0.20589 * (a / (2 * c)) - 11.493 * (a / (2 * c)) ^ 2 + 29.586 *
  (a / (2 * c)) ^ 3 - 23.584 * (a / (2 * c)) ^ 4

f_3 <- (3.85 * (1 - 1 / n) * sqrt(n) + (pi / n)) * (1 + 1 / n) #f3_1/n factor

sig_f <- sig_y+10 #ksi
sig_f2 <- sig_y+0.5*(sig_u-sig_y)#ksi

M_T1 <- sqrt(1 + 1.255 * (c / sqrt(R * t)) ^ 2 -
               0.0135 * (c / sqrt(R * t)) ^ 4) #for c^2/Rt <= 25


M_T2 <- 0.064 * (c / sqrt(R * t)) #for c^2/Rt > 25

M_T <- ifelse(c ^ 2 / (R * t) > 25, M_T2, M_T1)
