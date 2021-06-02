
# CorLas Equations 05202021 -----------------------------------------------

## These equations were taken from a published paper "Applicability of Existing
## Fracture Initiation Models to Modern Line Pipe Steels"

# All units are ksi, inches and in^2
library(readxl)
library(tidyverse)
library(tidymodels)
library(patchwork)
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
    n = mean(n, na.rm = T),
    .groups = "drop"
  ) %>%
  select(group, feature, ys, uts, n) %>%
  drop_na() %>%
  mutate(ys_uts = ys / uts)

comp <- read_excel(path = "C:\\Users\\Joel\\OneDrive - RSI Pipeline Solutions\\PGE\\charpy\\MasterDB-SQL-2020-09-24.xlsx",
                   sheet = "Composition") %>%
  janitor::clean_names() %>%
  mutate(across(fe:pb, as.numeric)) %>%
  group_by(group, feature) %>%
  summarise(c = mean(c, na.rm = T), mn = mean(mn, na.rm=T), .groups = "drop")

join <- left_join(tensile, comp, by=c("group", "feature"))

tensile %>%
  ggplot(aes(uts, ys)) +
  geom_point(col = 'orangered', alpha = 0.5) +
  geom_smooth(method = "lm", se = F) +
  theme_minimal(14, "serif") +
  labs(title = "Strain Hardening Exponent to YS/UTS",
       x = "YS/UTS")+
  scale_x_log10()+scale_y_log10()

uts_ys_mod <- lm(uts ~ ys, tensile)
summary(uts_ys_mod)


n_mod <- lm(n ~ ys/uts, data = tensile)

n_mod2 <- lm(n ~ ys + uts, data = tensile)

n_mod3 <- lm(n ~ ys_uts, data = tensile)

n_mod4 <- lm(n ~ ys*uts, data = tensile) #interaction not significant

mod4_fit <- predict(n_mod4)


mod2_chk <- performance::check_model(n_mod2)
mod3_chk <- performance::check_model(n_mod3)

mod2_chk
mod3_chk


g1 <- glance(n_mod)
g2 <- glance(n_mod2)
g3 <- glance(n_mod3)
g4 <- glance(n_mod4)

bind_rows(g1, g2, g3, g4)

aug1 <- augment(n_mod)
aug2 <- augment(n_mod2)%>% mutate(model = "Mod2")
aug3 <- augment(n_mod3) %>% mutate(model = "Mod3")
aug4 <- augment(n_mod4)%>% mutate(model = "Mod4")

bind_rows(aug2, aug3) %>%
  ggplot(aes(.fitted, n)) +
  geom_point(alpha = 0.6, aes(col=model)) +
  # geom_smooth(method = "lm", se = F) +
  labs(title = "Observed to Predicted Strain Hardening Exponent (n)",
       x = "Predicted n",
       y = "Observed n") +
  geom_abline(lty=2, col='grey50')+
  coord_obs_pred()




ssq1 <- sum(aug4$.resid^2)
ssq2 <- sum(aug3$.resid^2)

ssq1
ssq2

lm_eqn <- function(df) {
  n_mod3 <- lm(n ~ ys_uts, data = df)
  eq <-
    substitute(
      italic(n) == a ~ b %.% ~ italic(YS/UTS) * "," ~  ~ italic(r) ^2 ~ "=" ~ r2,
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
  # step_log(all_predictors()) %>%
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
  geom_abline(col='grey50', lty=2)+
  labs(title = "Linear Model and IPC Paper correlations for n",
       x = "Predicted n",
       y =" Observed n")



tidy(fit_n)

aug %>%
  ggplot(aes(lm_model, n)) +
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

E =30e3 # Young's Modulus (ksi)

sig_y <- 43.4 # ksi - YS
sig_u <- 65 # ksi - UTS

# n <- 0.537 - 0.526 * sig_y / sig_u #fit equation to strain hardening data
#Do not use the equation given in the paper, horrible fit to data and only based
#on 4 data points

D <- 24 # Diameter (in.)

R <- D/2 # Radius (in.)

n <- 0.167 # from testing data

a <- 0.15 # crack depth (in.)

L <- 15 # crack length (in.)

a/L # depth/Length ratio

c <- L/2 # c = 1/2 * L (in.)

t <- 0.313 # wall thickness (in.)

at <- a/t # depth to wt ratio

F_sf1 <- 2*t/(pi*a) * tan(pi * a / (2 * t))*(1-2*a/L) # for a/t <= 0.95

F_sf2 <- a / c + (1 - a / c) * (8.515 + (162 / t) * (a / t - 0.95)) # for 0.95 < a/t < 1.0

F_s <-  ifelse(at <= 0.95, F_sf1, F_sf2)
# choose right equation based on a/t ratio

Q_f <-
  1.2581 - 0.20589 * (a / (2 * c)) - 11.493 * (a / (2 * c)) ^ 2 + 29.586 *
  (a / (2 * c)) ^ 3 - 23.584 * (a / (2 * c)) ^ 4 #elliptical shape factor

f_3 <-  (3.85 * sqrt (1 / n) * (1 - n) + (pi * n)) * (1 + n) #f3(n) factor
# Shih and Hutchinson equation

eps_p <- 0.005 - sig_y / E #plastic strain @ yield


# Charpy Data -------------------------------------------------------------

CVN_ss <- 10 #Sub Size Charpy Energy (Ft-lbs)

t_ss <- 6.67 # Subsize charpy thickness in mm
t_fs <- 0.394 # Full size charpy in inches
SA_ss <- 0.99 # Shear area in decimal percent


CVN_us_ss <- CVN_ss/(0.9*SA_ss+.1) # CVN upper shelf sub size

CVN_us_fs <- CVN_us_ss*(t_fs/(t_ss/25.4)) # CVN upper shelf full size

# J Integral --------------------------------------------------------------
# A conservative lower-bound value for submerged arc welds
# in many structural and low-alloy steels is 500 in-lb/in2

J1c <- 10e-3* CVN_us_fs #Critical toughness (in-Klb/in2)
# converted to Kips so units are consistent with YS, etc.

J1c #print the J1c
# Other correlations suggested are 12CVN/Ac (Ac = area of FS charpy 0.124 in^2)


# J integral function used to find stress when Jc = J
J <-  function(S) {
  Q_f * F_s * a*(S ^ 2 * pi  / E + f_3 * n * eps_p * S) - J1c
}

S1 <- uniroot(J,c(1,100), extendInt = "yes")$root
S1

# keep just the root solution (ksi)
# failure stress for J-integral, finds the root when J1c = J, solves for stress
# wider search intervals don't effect the time to find root in any practical
# way - set search interval between 1 and 100 to find root

# Folias Factor -----------------------------------------------------------
M_T1 <- sqrt(1 + 1.255 * (c / sqrt(R * t)) ^ 2 -
               0.0135 * (c / sqrt(R * t)) ^ 4) #for c^2/Rt <= 25


M_T2 <- 0.064 * (c / sqrt(R * t)) # for c^2/Rt > 25

M_T <- ifelse(c ^ 2 / (R * t) > 25, M_T2, M_T1) # Select appropriate Folias factor


# Flow stress -------------------------------------------------------------


sig_f1 <- sig_y + 10 # ksi

sig_f2 <- sig_y + 0.5 * (sig_u - sig_y)
#ksi - Paper suggests that Cf = 0.5 is reasonable fit

sig_flow <- min(sig_f1, sig_f2) #flow stress = min between the two

S2 <- sig_flow * ((1 - L * a / (L * t)) / (1 - L * a / (M_T * L * t)))
# failure pressure based on effective area and folias factor

S1 # print S1 as a check
S2 # print S2 as a check

fail <- min(S1, S2) #find min between S1 and S2

fail # print failure pressure


#Quick plot to show J vs. Sigma
x <- 1:100
J_x <- 1:100 %>%  map_dbl( ~ Q_f * F_s * (.x ^ 2 * pi * a / E + f_3 * n * a * eps_p * .x))

 j_df <- bind_cols(x=x, J =J_x) # bind into single data frame

 #Quick plot for sanity check
 plot(J ~ x , j_df, main = "J vs. Sigma", xlab = "Sigma", type="l", lwd=2, col='red')
