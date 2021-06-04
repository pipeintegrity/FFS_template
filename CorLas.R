library(purrr)
library(dplyr)

# Corlas Functions --------------------------------------------------------

E =30e3 # Young's Modulus (ksi)

sig_y <- 46 # ksi - YS
sig_u <- 1.49*sig_y # ksi - UTS
# UTS ~ 1.49* YS on average if unknown

n <- 0.537 - 0.526 * sig_y / sig_u
# fit equation to strain hardening data
# Do not use the equation given in the paper,
# horrible fit to data and only based on 4 data points

D <- 24 # Diameter (in.)

R <- D/2 # Radius (in.)

# n <- 0.167 # from testing data if available

a <- 0.250 # crack depth (in.)

L <- 10 # crack length (in.)

a/L # depth/Length ratio

c <- L/2 # c = 1/2 * L (in.)

t <- 0.375 # wall thickness (in.)

at <- a/t # depth to wt ratio

# Charpy Data -------------------------------------------------------------


CVN_ss <- 5 #Sub Size Charpy Energy (Ft-lbs)

t_ss_mm <- 10 # Sub size charpy thickness in mm

t_ss_in <- t_ss_mm / 25.4 # convert t_ss from mm to inches

t_fs_in <- 0.394 # Full size charpy in inches

SA_ss <- 0.9 # Shear area in decimal percent

CVN_us_ss <- CVN_ss / (0.9 * SA_ss + 0.1)
# CVN upper shelf sub size (ft-lbs)

CVN_us_fs <- CVN_us_ss * (t_fs_in / t_ss_in)
# CVN upper shelf full size (ft-lbs)


# Free surface factor -----------------------------------------------------

F_sf1 <- a / c + 2 * t / (pi * a) * (1 - a / c) * tan(pi * a / (2 * t))
# for at <= 0.95

F_sf2 <- a / c + (1 - a / c) * (8.515 + (162 / t) * (a / t - 0.95))
# for 0.95 < at < 1.0

F_s <-  ifelse(at <= 0.95, F_sf1, F_sf2)
# choose right equation based on a/t ratio


# Shape Factor ------------------------------------------------------------
Q_f <-
  1.2581 - 0.20589 * (a / (2 * c)) - 11.493 * (a / (2 * c)) ^ 2 + 29.586 *
  (a / (2 * c)) ^ 3 - 23.584 * (a / (2 * c)) ^ 4 #elliptical shape factor

f_3 <-  (3.85 * sqrt (1 / n) * (1 - n) + (pi * n)) * (1 + n) #f3(n) factor
# Shih and Hutchinson equation

eps_y <- 0.005 - sig_y / E #plastic strain @ 0.5% yield

Ks <- (sig_y / (0.005 - sig_y / E) ^ n) #eq(14) in IPC paper
# This relates the plastic strain at an arbitrary stress
# epsilon_p = (S/K_s)^1/n eq(15): IPC paper use eq(15)in J calc in replace of
# The plastic strain variable since it is dependent on S


# J Integral --------------------------------------------------------------
# A conservative lower-bound value for submerged arc welds
# in many structural and low-alloy steels is 500 in-lb/in2

J1c <- 10e-3 * CVN_us_fs#Critical toughness (in-Klb/in2)
# converted to Kips so units are consistent with YS, etc.
# Other correlations suggested are 12CVN/Ac (Ac = area of FS charpy 0.124 in^2)

J1c #print the J1c

# J integral function used to find stress when J1c = J
# J1c is subtracted from J equation since it is being used to find the root
# uniroot function assumes equation is equal to zero to find root

J <-  function(S) {
  Q_f * F_s * a * (S ^ 2 * pi  / E + f_3 * (S / Ks) ^ (1 / n)  * S) - J1c
}

 # Failure stress based on J integral
S1 <- uniroot(J,c(1,sig_y), extendInt = "no")$root
S1
# keep just the root solution (ksi)
# failure stress for J-integral, finds the root when J1c = J, solves for
# stress wider search intervals don't effect the time to find root is no
# practical difference - set search interval between 1 and 100 to find root

# Folias Factor -----------------------------------------------------------
M_T1 <- sqrt(1 + 1.255 * (c / sqrt(R * t)) ^ 2 -
               0.0135 * (c / sqrt(R * t)) ^ 4) #for c^2/Rt <= 25


M_T2 <- 0.064 * (c / sqrt(R * t)) # for c^2/Rt > 25

M_T <- ifelse(c ^ 2 / (R * t) > 25, M_T2, M_T1)
# Select appropriate Folias factor


# Flow stress -------------------------------------------------------------
sig_f1 <- sig_y + 10 # ksi

sig_f2 <- sig_y + 0.5 * (sig_u - sig_y)
#ksi - Paper suggests that Cf = 0.5 is reasonable fit

sig_flow <- min(sig_f1, sig_f2) #flow stress = min between the two

S2 <- sig_flow * ((1 - L * a / (L * t)) / (1 - L * a / (M_T * L * t)))
# failure pressure based on effective area and folias factor

S1 # print S1 as a check
S2 # print S2 as a check

fail <- min(S1, S2) # find min between S1 and S2

paste("The failure Pressure is ",round(fail,1)) # print failure pressure


#Quick plot to show J vs. Sigma

x <- 1:50

eps_p <- (S1 / Ks) ^ (1 / n)

J_x <- 1:50 %>%
  map_dbl( ~ Q_f * F_s * a * (.x ^ 2 * pi / E + f_3 * (.x / Ks) ^ (1 / n) * .x))

j_df <- bind_cols(x= x, J =J_x) #bind into single data frame

#Quick plot
plot(
  J ~ x ,
  j_df,
  main = "J vs. Sigma",
  xlab = "Sigma",
  type = "l",
  lwd = 2,
  col = 'red'
)

F_s
Q_f
S1
f_3
J1c
Ks

j_df %>%
  ggplot(aes(x, J))+
  geom_line(col='red', lwd=1)+
  theme_minimal(14)+
  scale_y_continuous(breaks = seq(0,1.4, length.out=15), limits = c(0,0.4))+
  xlim(0,40)+
  labs(title = "J vs. Sigma",
       x = "Sigma (ksi)",
       y = "J (1,000 lb/in)")+
  theme(plot.margin = margin(0.6,0.6,0.6,0.6,"cm"))
