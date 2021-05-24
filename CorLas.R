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

F_sf1 <- 2*t/(pi*a) * tan(pi * a / (2 * t))*(1-2*a/L) # for at <= 0.95

F_sf2 <- a / c + (1 - a / c) * (8.515 + (162 / t) * (a / t - 0.95)) # for 0.95 < at < 1.0

F_s <-  ifelse(at <= 0.95, F_sf1, F_sf2)
# choose right equation based on a/t ratio

Q_f <-
  1.2581 - 0.20589 * (a / (2 * c)) - 11.493 * (a / (2 * c)) ^ 2 + 29.586 *
  (a / (2 * c)) ^ 3 - 23.584 * (a / (2 * c)) ^ 4 #elliptical shape factor

f_3 <-  (3.85 * sqrt (1 / n) * (1 - n) + (pi * n)) * (1 + n) #f3(n) factor
# Shih and Hutchinson equation

eps_p <- 0.005 - sig_y / E #plastic strain @ yield


# Charpy Data -------------------------------------------------------------

CVN_ss <- 20 #Sub Size Charpy Energy (Ft-lbs)

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
library(purrr)

x <- 1:100
J_x <- 1:100 %>%  map_dbl( ~ Q_f * F_s * (.x ^ 2 * pi * a / E + f_3 * n * a * eps_p * .x) )

j_df <- bind_cols(x= x, J =J_x) #bind into single data frame

#Quick plot
plot(J ~ x , j_df, main = "J vs. Sigma", xlab = "Sigma", type="l", lwd=2, col='red')

