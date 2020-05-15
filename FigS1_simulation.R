# For Briskin, Wang, & Bartel. The biochemical basis for the cooperative action of microRNAs. PNAS (2020).
# Simulate the rate of approach to equilibrium (k_obs) under different regimes of k_on and k_off.
# Fig. S1
# Peter Y. Wang 2020
# R

library(tidyr)
library(dplyr)
library(broom)
library(nlstools)
library(ggplot2)
library(viridis)
library(ggnewscale)
library(deSolve)
library(FME)
library(gridExtra)

# Set ggplot2 theme
theme1 = theme(
  panel.grid.major  = element_blank(),
  panel.grid.minor  = element_blank(),
  panel.background  = element_blank(),
  axis.line         = element_line(colour = "black"),
  axis.ticks        = element_line(colour = "black"),
  axis.text         = element_text(color = "black", size = 11),
  axis.title        = element_text(color = "black", size = 12),
  legend.background = element_blank(),
  legend.key        = element_blank(),
  strip.background  = element_blank(),
  strip.text.x      = element_text(size = 11, colour = "black"),
  strip.text.y      = element_text(size = 11, colour = "black", angle = 0)
)

# Initialize df
df.simul = data.frame(row.names = F)
df.nls   = data.frame(row.names = F)

# ks to try
# Facets:
koffSet = c(0.1, 0.05, 0.02, 0.01, 0.005)
# Lines:
konSet  = c(1, 0.6, 0.3, 0.2, 0.1, 0.06, 0.03, 0.02, 0.01, 0.006, 0.003, 0.002)

# Initialize conditions
A_T.dil    = 1  # pM
R_T.dil    = 1  # pM (by sites)
dil.factor = 20

# Function to calculate theoretical fraction bound by KD and A_T, assuming AGO:RNA = 1:1, where KD = koff/kon
ththeta = function(KD, AT){(-sqrt(KD) + sqrt(KD + 4*AT))/(sqrt(KD) + sqrt(KD + 4*AT)) }
# Function to calculate concs
concs = function(A_T, R_T, Fb){
  c(
    AR = R_T * Fb,
    A  = A_T * (1-Fb),
    R  = R_T * (1-Fb)
  )
}

# The ODE
DEmodel = function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    dAR =  kon*R*A - koff*AR
    dA  = -kon*R*A + koff*AR
    dR  = -kon*R*A + koff*AR
    list(c(dAR, dA, dR))})}
t = seq(0, 200, by=0.1)

# Exponential setup
eqn = y ~ ff + fd * exp(-x * kobs)
expCoef = function(df, koffv, konv){
  finalTheta = ththeta(koffv/konv, A_T.dil)
  initTheta = ththeta(koffv/konv, A_T.dil * dil.factor)
  nlsFit = df %>%
    group_by(koff, kon) %>%
    do(data.frame(model = tidy(
      nls(
        fracbound ~ ff + fd * exp(-time * kobs),
        algorithm = "port",
        start     = c(ff = finalTheta, fd = initTheta - finalTheta, kobs = koffv + konv * A_T.dil*(1-finalTheta)),
        upper     = c(ff = 1, fd = 1, koff = Inf),
        data      = .)
    ))) %>%
    ungroup() %>%
    transmute(
      koff  = koff,
      kon   = kon,
      coeff = model.term,
      val   = model.estimate) %>%
    pivot_wider(id_cols = c(koff, kon), names_from = coeff, names_sep = "_", values_from = val) %>%
    transmute(
      koff = koff,
      kon  = kon,
      ff   = ff,
      kobs = kobs,
      fi   = ff + fd)
  nlsFit
}

# Running everything per kon-koff pair
for (i in 1:length(koffSet)) {
  for (j in 1:length(konSet)) {

    # Parameters
    koff   = koffSet[i]
    kon    = konSet[j]
    params = c(koff = koff, kon = kon)

    initC = concs(A_T.dil, R_T.dil, ththeta(koff/kon, A_T.dil * dil.factor))
    state = initC

    # ODE simulation
    simul = ode(y = state, times = t, func = DEmodel, parms = params, maxsteps = 100000) %>%
      as.data.frame() %>%
      transmute(koff = koff, kon = kon, time = time, fracbound = (AR)/(AR+R))
    df.simul = df.simul %>%
      rbind(simul)

    # Exp fit
    thisnls = expCoef(simul, koff, kon)
    df.nls = df.nls %>%
      rbind(thisnls)
  }
}

# Make k cols factors for discrete plotting
df.simul = df.simul %>%
  mutate(
    kon  = factor(kon),
    koff = factor(koff))
df.nls = df.nls %>%
  mutate(
    kon.n  = kon,
    koff.n = koff,
    kon    = factor(kon),
    koff   = factor(koff))

# PLOT TIME
gg.simul = df.simul %>%
  ggplot(aes(x = time, y = fracbound, color = kon)) +
  geom_line() +
  facet_wrap("koff", ncol = 5, scales = "free") +
  xlab("Time (min)") +
  ylab("Fraction bound") +
  scale_color_viridis_d(direction = -1, name = "kon (pM-1 min-1)") +
  scale_y_continuous(limits = c(0, 1.05), expand = c(0, 0), breaks = seq(0,1,0.2)) +
  theme1

gg.deviation1 = df.nls %>%
  ggplot(aes(x = koff.n, y = kobs, color = kon)) +
  geom_point(shape = 16) +
  geom_line(aes(group = kon)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  xlab("koff (min-1, log-scale)") +
  ylab("kobs (min-1, log-scale)") +
  scale_color_viridis_d(direction = -1, name = "kon (pM-1 min-1)") +
  scale_x_log10(breaks = c(1, 0.3, 0.1, 0.03, 0.01, 0.003)) +
  scale_y_log10(breaks = c(1, 0.3, 0.1, 0.03, 0.01, 0.003)) +
  theme1

gg.deviation2 = df.nls %>%
  group_by(kon.n) %>%
  ggplot(aes(x = koff.n, y = koff.n/kobs, color = kon)) +
  geom_point(shape = 16) +
  geom_line(aes(group = kon)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  xlab("koff (min-1, log-scale)") +
  ylab("koff/kobs") +
  scale_color_viridis_d(direction = -1, name = "kon (pM-1 min-1)") +
  scale_x_log10(breaks = c(1, 0.3, 0.1, 0.03, 0.01, 0.003)) +
  theme1

out = arrangeGrob(gg.simul, arrangeGrob(gg.deviation1, gg.deviation2, nrow = 1), nrow = 2)
ggsave(file = paste0("kineticModel.pdf"), out, width = 10, height = 8, units = "in")
