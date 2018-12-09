library(deSolve)
library(ggplot2)

# system("R CMD SHLIB siir_ode.c")
# https://doi.org/10.1038/s41564-018-0240-5
dyn.load("siir_ode.so")

iarToR0 <- function(IAR, resistantProp = 0) {
  log((1 - resistantProp) / (1 - IAR - resistantProp)) / IAR
}

plauge <- function(pop = 66000000,
                   infA = 1000,
                   infB = 2.4 * 1000000, # http://dx.doi.org/10.1155/2015/504831 enough for all uni students
                   R0A = 2.011797, # spi-m
                   R0B = 2.011797, # spi-m
                   incA = 2,
                   incB = 2,
                   durA = 4,
                   durB = 4,
                   CRA = 0.625, # spi-m RWC
                   CRB = 0.625, # spi-m RWC
                   CFRA = 0.025, # spi-m RWC
                   CFRB = 0.004, # spi-m 1968
                   time_to_B = 14,
                   time_to_rep = 500 # assuming one dose
) {
  Y <- c(S = pop - infA, EA = infA, EB = 0, IA = 0, IB = 0, RA = 0, RB = 0)
  parms <- c(
    bA = R0A / durA, bB = R0B / durB, fA = 1 / incA, fB = 1 / incB, gA = 1 / durA,
    gB = 1 / durB, pop = sum(Y)
  )

  times <- 0:time_to_B
  out <- list()
  out[[1]] <- ode(Y, times,
    func = "derivs", parms = parms,
    jacfunc = "jac", dllname = "siir_ode",
    initfunc = "initmod", nout = 0
  )
  while (times[length(times)] <= 200) {
    times <- times[length(times)]:(times[length(times)] + time_to_rep)
    Y <- out[[length(out)]][nrow(out[[length(out)]]), 2:8]
    out[[length(out)]] <- out[[length(out)]][-nrow(out[[length(out)]]), ]
    infC <- infB * Y[1] / sum(Y)
    Y <- Y + c(-infC, 0, infC, 0, 0, 0, 0)
    out[[length(out) + 1]] <- ode(Y, times,
      func = "derivs", parms = parms,
      jacfunc = "jac", dllname = "siir_ode",
      initfunc = "initmod", nout = 0
    )
  }
  out <- do.call(rbind, out)
  deaths <- sum(out[nrow(out), c(7, 8)] * c(CFRA, CFRB) * c(CRA, CRB))
  cases <- rowSums(t(t(out[, c(5, 6)]) * c(CRA, CRB)))
  peak_time <- out[which.max(cases), 1]
  peak_case <- cases[peak_time]
  c(deaths = deaths, cases = sum(cases), peak_time = peak_time, peak_case = peak_case)
}

senarios <- expand.grid(delay = seq(7, 7 * 8, 1), doses = seq(1.8 * 10^6, 2.8 * 10^6, 10^5))
tablename <- c(colnames(senarios), names(plauge(infB = 0)))
impact <- sapply(1:nrow(senarios), function(i) {
  c(as.numeric(senarios[i, ]), as.numeric(plauge(infB = 0) -
    plauge(time_to_B = senarios[i, "delay"], infB = senarios[i, "doses"])))
})
impact <- as.data.frame(t(impact))
colnames(impact) <- tablename

impact$persentageDeathsAvoided <- 100 * impact$deaths / (plauge(infB = 0)["deaths"])

ggplot(data=impact)+geom_contour(aes(x=delay, y=doses, z = deaths, colour=stat(level))) + scale_colour_gradientn(colours=rainbow(5))
