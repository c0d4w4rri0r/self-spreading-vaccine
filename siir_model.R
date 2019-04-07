library(deSolve)
library(ggplot2)
library(data.table)

# system("R CMD SHLIB siir_ode.c")
# https://doi.org/10.1038/s41564-018-0240-5
dyn.load("siir_ode.so")

iarToR0 <- function(IAR, resistantProp = 0) {
  log((1 - resistantProp) / (1 - IAR - resistantProp)) / IAR
}

plauge <- function(pop = 66900000,
                   infA = 1,
                   infB = 2.4 * 1000000, # http://dx.doi.org/10.1155/2015/504831 enough for all uni students
                   R0A = 2.011797, # spi-m
                   R0B = R0A, # spi-m
                   incA = 0.64, # http://dx.doi.org/10.1016/j.epidem.2012.06.001 # used to be 0.8 no idea where from
                   incB = incA, # http://dx.doi.org/10.1016/j.epidem.2012.06.001
                   durA = 1.27, # http://dx.doi.org/10.1016/j.epidem.2012.06.001 used to be 1.8 no idea where from
                   durB = durA, # http://dx.doi.org/10.1016/j.epidem.2012.06.001
                   CRA = 0.625, # spi-m RWC
                   CRB = CRA, # spi-m RWC
                   CFRA = 0.025, # spi-m RWC
                   CFRB = 0.004, # spi-m 1968
                   time_to_B = 14,
                   deaths_till_noticed = 1,
                   outputCurves = FALSE
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
  times <- times[length(times)]:250
  Y <- out[[1]][nrow(out[[1]]), 2:8]
  out[[1]] <- out[[1]][-nrow(out[[1]]), ]
  infC <- infB * Y[1] / sum(Y)
  Y <- Y + c(-infC, 0, infC, 0, 0, 0, 0)
  out[[2]] <- ode(Y, times,
                  func = "derivs", parms = parms,
                  jacfunc = "jac", dllname = "siir_ode",
                  initfunc = "initmod", nout = 0
  )
  out <- do.call(rbind, out)
  deaths <- sum(out[nrow(out), c(7, 8)] * c(CFRA, CFRB) * c(CRA, CRB))
  endCases <- sum(out[nrow(out), c(7, 8)] * c(CRA, CRB))
  cases <- rowSums(t(t(out[, c(5, 6)]) * c(CRA, CRB)))
  time_noticed<-which.min(abs(deaths_till_noticed-rowSums(out[, c(5, 7)]) * CFRA * CRA))
  if (infB!=0) time_noticed<-0
  peak_time <- out[which.max(cases), 1]
  peak_case <- cases[which.max(cases)]
  if (outputCurves) return(out)
  else return(c(deaths = deaths, cases = endCases, peak_time = peak_time,
                peak_case = peak_case, time_noticed = time_noticed))
}

senarios <- expand.grid(
  delay = seq(7, 7 * 20, 1),
  doses = seq(1.8 * 10^6, 2.8 * 10^6, 10^5),
  latentPeriod = seq(0.64, 2, (2-0.64)/3),
  infectiousPeriod = seq(1.27, 4, (4-1.27)/3)
)
tablename <- c(colnames(senarios), names(plauge(infB = 0)))
impact <- sapply(1:nrow(senarios), function(i) {
  c(as.numeric(senarios[i, ]), as.numeric(plauge(
    infB = 0, incA = senarios[i, "latentPeriod"],
    durA = senarios[i, "infectiousPeriod"]
  ) -
    plauge(time_to_B = senarios[i, "delay"], infB = senarios[i, "doses"],
  incA = senarios[i, "latentPeriod"], durA = senarios[i, "infectiousPeriod"])
  ))
})
impact <- as.data.frame(t(impact))
colnames(impact) <- tablename

impact$delay<-impact$delay-impact$time_noticed
impact$latentPeriod<-as.factor(impact$latentPeriod)
impact$infectiousPeriod<-as.factor(impact$infectiousPeriod)

refDeaths<-plauge(
  infB = 0, incA = senarios[1, "latentPeriod"],
  durA = senarios[1, "infectiousPeriod"]
)["deaths"]

#impact$persentageDeathsAvoided <- 100 * impact$deaths / (plauge(infB = 0)["deaths"])

graphWidth<-8.3
graphHeight<-5.3

#contours of delay v doses for deaths
pdf("output/doses-delay-contours.pdf",width = graphWidth,height=graphHeight)
ggplot(data = impact[impact$latentPeriod==0.64 & impact$infectiousPeriod==1.27,c("delay","doses","deaths")]) + geom_contour(aes(x = delay, y = doses, z = deaths, colour = stat(level))) + scale_colour_gradientn(colours = rainbow(5), name="deaths\nprevented")
dev.off()
#contours of delay v doses for peek cases
pdf("output/doses-delay-extra-cases.pdf",width = graphWidth,height=graphHeight)
ggplot(data = impact[impact$latentPeriod==0.64 & impact$infectiousPeriod==1.27,]) + geom_contour(aes(x = delay, y = doses, z = -peak_case, colour = stat(level))) + scale_colour_gradientn(colours = rainbow(5), name="Extra\ncases\nat peek")
dev.off()
#contours of delay v doses for peek time
pdf("output/doses-delay-contours-shift.pdf",width = graphWidth,height=graphHeight)
ggplot(data = impact[impact$latentPeriod==0.64 & impact$infectiousPeriod==1.27,]) + geom_contour(aes(x = delay, y = doses, z = peak_time.time, colour = stat(level))) + scale_colour_gradientn(colours = rainbow(5), name="peek\nshift\n(days)")
dev.off()
# deaths by delay graph for 3 difrent does levels
pdf("output/deaths-delay-curves-by-doses.pdf",width = graphWidth,height=graphHeight)
ggplot(data = impact[impact$doses %in% (c(1.8, 2.4, 2.8) * 10^6) &
                       impact$latentPeriod == 0.64 &
                       impact$infectiousPeriod == 1.27 &
                       impact$delay>=0, ]) + 
  geom_path(aes(delay, deaths, colour = factor(doses / 10^6),
                group=factor(doses / 10^6))) +
  labs(color = "Doses") +
  xlab("Days till vaccination") +
  ylab("Deaths averted") +
  scale_x_continuous(limits=c(0,50)) +
  scale_y_continuous(sec.axis = sec_axis(~./refDeaths, labels = scales::percent)) +
  scale_colour_discrete(name="Doses", breaks=c("1.8","2.4","2.8"), labels=c("1.8m","2.4m","2.8m"))
dev.off()
# deaths by delay graph for one does but full range of latent / infectious periods
impact2<-impact[impact$doses==2.4 * 10^6, ]
impact2<-impact2[impact2$infectiousPeriod %in% c("1.27","4"),]
impact2<-impact2[impact2$latentPeriod %in% c("0.64","2"),]
pdf("output/deaths-delay-curves-by-periods.pdf",width = graphWidth,height=graphHeight)
ggplot(data = impact2[impact2$delay>=0,]) + 
  geom_path(aes(delay, deaths, colour = infectiousPeriod, linetype = latentPeriod,
                group=interaction(latentPeriod,infectiousPeriod))) +
  scale_color_discrete(name="Infectious\nPeriod", breaks=c("1.27","4"), labels=c("1.27 days","4 days")) +
  scale_linetype_discrete(name="Latent\nPeriod", breaks=c("0.64","2"), labels=c("0.64 days","2 days")) +
  scale_y_continuous(name="deaths averted", sec.axis = sec_axis(~./refDeaths, labels = scales::percent)) +
  scale_x_continuous(limits=c(0,120))
dev.off()

impact<-as.data.table(impact)

impact[,lapply(.SD, sd),by=.(delay,doses),.SDcols=names(plauge(infB = 0))]
