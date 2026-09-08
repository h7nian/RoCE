files <- list.files("diagnosis/face_probe/validation/chunks", pattern="^aggdiag_rho.*\\.csv$", full.names=TRUE)
raw <- do.call(rbind, lapply(files, read.csv))
cat(sprintf("aggdiag: %d files, %d sims\n", length(files), nrow(raw)))
cat(sprintf("%5s | %7s %7s | %9s %9s | %8s | %7s %7s\n","rho","eta_s1","eta_s2","disc_s1","disc_s2","lambda","bias","cover"))
for (r in sort(unique(raw$rho))) { d<-raw[raw$rho==r,]
  cat(sprintf("%5.1f | %7.3f %7.3f | %9.4f %9.4f | %8.4f | %+7.4f %7.3f\n",
    r, mean(d$eta_s1,na.rm=T), mean(d$eta_s2,na.rm=T),
    mean(abs(d$mu_ts1-d$mu_ot),na.rm=T), mean(abs(d$mu_ts2-d$mu_ot),na.rm=T),
    mean(d$lambda,na.rm=T), mean(d$estimate-d$truth,na.rm=T),
    mean(d$truth>=d$estimate-1.96*d$se & d$truth<=d$estimate+1.96*d$se,na.rm=T))) }
