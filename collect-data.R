cegores <- NULL
for(i in 1:30) {
  cegores <- rbind(cegores,
                   readRDS(paste0("cegores-er0-r", i, ".rds")),
                   readRDS(paste0("cegores-er1-r", i, ".rds")))
}
saveRDS(cegores, file="cegores.rds")
write.csv(cegores, file="cegores.csv", row.names=FALSE)
aggregate(cegores$time, list(er=cegores$eval_ranks), mean)
print(cegores)
