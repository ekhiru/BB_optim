cegores <- NULL
for(i in 1:30) {
  cegores <- rbind(cegores,
                   readRDS(paste0("cegores-er0-r", i, ".rds")),
                   readRDS(paste0("cegores-er1-r", i, ".rds")))
}
saveRDS(cegores, file="cegores.rds")
write.csv(cegores, file="cegores.csv", row.names=FALSE)
aggregate(cegores$time, list(er=cegores$eval_ranks), mean)
head(cegores[cegores$eval_ranks == 0,])
head(cegores[cegores$eval_ranks == 1,])

best_fitness <- aggregate(cegores$Fitnes, cegores[c("eval_ranks", "seed")], min)

library(ggplot2)
best_fitness$eval_ranks <- factor(best_fitness$eval_ranks)
gg <- ggplot(best_fitness, aes(x=eval_ranks, group=eval_ranks, y=x)) +
  geom_boxplot() + geom_jitter()

ggsave(gg, file="boxplot.pdf")

aggregate(best_fitness$x, best_fitness["eval_ranks"], median)
