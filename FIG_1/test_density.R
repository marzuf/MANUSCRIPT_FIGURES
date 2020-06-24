hicds <- "LG1_40kb"
exprds <- "TCGAluad_mutKRAS_mutEGFR"
fcc <-   get(load(file.path("../../v2_Yuanlong_Cancer_HiC_data_TAD_DA//PIPELINE/OUTPUT_FOLDER", hicds, exprds, "8cOnlyFCC_runAllDown", "all_obs_prodSignedRatio.Rdata")))

plot_dt <- data.frame(fcc=fcc)

plot(density(fcc))
ggplot(plot_dt, aes(x=fcc)) + geom_density()
plot(density(fcc, from=min(fcc), to=max(fcc)))
plot(density(fcc, from=-1, to=1))

x_v1 <- density(fcc)$x
x_v2 <- density(fcc, from=-1, to=1)$x
x_v3 <- density(fcc, from=min(fcc), to=max(fcc))$x

y_v1 <- density(fcc)$y
y_v2 <- density(fcc, from=-1, to=1)$y
y_v3 <- density(fcc, from=min(fcc), to=max(fcc))$y

x_range <- range(c(x_v1, x_v2, x_v3))
y_range <- range(c(y_v1, y_v2, y_v3))

plot(NULL,
     xlim=x_range,
     ylim=y_range)
lines(x=x_v1, y=y_v1, col=1)
lines(x=x_v2, y=y_v2, col=2)
lines(x=x_v3, y=y_v3, col=3)