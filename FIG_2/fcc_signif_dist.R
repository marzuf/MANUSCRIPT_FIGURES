
# Rscript fcc_signif_dist.R

plotType <- "svg"
source("../settings.R")


SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")

plotType <- "svg"
myHeight <- ifelse(plotType == "png", 480 , 7)
myWidth <- ifelse(plotType == "png", 600, 10)
plotCex <- 1.4

require(flux)
require(foreach)
require(doMC)
registerDoMC(nCpu)
require(reshape2)
require(ggpubr)
require(ggsci)

library(ggridges)

ggsci_pal <- "aaas"
ggsci_subpal <- ""

tadSignifThresh <- 0.01


mycols <- pal_uchicago()(5)[4:5]

outFolder <- file.path("FCC_SIGNIF_DIST")
dir.create(outFolder, recursive=TRUE)

options(scipen=100)

startTime <- Sys.time()

myHeightGG <- 7
myWidthGG <- 9

buildData <- TRUE

hicds <- "ENCSR489OCU_NCI-H460_40kb"
exprds <- "TCGAlusc_norm_lusc"
# all_obs_hicds=all_obs_hicds[1]
if(buildData){
  
  all_dt <- foreach(hicds = all_obs_hicds, .combine='rbind') %dopar% {
    exprds = all_exprds[[paste0(hicds)]][1]
    hicds_dt <- foreach(exprds = all_obs_exprds[[paste0(hicds)]], .combine='rbind') %do% {
      
      tad_fcc <- get(load(file.path(pipFolder, hicds, exprds, "8cOnlyFCC_runAllDown", "all_obs_prodSignedRatio.Rdata")))
      
      tad_pval <- get(load(file.path(pipFolder, hicds, exprds, "11sameNbr_runEmpPvalCombined", "emp_pval_combined.Rdata")))
      tad_adjPval <- p.adjust(tad_pval, method="BH")
      
      stopifnot(setequal(names(tad_fcc), names(tad_adjPval)))
      
      tad_adjPval <- tad_adjPval[names(tad_fcc)]
      
      stopifnot(names(tad_adjPval)==names(tad_fcc))
      
      data.frame(
        hicds = hicds, 
        exprds=exprds,
        region = names(tad_fcc),
        FCC = as.numeric(tad_fcc), 
        adjPval = as.numeric(tad_adjPval),
        stringsAsFactors = FALSE
      )
    }
  }  
  outFile <- file.path(outFolder, paste0("all_dt.Rdata"))
  auc_fract_file <- outFile
  save(all_dt, file = outFile, version=2)
  cat(paste0("... written: ", outFile, "\n"))
} else {
  inFile <- file.path(outFolder, paste0("all_dt.Rdata"))
  all_dt <- get(load(inFile))
  # load("FCC_SIGNIF_DIST/all_dt.Rdata")
}

all_tads_dt <- all_dt
all_tads_dt$tad_label <- "all"
signif_tads_dt <- all_dt[all_dt$adjPval <= tadSignifThresh,] 
signif_tads_dt$tad_label <- "signif. only"

plot_dt <- rbind(all_tads_dt, signif_tads_dt)

###### FCC distribution - signif vs. not signif.

plot_dt2 <- all_dt
plot_dt2$tad_label <- ifelse(all_dt$adjPval <= tadSignifThresh, "signif.", "not signif.")

plotTit <- "Distribution of TAD FCC"
subTit <- paste0("signif. TADs: adj. p-val <= ", tadSignifThresh)

dist_p <- ggplot(plot_dt2, aes(x=FCC, fill=tad_label, col=tad_label)) + 
  ggtitle(plotTit, subtitle = subTit)+
  labs(fill="TADs", x = paste0("TAD FCC"), y ="Density") + 
  guides(color=FALSE)+
  # geom_histogram(aes(y=..density..), alpha=0.5, 
  #                position="identity")+
  geom_density(alpha=.2) +
  # geom_histogram(aes(y=..density..))+
  # geom_density(alpha=.3)+
  eval(parse(text=paste0("scale_fill_", ggsci_pal, "(", ggsci_subpal, ")"))) + 
  eval(parse(text=paste0("scale_color_", ggsci_pal, "(", ggsci_subpal, ")"))) + 
  scale_x_continuous(breaks = scales::pretty_breaks(n = 8))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 8))+
  my_box_theme+
  theme(
    axis.line = element_line(),
    axis.text.x = element_text(size=14),
    axis.text.y = element_text(size=14),
    legend.text =element_text(size=12)
  )

outFile <- file.path(outFolder, paste0("FCC_signifNot_density.", plotType))
ggsave(dist_p, filename = outFile, height=myHeightGG, width=myWidthGG)
cat(paste0("... written: ", outFile, "\n"))



ridge1_p <- ggplot(plot_dt2, aes(x = FCC, y = tad_label, fill = tad_label)) +
  ggtitle(plotTit, subtitle = subTit)+
  labs(fill="TADs", x = paste0(" TAD FCC"), y="") + 
  geom_density_ridges_gradient(scale = 0.95, size = 0.1, rel_min_height = 0.001) +
  eval(parse(text=paste0("scale_fill_", ggsci_pal, "(", ggsci_subpal, ")"))) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 8))+
  my_box_theme+
  theme(
    axis.line = element_line(),
    axis.text.x = element_text(size=14),
    axis.text.y = element_text(size=14),
    legend.text =element_text(size=12)
  )
outFile <- file.path(outFolder, paste0("FCC_signifNot_density_ggridge_v1.", plotType))
ggsave(ridge1_p, filename = outFile, height=myHeightGG, width=myWidthGG)
cat(paste0("... written: ", outFile, "\n"))


ridge2_p <- ggplot(plot_dt2, aes(x = FCC, y = tad_label, fill = tad_label)) +
  ggtitle(plotTit, subtitle = subTit)+
  labs(fill="TADs", x = paste0("TAD FCC"), y="") + 
  geom_density_ridges_gradient(scale = 1.5, size = 0.1, rel_min_height = 0.001) +
  eval(parse(text=paste0("scale_fill_", ggsci_pal, "(", ggsci_subpal, ", alpha=0.7)"))) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 8))+
  my_box_theme+
  theme(
    axis.line = element_line(),
    axis.text.x = element_text(size=14),
    axis.text.y = element_text(size=14),
    legend.text =element_text(size=12)
  )


outFile <- file.path(outFolder, paste0("FCC_signifNot_density_ggridge_v2.", plotType))
ggsave(ridge2_p, filename = outFile, height=myHeightGG, width=myWidthGG)
cat(paste0("... written: ", outFile, "\n"))












###### FCC distribution - signif vs. all

plotTit <- "Distribution of TAD FCC"
subTit <- paste0("signif. TADs: adj. p-val <= ", tadSignifThresh)

dist_p <- ggplot(plot_dt, aes(x=FCC, fill=tad_label, col=tad_label)) + 
  ggtitle(plotTit, subtitle = subTit)+
  labs(fill="TADs", x = paste0("TAD FCC"), y ="Density") + 
  guides(color=FALSE)+
  # geom_histogram(aes(y=..density..), alpha=0.5, 
  #                position="identity")+
  geom_density(alpha=.2) +
  # geom_histogram(aes(y=..density..))+
  # geom_density(alpha=.3)+
  eval(parse(text=paste0("scale_fill_", ggsci_pal, "(", ggsci_subpal, ")"))) + 
  eval(parse(text=paste0("scale_color_", ggsci_pal, "(", ggsci_subpal, ")"))) + 
  scale_x_continuous(breaks = scales::pretty_breaks(n = 8))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 8))+
  my_box_theme+
  theme(
    axis.line = element_line(),
    axis.text.x = element_text(size=14),
    axis.text.y = element_text(size=14),
    legend.text =element_text(size=12)
  )

outFile <- file.path(outFolder, paste0("FCC_signifAll_density.", plotType))
ggsave(dist_p, filename = outFile, height=myHeightGG, width=myWidthGG)
cat(paste0("... written: ", outFile, "\n"))



ridge1_p <- ggplot(plot_dt, aes(x = FCC, y = tad_label, fill = tad_label)) +
  ggtitle(plotTit, subtitle = subTit)+
  labs(fill="TADs", x = paste0(" TAD FCC"), y="") + 
  geom_density_ridges_gradient(scale = 0.95, size = 0.1, rel_min_height = 0.001) +
  eval(parse(text=paste0("scale_fill_", ggsci_pal, "(", ggsci_subpal, ")"))) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 8))+
  my_box_theme+
  theme(
    axis.line = element_line(),
    axis.text.x = element_text(size=14),
    axis.text.y = element_text(size=14),
    legend.text =element_text(size=12)
  )
outFile <- file.path(outFolder, paste0("FCC_signifAll_density_ggridge_v1.", plotType))
ggsave(ridge1_p, filename = outFile, height=myHeightGG, width=myWidthGG)
cat(paste0("... written: ", outFile, "\n"))


ridge2_p <- ggplot(plot_dt, aes(x = FCC, y = tad_label, fill = tad_label)) +
  ggtitle(plotTit, subtitle = subTit)+
  labs(fill="TADs", x = paste0("TAD FCC"), y="") + 
  geom_density_ridges_gradient(scale = 1.5, size = 0.1, rel_min_height = 0.001) +
  eval(parse(text=paste0("scale_fill_", ggsci_pal, "(", ggsci_subpal, ", alpha=0.7)"))) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 8))+
  my_box_theme+
  theme(
    axis.line = element_line(),
    axis.text.x = element_text(size=14),
    axis.text.y = element_text(size=14),
    legend.text =element_text(size=12)
  )


outFile <- file.path(outFolder, paste0("FCC_signifAll_density_ggridge_v2.", plotType))
ggsave(ridge2_p, filename = outFile, height=myHeightGG, width=myWidthGG)
cat(paste0("... written: ", outFile, "\n"))








