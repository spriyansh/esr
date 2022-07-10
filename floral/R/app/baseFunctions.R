setwd("/home/vega/James/esr/floral/R")
# styler::style_dir()

# Load the library
suppressPackageStartupMessages(library(flowCore))
library(flowViz)
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(ggcyto))
suppressPackageStartupMessages(library(flowStats))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggridges))

# Load the dataset
fcsObject <- read.FCS("exampleInputFile.fcs", transformation = FALSE)

# Exploration
summary(fcsObject)
summary(fcsObject) # Set up in table
colnames(fcsObject) # Drop down
slotNames(fcsObject) # Always fixed

exprsData <- as.data.frame(fcsObject@exprs)

exprsData <- exprsData[, c("Time", "FSC-A", "FSC-W")]



p <- ggplot(exprsData) + geom_hex(aes(x = `FSC-A`, y =  `FSC-W`),bins =200) + theme_bw()+
  theme(
    plot.title = element_text(size = 15),
    strip.background = element_blank(),
    strip.text = element_blank(),
    panel.background = element_blank(),
    panel.grid.major = element_line(size = rel(0.3), linetype = 2, colour = "#e65320"),
    panel.grid.minor = element_line(size = rel(0.1), linetype = 1, colour = "#e65320")
  )
p
lg <- flowStats::lymphGate(fcsObject, channels=c("FSC-A", "SSC-H"), scale=16)
fres <- filter(fcsObject, lg)
p <- p + geom_gate(lg)+ geom_stats()
p

library(hexbin)


hex <- hexbin(exprsData$`FSC-A`, exprsData$`FSC-W`, xbins = 200,
              colramp = colorRampPalette(hcl.colors(12, "GnBu")))
plot(hex)



xyplot("SSC-H" ~ "FSC-A" | Visit, data=fcsObject)
