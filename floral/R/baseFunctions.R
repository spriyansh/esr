setwd("/home/vega/James/esr/floral/R")
#styler::style_dir()

# Load the library
suppressPackageStartupMessages(library(flowCore))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(ggcyto))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(hrbrthemes))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggridges))

# Load the dataset
fcsObject <- read.FCS("../testData/exampleInputFile.fcs", transformation = FALSE)

# Exploration
summary(fcsObject)
summary(fcsObject) # Set up in table
colnames(fcsObject) # Drop down
slotNames(fcsObject) # Always fixed

# Check and Load
for (i in slotNames(fcsObject)) {
  if (!is.null(slot(fcsObject, name = i))) {
    if (i == "exprs") {
      exprsDataSlot <- as.data.frame(slot(fcsObject, name = i))
    } else if (i == "parameters") {
      parametersSlot <- slot(fcsObject, name = i)
      parametersDataSlot <- slot(parametersSlot, name = "data")
      parametersVarMetaDataSlot <- slot(parametersSlot, name = "varMetadata")
    } else if (i == "description") {
      descriptionDataSlot <- as.data.frame(slot(fcsObject, name = i))
    }
  } else {
    print(paste(i, "slot is empty"))
  }
}

# Melt Expression Data
exprsDataSlotMelted <- as.data.frame(melt(exprsDataSlot, id = "Time"))


ggplot(data = exprsDataSlotMelted,  aes(x=value, color=variable, fill=variable)) +
  geom_histogram(alpha=0.6, binwidth = 5) +
  scale_fill_viridis(discrete=TRUE) +
  scale_color_viridis(discrete=TRUE) +
  theme_ipsum() +
  theme(
    legend.position="none",
    panel.spacing = unit(0.1, "lines"),
    strip.text.x = element_text(size = 8)
  ) +
  xlab("") +
  ylab("Assigned Probability (%)") +
  facet_wrap(~variable)


exprs_FSC_A <- exprsDataSlotMelted[exprsDataSlotMelted$variable == "FSC-A",]

K = round(1 + 3.322*log(nrow(exprs_FSC_A)))
R = round(max(exprs_FSC_A$value)) - round(min(exprs_FSC_A$value))
B = round(R/K)

ggplot(data = exprs_FSC_A, aes(x=value)) +
  geom_histogram( binwidth=B, fill="#69b3a2", color="#e9ecef", alpha=0.9) + ggtitle("FSC_A")+
  theme_ipsum() +
  theme(
    plot.title = element_text(size=15)
  )



ggplot(data = exprs_FSC_A, aes(x=value)) +
  geom_density(fill="#69b3a2", color="#e9ecef", alpha=0.2) +
  theme_ipsum() +
  theme(
    plot.title = element_text(size=15)
  )
