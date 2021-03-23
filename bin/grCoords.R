#!/usr/bin/env Rscript

## grCoords: Group coordinates
## Locate chromatin diminution from long reads with trimmed 
## telomeric repeats mapped to the assembly
##
## This is the automated script for grouping coordinates 
## from a mapping coordinate file

library(dplyr)

getMode <- function(x){
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}


args <- commandArgs(TRUE)
inFile <- args[[1]]
outFile <- args[[2]]

# Read input
endReadSites <- read.delim(inFile,
                           header = F, sep = "\t",
                         col.names = c("query_name", "flag",
                                       "target_name", "target_start",
                                       "mapq",  "target_end")) 
# Filter by MAPQ and SAM flags
hqReadSites <- filter(endReadSites, !(bitwAnd(256, flag) | 
                                        bitwAnd(2048, flag)),
                      mapq > 10)

# Group coordinates
teloBlocks <- mutate(hqReadSites, strand = ifelse(bitwAnd(16, flag), "-", "+"),
                     rstart = ifelse(strand == "-",
                                     target_end, target_start)) %>%
  arrange(target_name, rstart) %>%
  group_by(target_name, strand) %>%
  mutate(block = cumsum(c(1, diff(rstart) > 6))) %>% 
  group_by(target_name, strand, block) %>%
  summarise(strand = unique(strand),
            teloStart = getMode(ifelse(strand == "+", target_start,
                                       target_end)),
            regSupport = n(),
            .groups = "drop") %>%
  group_by(target_name) %>%
  filter(regSupport > max(regSupport) * 0.05) %>%
  filter(regSupport > 1) %>%
  ungroup() %>%
  select(target_name, teloStart, strand, regSupport)

# Write output
write.table(teloBlocks, outFile, sep = "\t", row.names = F,
            quote = F)
