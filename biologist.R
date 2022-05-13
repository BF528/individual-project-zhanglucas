library(GSEABase)
library(GO.db)
library(tidyverse)
library(BiocManager)
library(biomaRt) #needed for ensembl
library(readr)
library(dplyr)
library(hgu133plus2.db)

path = "/project/bf528/project_1/data/differential_expression_results.csv"
rawData <- read.csv(file = path)
diff_exp <- setNames(cbind(rownames(rawData), rawData, row.names = NULL), 
         c("PID", "t", "p", "padj"))

#Add additional column to contain Symbol 
keys <- pull(diff_exp['PID'])
match <- select(hgu133plus2.db, keys, columns = c('PROBEID', 'SYMBOL'))
merged <- merge(match, diff_exp, by.x = "PROBEID", by.y = "PID")
merged <- merged %>%
  filter(!is.na(SYMBOL)) %>%
  group_by(SYMBOL) %>% 
  filter(padj == min(padj)) %>% 
  ungroup(SYMBOL)

#Top 10 up and down regulated genes 
upreg <- top_n(merged, 1000, t)
downreg <- top_n(merged, -1000, t)
tenup <- top_n(upreg, 10, t)
tendown <- top_n(downreg, -10, t)

#Get .GMT gene sets 
go <- getGmt('~/individual-project-zhanglucas/samples/go.gmt')
hall <- getGmt('~/individual-project-zhanglucas/samples/hall.gmt')
kegg <- getGmt('~/individual-project-zhanglucas/samples/kegg.gmt')

#10402 Elements in GO 
#50 Elements in Hall 
#186 Elements in KEGG

#Contigency Table Creation: 
cont_table <- function(gs, ngs, total) {
  ex_gs <- length(intersect(gs, ngs))
  not_ex_gs <- length(intersect(ngs, total))
  ex_notgs <- length(ngs) - ex_gs
  not_ex_notgs <- length(total) - not_ex_gs
  return(c(ex_gs, ex_notgs, not_ex_gs, not_ex_notgs))
}

up_no_exp <- subset(merged, !merged$SYMBOL %in% upreg$SYMBOL)
down_no_exp <- subset(merged, !merged$SYMBOL %in% downreg$SYMBOL)

#Function to take a gene set and produce the top results after fischer's test
stats_table <- function(set) {
  stats <- data.frame(gs_name = character(), stat_est = as.numeric(), pval = as.numeric())
  for (i in 1:length(set)){
    id <- geneIds(set[i])
    
    up_test <- fisher.test(matrix(cont_table(upreg $ SYMBOL, 
                                                  id[[names(id)]], up_no_exp$SYMBOL),nrow=2))
    down_test <- fisher.test(matrix(cont_table(downreg $ SYMBOL, 
                                                    hall_id[[names(id)]], down_no_exp$SYMBOL), nrow=2))
    stats[nrow(stats) + 1,] <- c(names(id), as.numeric(up_test$estimate), up_test$p.value)
    stats[nrow(stats) + 1,] <- c(names(id), as.numeric(down_test$estimate), down_test$p.value)
    
  }
  stats$adjp <- p.adjust(stats$pval, method = "BH",
                         n = length(stats$pval))
  
  sig <- stats%>%
    filter(adjp < 0.05) 
  
  print(count(sig))
  
  stats$statnum <- as.numeric(stats$stat_est)
  
  stats <- stats %>%
    filter(!(statnum == "Inf"))
  
  top <- top_n(stats, 3, statnum)
  return(top)
  
}

hallresult <- stats_table(hall)
keggresult <- stats_table(kegg)
goresult <- stats_table(go)



