# redirect R output to log
log <- file(snakemake@log[[1]], open="wt")
sink(log, type = "output")
sink(log, type = "message")

# load required libraries
library(tidyverse)
library(cowplot)
library(reshape2)

# plot data
plot <- read.csv(snakemake@input[[1]], header = TRUE) %>%
  dplyr::select(-c("total_reads")) %>%
  melt(id.vars = "sample",
       value.name = "Number of reads") %>%
  mutate("Number of reads" = `Number of reads` / 1000000) %>%
  ggplot(aes(x = sample, 
             y = `Number of reads`, 
             fill = variable)) +
  geom_bar(stat = "identity", 
           position = position_stack(reverse = TRUE),
           colour = "black") +
  theme_cowplot(16) +
  labs(x = NULL, y = bquote('Number of reads '~(x10^6))) +
  theme(legend.title = element_blank()) +
  scale_fill_manual(values = c("black","grey30","grey60","grey90", "forestgreen")) +
  scale_y_continuous(expand = c(0,0)) +
  guides(fill = guide_legend(reverse = TRUE)) +
  scale_x_discrete(guide = guide_axis(angle = 90))
  
#save plot
ggsave(plot, 
       filename = snakemake@output[[1]]
       )


# close redirection of output/messages
sink(log, type = "output")
sink(log, type = "message")
