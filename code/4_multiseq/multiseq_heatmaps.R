# this script generates heatmaps in Figures S14 and S15

library(dplyr)
library(DT)
library(purrr)
library(ggplot2)
library(ggpubr)
library(scales)
library(tibble)

# data are not included in the repo
df_sieve <- read.csv(here::here("data/HVTN_705_multiseq_sieve_data_v3.csv"))
df_sieve_primary_ids <- read.csv(here::here("data/HVTN705_subjid_prim_endpts_with_seq_sieve.csv"))

# filter for first visit only 
df_sieve_visit_to_keep <- df_sieve %>% 
  arrange(subjid, visit) %>% 
  select(subjid, visit) %>% 
  distinct(subjid, .keep_all = T)

df_sieve <- df_sieve %>% merge(df_sieve_visit_to_keep, by = c("subjid", "visit"))

df_sieve <- df_sieve %>% 
  filter(subjid %in% df_sieve_primary_ids$subjid)

##############################################################################
# Heatmap 1 
##############################################################################

screened_in_residues <- c(".51.is.P", ".51.is.T", ".57.is.D", ".57.is.N", 
                          ".59.is.K", ".60.is.A", ".60.is.G", ".60.is.S", 
                          ".61.is.H", ".72.is.H", ".78.is.D", ".130.is.E", 
                          ".130.is.H", ".130.is.K", ".130.is.N", ".130.is.S", 
                          ".130.is.T", ".155.is.K", ".155.is.Q", ".155.is.R", 
                          ".155.is.T", ".158.is.S", ".160.is.K", ".160.is.N", 
                          ".161.is.A", ".161.is.I", ".161.is.M", ".161.is.T", 
                          ".161.is.V", ".162.is.T", ".163.is.S", ".163.is.T", 
                          ".164.is.E", ".165.is.I", ".165.is.L", ".165.is.V", 
                          ".166.is.K", ".166.is.R", ".166.is.T", ".167.is.D", 
                          ".167.is.N", ".168.is.K", ".168.is.R", ".169.is.E", 
                          ".169.is.I", ".169.is.K", ".169.is.Q", ".169.is.R", 
                          ".169.is.T", ".170.is.K", ".170.is.Q", ".170.is.R", 
                          ".171.is.K", ".171.is.N", ".171.is.Q", ".171.is.R", 
                          ".172.is.A", ".172.is.E", ".172.is.M", ".172.is.T", 
                          ".172.is.V", ".173.is.H", ".173.is.R", ".173.is.S", 
                          ".173.is.Y", ".175.is.L", ".177.is.Y", ".178.is.K", 
                          ".178.is.R", ".179.is.L", ".179.is.P", ".179.is.S", 
                          ".179.is.V", ".181.is.I", ".181.is.L", ".181.is.T", 
                          ".181.is.V", ".182.is.E", ".182.is.I", ".182.is.T", 
                          ".182.is.V", ".183.is.A", ".183.is.K", ".183.is.P", 
                          ".183.is.Q", ".183.is.S", ".184.is.I", ".184.is.L", 
                          ".186.is.D", ".186.is.E", ".186.is.G", ".186.is.K", 
                          ".186.is.N", ".186.is.S", ".186.is.T", ".186.is.gap", 
                          ".200.is.A", ".200.is.T", ".200.is.V", ".230.is.D", 
                          ".230.is.N", ".234.is.N", ".234.is.S", ".279.is.D", 
                          ".279.is.N", ".280.is.N", ".280.is.S", ".281.is.A", 
                          ".281.is.I", ".281.is.T", ".281.is.V", ".282.is.K", 
                          ".350.is.A", ".350.is.E", ".350.is.G", ".350.is.K", 
                          ".350.is.Q", ".350.is.R", ".364.is.A", ".364.is.H", 
                          ".364.is.P", ".364.is.S", ".429.is.E", ".429.is.G", 
                          ".429.is.K", ".429.is.R", ".432.is.K", ".432.is.Q", 
                          ".432.is.R", ".456.is.R", ".456.is.W", ".471.is.A", 
                          ".471.is.E", ".471.is.G", ".471.is.I", ".471.is.Q", 
                          ".471.is.T", ".471.is.V")

sequon_vars <- c("130.is.sequon", "156.is.sequon", 
                 "160.is.sequon", "234.is.sequon", 
                 "332.is.sequon", "392.is.sequon")

heatmap1_pattern <- paste(gsub("\\.", "\\\\.", 
                               c(screened_in_residues, sequon_vars)), collapse = "|")

cols_to_summarize <- grep(heatmap1_pattern, colnames(df_sieve), value = TRUE)

heatmap1_summary <- df_sieve %>%
  group_by(subjid) %>% 
  summarise(across(all_of(cols_to_summarize), ~ mean(.x, na.rm = TRUE)))

heatmap1_summary$proportion_not_0_or_1 <- 
  apply(heatmap1_summary[ , -which(names(heatmap1_summary) == "subjid")], 1, 
        function(row) {mean(row != 0 & row != 1)})

# subjid order
subjid_order <- (heatmap1_summary %>% arrange(proportion_not_0_or_1))$subjid

# column order 
proportion_not_0_or_1 <- sapply(heatmap1_summary %>% select(-subjid, -proportion_not_0_or_1), 
                                function(column) {
                                  mean(column != 0 & column != 1)
                                })

# change names of columns
make_columns_readable <- function(vec){
  vec <- gsub("^[^.]+\\.(\\d+)\\.is\\.([^.]+)\\..+$","\\2\\1",  vec)
  vec <- gsub("^[^.]+\\.(\\d+)\\.is\\.(.+)$",  "\\2\\1", vec)
  vec <- gsub("gap",  "Gap", vec)
  vec <- gsub("sequon",  "PNG", vec)
}

ordered_columns <- names(sort(proportion_not_0_or_1))
ordered_columns <- make_columns_readable(ordered_columns)

heatmap1_to_graph <- heatmap1_summary %>% 
  select(-proportion_not_0_or_1) %>%
  mutate(subjid = factor(subjid, levels = subjid_order)) %>% 
  reshape2::melt(id = "subjid")  %>% 
  mutate(value = as.numeric(value)) %>% 
  mutate(variable = make_columns_readable(variable),
         variable = factor(variable, levels = ordered_columns))

plt1 <- ggplot(aes(y = factor(subjid), x = variable, fill = value), 
               data = heatmap1_to_graph) + 
  geom_tile() +
  scale_fill_gradientn(
    colors = c("white", "#ebd2fc", "purple", "#283bc2", "#041483"),
    values = scales::rescale(c(0, 0.001, 0.5, 0.999, 1)),
    limits = c(0, 1)
  ) +
  ylab("")+
  xlab("") + 
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
    axis.text.y = element_blank(),
    axis.title.x = element_blank(),
    legend.text = element_text(size = 9),   
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank()
    
  )  +
  labs(fill = "") +
  theme(legend.position = "bottom")

ggsave("heatmap1.pdf", plot = plt1, height = 8, width = 15, device = "pdf")

# some numbers for paper
# length(proportion_not_0_or_1[proportion_not_0_or_1 == 0])
# length(proportion_not_0_or_1[proportion_not_0_or_1 <= 0.05])
# length(proportion_not_0_or_1[proportion_not_0_or_1 <= 0.2])
# 
# heatmap1_summary_not_dominant <- heatmap1_summary %>%
#   mutate(across(
#     where(is.numeric) & !all_of("subjid"),
#     ~ if_else(. > 0.5, 1 - ., .)
#   ))
# 
# heatmap1_summary_not_dominant_max <- heatmap1_summary_not_dominant %>%
#   select(-subjid) %>%
#   summarise(across(everything(), max, na.rm = TRUE)) %>% 
#   reshape2::melt() %>% 
#   arrange(value)
# 
# heatmap1_summary_not_dominant_max


##############################################################################
# Heatmap 2 
##############################################################################

# change this line to modify the number of standard deviations around the mode
num_sd <- 0.5

mode <- function(x, na.rm = TRUE) {
  if (na.rm) {
    x <- na.omit(x)
  }
  unique_x <- unique(x)
  unique_x[which.max(tabulate(match(x, unique_x)))]
}

continuous_marks <- c("hdist.zspace.c97za.v2", "hdist.zspace.c97za.v2_ab", 
                      "hdist.zspace.c97za.adcc", "hdist.zspace.c97za.c_ab", 
                      "hdist.zspace.c97za.hvtn505.cd4bs.kmer",
                      "hdist.zspace.c97za.hvtn505.cd4bs.antibody", 
                      "hdist.zspace.1428v1v2", "length.v1v2", "num.sequons.v1v2", 
                      "charge.v2", "cysteine.count")

mode_list <- df_sieve %>%
  group_by(subjid) %>% 
  summarise(across(all_of(continuous_marks), ~ mode(.x, na.rm = TRUE))) %>%
  rowwise() %>%
  mutate(
    nested_list = list(setNames(c_across(all_of(continuous_marks)), continuous_marks))
  ) %>%
  select(subjid, nested_list) %>%
  deframe()

mean_sd_list <- df_sieve %>%
  group_by(subjid) %>% 
  summarise(across(all_of(continuous_marks), ~ sd(.x, na.rm = TRUE))) %>% 
  summarise(across(all_of(continuous_marks), mean, na.rm = TRUE)) %>%
  as.list()

heatmap2_summary <- df_sieve %>%
  group_by(subjid) %>% 
  summarise(across(all_of(continuous_marks), 
                   ~ mean((.x < flatten(mode_list[as.character(subjid)])[[cur_column()]] - 
                             (num_sd * mean_sd_list[[cur_column()]])) | 
                            (.x > flatten(mode_list[as.character(subjid)])[[cur_column()]] +
                               (num_sd * mean_sd_list[[cur_column()]]))))) 

heatmap2_summary$avg_proportion <- apply(heatmap2_summary[ , -which(names(heatmap2_summary) == "subjid")], 
                                         1, function(row) {mean(row)})

avg_proportion <- sapply(heatmap2_summary %>% select(-subjid, -avg_proportion), 
                         function(column) {
                           mean(column)
                         })

subjid_order <- (heatmap2_summary %>% arrange(avg_proportion))$subjid
ordered_columns <- names(sort(avg_proportion))

heatmap2_to_graph <- heatmap2_summary %>%
  select(-avg_proportion) %>%
  mutate(subjid = factor(subjid, levels = subjid_order)) %>% 
  reshape2::melt(id = "subjid")  %>% 
  mutate(value = as.numeric(value)) %>% 
  mutate(variable = factor(variable, levels = ordered_columns))


plt2 <- ggplot(aes(y = factor(subjid), x = variable, fill = value), 
               data = heatmap2_to_graph) + 
  geom_tile() +
  scale_fill_gradient2(low = "white", high = "blue") + 
  ylab("")+
  theme(
    axis.text.x = element_text(size = 7.5),
    axis.text.y = element_blank(),
    legend.text = element_text(size = 13),   
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.title.x = element_blank()
  )  +
  labs(fill = "") + 
  xlab("") +
  scale_x_discrete(labels = c(
    "charge.v2" = "V2 charge", 
    "length.v1v2" = "V1V2 length", 
    "num.sequons.v1v2" = "Number of PNG sites\nin V1V2", 
    "cysteine.count" = "Cysteine count\nin Env gp160", 
    "hdist.zspace.c97za.v2" = "PC-weighted Hamming\ndistance to C97ZA in\nV2 hotspot in NHP\nchallenge model", 
    "hdist.zspace.c97za.v2_ab" = "PC-weighted Hamming\ndistance to C97ZA in\nIgG3 V1V2 correlates set", 
    "hdist.zspace.c97za.adcc" = "PC-weighted Hamming\ndistance to C97ZA in\nC1 epitope of ADCC\ncorrelate in RV144", 
    "hdist.zspace.c97za.c_ab" = "PC-weighted Hamming\ndistance to C97ZA in\nclade C bnAb signatures", 
    "hdist.zspace.c97za.hvtn505.cd4bs.kmer" = "PC-weighted Hamming\ndistance to C97ZA in\nCD4bs-overlapping\nsignature k-mers in\nHVTN 505",
    "hdist.zspace.c97za.hvtn505.cd4bs.antibody" = "PC-weighted Hamming\ndistance to C97ZA in\nsignature CD4bs\nantibody footprints\nin HVTN 505", 
    "hdist.zspace.1428v1v2" = "PC-weighted Hamming\ndistance to gp70-001428.2.42\nV1V2 in IgG3 V1V2\ncorrelates set"
  )) +
  theme(legend.position = "bottom")

ggsave(here::here("figures/heatmap2.pdf"), plot = plt2, height = 10, width = 15, device = "pdf")

