###
# pub_figures.r
# The purpose of this script is to generate all figures for the writeup
# It uses post-processed data in results so that it can be run without
# signalP and phobius related dependencies.
###

library(here)
library(tidyverse)
library(cowplot)
library(ggExtra)
library(ggplotify)
library(ggthemes)
library(vcd)
library(Biostrings)
source(here("src", "utils.r"))
theme_set(theme_cowplot(font_size = 10) +
  theme(
    strip.background = element_blank(),
    plot.margin = unit(c(0, 0, 0, 0), units = "inches")
  ))


# Load data with S. cerevisiae hydropathy calculation

hydropathy_df <- read_csv(here("results", "figures", "SC_first_60.csv"))

screened_non_srp <- read_lines(here("data", "SC_screened.txt"))
verified_srp <- read_lines(here("data", "SC_SRP.txt"))
verified_non_srp <- read_lines(here("data", "SC_non_SRP.txt"))

labelled_df <- hydropathy_df %>%
  mutate(
    `Experimental label` =
      case_when(
        seqid %in% verified_non_srp ~ "Sec63-dependent",
        seqid %in% screened_non_srp ~ "Sec63-dependent",
        seqid %in% verified_srp ~ "SRP-dependent",
        TRUE ~ "Unverified"
      )
  )

verified_df <- labelled_df %>%
  filter(`Experimental label` != "Unverified")


# plot Kyte-Doolittle hydrophobicity axis
rough_KD_limits <- c(
  min(labelled_df$KD_max_hydropathy),
  max(labelled_df$KD_max_hydropathy)
)
rough_KD_limits

scale_y_KD_hydropathy <-
  scale_y_continuous("Max. hydropathy (Kyte-Doolittle)",
    breaks = 0:4,
    limits = c(0, 4.5),
    expand = c(0, 0)
  )

# Figure 0

# Plot of contingency table of helix length predicted by phobius
# and those verified experimentally

# make contingency table for display
contingency_df_length <- verified_df %>%
  mutate(
    `Experimental label` =
      factor(`Experimental label`,
        levels = c("Sec63-dependent", "SRP-dependent")
      )
  ) %>%
  group_by(`Experimental label`) %>%
  summarise(
    short = sum(window_length < 14),
    long = sum(window_length >= 14)
  )

contingency_table_length <- as.table(as.matrix(contingency_df_length[, 2:3]))

# make table pretty for display
names(dimnames(contingency_table_length)) <- c("Experimental label", "Helix length")
rownames(contingency_table_length) <- c("Sec63", "SRP")
# For a two-way table, mosaic() fits a model of independence, [A][B] or ~A+B as an R formula
# https://www.datavis.ca/courses/VCD/vcd-tutorial.pdf
ScHydropathy_contingency_plot_length <- as.grob(~ vcd::mosaic(contingency_table_length, shade = TRUE, legend = TRUE, main = "Verified SP/TM regions vs Phobius predictions"))


# Contingency table of SP/TM regions predicted by phobius
# and those verified experimentally

# run chi-squared independence test and extract p-value
# the null hypothesis is that the two categorical variables are independent
# the p-value rejects this, so there is a high association between the phobius label and the experimental label
test_length <- chisq.test(contingency_table_length)
p_value_length <- test_length$p.value
p_value_length

# make contingency table for display
contingency_df_label <- verified_df %>%
  mutate(
    `Experimental label` =
      factor(`Experimental label`,
        levels = c("Sec63-dependent", "SRP-dependent")
      )
  ) %>%
  group_by(`Experimental label`) %>%
  summarise(
    SP = sum(window_type == "SP"),
    TM = sum(window_type == "TM")
  )

contingency_table_label <- as.table(as.matrix(contingency_df_label[, 2:3]))

# make table pretty for display
names(dimnames(contingency_table_label)) <- c("Experimental label", "Phobius label")
rownames(contingency_table_label) <- c("Sec63", "SRP")

# run chi-squared independence test and extract p-value
# the null hypothesis is that the two categorical variables are independent
# the p-value rejects this, so there is a high association between the phobius label and the experimental label
test_label <- chisq.test(contingency_table_label)
p_value_label <- test_label$p.value



# For a two-way table, mosaic() fits a model of independence, [A][B] or ~A+B as an R formula
# https://www.datavis.ca/courses/VCD/vcd-tutorial.pdf
ScHydropathy_contingency_plot_label <- as.grob(~ vcd::mosaic(contingency_table_label, shade = TRUE, legend = TRUE, main = "Verified SP/TM regions vs Phobius predictions"))



# Plot of contingency table of compound hydropathy
# and those verified experimentally

# make contingency table for display
contingency_df_compound <- verified_df %>%
  mutate(
    `Experimental label` =
      factor(`Experimental label`,
        levels = c("Sec63-dependent", "SRP-dependent")
      )
  ) %>%
  mutate(compound = window_length * KD_max_hydropathy) %>%
  group_by(`Experimental label`) %>%
  summarise(
    hi = sum(compound < 40),
    lo = sum(compound >= 40)
  )

contingency_table_compound <- as.table(as.matrix(contingency_df_compound[, 2:3]))

# make table pretty for display
names(dimnames(contingency_table_compound)) <- c("Experimental label", "Compound hydropathy")
rownames(contingency_table_compound) <- c("Sec63", "SRP")

# For a two-way table, mosaic() fits a model of independence, [A][B] or ~A+B as an R formula
# https://www.datavis.ca/courses/VCD/vcd-tutorial.pdf
ScHydropathy_contingency_plot_compound <- as.grob(~ vcd::mosaic(contingency_table_compound, shade = TRUE, legend = TRUE, main = "Verified SP/TM regions vs Phobius predictions"))

# run chi-squared independence test and extract p-value
# the null hypothesis is that the two categorical variables are independent
# the p-value rejects this, so there is a high association between the phobius label and the experimental label
test_compound <- chisq.test(contingency_table_compound)
p_value_compound <- test_compound$p.value
p_value_compound



# Scatter plot of all proteins with SP/TM regions found by phobius
# Highlighting verified proteins

breaks_explabel <- c(
  "Sec63-dependent",
  "SRP-dependent",
  "Unverified"
)
colour_explabel <- c(
  "Sec63-dependent" = "blue",
  "SRP-dependent" = "red",
  "Unverified" = "plum4"
)
size_explabel <- c(
  "Sec63-dependent" = 1.5,
  "SRP-dependent" = 1.5,
  "Unverified" = 0.5
)

compound_hydropathy_40linedf <-
  tibble(
    window_length = seq(4, 34, 0.2),
    KD_max_hydropathy = 40 / window_length
  )


# Figure 1 - Scatter plot of Kyte-Doolittle hydropathy vs helix length

base_ScScatMarg <-
  ggplot(labelled_df, aes(
    x = window_length,
    y = KD_max_hydropathy
  )) +
  geom_point(aes(
    colour = `Experimental label`,
    group = `Experimental label`,
    size = `Experimental label`
  )) +
  geom_vline(xintercept = 13.5, linetype = "dashed") +
  geom_line(
    data = compound_hydropathy_40linedf,
    # linetype = "dotted"
    colour = "grey50"
  ) +
  annotate(
    label = "compound\nhydropathy\n= 40", geom = "text",
    x = 34, y = 1.15, hjust = 1, vjust = 1,
    colour = "grey50"
  ) +
  # ggtitle("Phobius detected SP/TM regions") +
  scale_x_helix_length +
  scale_y_KD_hydropathy +
  scale_colour_manual(breaks = breaks_explabel, values = colour_explabel) +
  scale_size_manual(breaks = breaks_explabel, values = size_explabel) +
  theme(
    legend.box.background = element_rect(colour = "grey60"),
    legend.box.margin = margin(2, 2, 2, 2, unit = "pt"),
    legend.justification = c(0.5, 0.5),
    plot.margin = unit(c(.1, .1, .1, .1), "mm")
  )

side_ScScatMarg <-
  ggplot(
    labelled_df %>%
      dplyr::mutate(
        `Experimental label` =
          factor(`Experimental label`,
            levels = rev(breaks_explabel)
          )
      ),
    aes(
      y = KD_max_hydropathy,
      fill = `Experimental label`,
      group = `Experimental label`
    )
  ) +
  geom_histogram(binwidth = 0.2) +
  scale_y_KD_hydropathy +
  scale_fill_manual(breaks = breaks_explabel, values = colour_explabel) +
  facet_grid(. ~ `Experimental label`, scales = "free_x") +
  labs(x = "Number of proteins") +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_blank(),
    # axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    plot.margin = unit(c(.1, .1, .1, .1), "mm")
  )

top_ScScatMarg <- ggplot(
  labelled_df %>%
    dplyr::mutate(
      `Experimental label` =
        factor(`Experimental label`,
          levels = breaks_explabel
        )
    ),
  aes(
    x = window_length,
    fill = `Experimental label`,
    group = `Experimental label`
  )
) +
  geom_histogram(binwidth = 1) +
  geom_vline(xintercept = 13.5, linetype = "dashed") +
  scale_x_helix_length +
  scale_fill_manual(values = colour_explabel) +
  facet_grid(`Experimental label` ~ ., scales = "free_y") +
  labs(y = "Number of proteins") +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_blank(),
    # axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    plot.margin = unit(c(.1, .1, .1, .1), "mm")
  )

ScHydropathy_scatter_marginals_plot <-
  plot_grid(top_ScScatMarg,
    get_legend(base_ScScatMarg),
    base_ScScatMarg +
      theme(legend.position = "none"),
    side_ScScatMarg,
    ncol = 2,
    align = "hv",
    axis = "bl",
    rel_heights = c(0.6, 1),
    rel_widths = c(1, 0.75)
  )
ScHydropathy_scatter_marginals_plot

# save plot
ggsave(
  filename = here("results", "figures", "ScHydropathy_scatter_marginals_KD.pdf"),
  plot = ScHydropathy_scatter_marginals_plot,
  width = 6.5, height = 5.5, dpi = 300
)

# save plot as .png for google docs input
ggsave(
  filename = here("results", "figures", "ScHydropathy_scatter_marginals_KD.png"),
  plot = ScHydropathy_scatter_marginals_plot,
  width = 6.5, height = 5.5, dpi = 300, bg = "white"
)

## Same figure but for Hessa hydrophobicity
base_ScScatHessa <-
  ggplot(labelled_df, aes(
    x = window_length,
    y = Hessa_max_hydropathy
  )) +
  geom_point(aes(
    colour = `Experimental label`,
    group = `Experimental label`,
    size = `Experimental label`
  )) +
  geom_vline(xintercept = 13.5, linetype = "dashed") +
  # ggtitle("Phobius detected SP/TM regions") +
  scale_x_helix_length +
  # scale_y_KD_hydropathy +
  scale_y_continuous("Hessa max. hydropathy") +
  scale_colour_manual(breaks = breaks_explabel, values = colour_explabel) +
  scale_size_manual(breaks = breaks_explabel, values = size_explabel) +
  theme(
    legend.box.background = element_rect(colour = "grey60"),
    legend.box.margin = margin(2, 2, 2, 2, unit = "pt"),
    legend.justification = c(0.5, 0.5),
    plot.margin = unit(c(.1, .1, .1, .1), "mm")
  )

side_ScScatHessa <-
  ggplot(
    labelled_df %>%
      dplyr::mutate(
        `Experimental label` =
          factor(`Experimental label`,
            levels = rev(breaks_explabel)
          )
      ),
    aes(
      y = Hessa_max_hydropathy,
      fill = `Experimental label`,
      group = `Experimental label`
    )
  ) +
  geom_histogram(binwidth = 0.1) +
  # scale_y_KD_hydropathy +
  scale_fill_manual(breaks = breaks_explabel, values = colour_explabel) +
  facet_grid(. ~ `Experimental label`, scales = "free_x") +
  labs(x = "Number of proteins") +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_blank(),
    # axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    plot.margin = unit(c(.1, .1, .1, .1), "mm")
  )


ScHessa_scatter_marginals_plot <-
  plot_grid(top_ScScatMarg,
    get_legend(base_ScScatMarg),
    base_ScScatHessa +
      theme(legend.position = "none"),
    side_ScScatHessa,
    ncol = 2,
    align = "hv",
    axis = "bl",
    rel_heights = c(0.6, 1),
    rel_widths = c(1, 0.75)
  )
ScHessa_scatter_marginals_plot

# save plot
ggsave(
  filename = here("results", "figures", "ScHydropathy_scatter_marginals_Hessa.pdf"),
  plot = ScHessa_scatter_marginals_plot,
  width = 6.5, height = 5.5, dpi = 300
)

ggsave(
  filename = here("results", "figures", "ScHydropathy_scatter_marginals_Hessa.png"),
  plot = ScHessa_scatter_marginals_plot,
  width = 6.5, height = 5.5, dpi = 300, bg = "white"
)


# Figure 4 - histograms of window lengths for each species
# --- along with GO term labels for each species --- #

phobius_composite_plot_fungi <- phobius_cladogram_plot("proteome_table_fungi.txt", "((((((((Sc:1,Ca:1):1,((Nc:1,Mg:1):1,(Zt:1,Af:1):1):1),Sp:1),(Pg:1,(Um:1,Cn:1):1):1):1):1,Rd:1):1,Bd:1):1);")
phobius_composite_plot_fungihuman <- phobius_cladogram_plot("proteome_table_fungihuman.txt", "(((((((((Sc:1,Ca:1):1,((Nc:1,Mg:1):1,(Zt:1,Af:1):1):1),Sp:1),(Pg:1,(Um:1,Cn:1):1):1):1):1,Rd:1):1,Bd:1):1,Hs:1):1);")

# save
ggsave(
  filename = here("results", "figures", "phobius_helix_length_fungi.pdf"),
  plot = phobius_composite_plot_fungi,
  width = 6, height = 8
)

ggsave(
  filename = here("results", "figures", "phobius_helix_length_fungi.png"),
  plot = phobius_composite_plot_fungi,
  width = 6, height = 8, bg = "white"
)

ggsave(
  filename = here("results", "figures", "phobius_helix_length_fungihuman.pdf"),
  plot = phobius_composite_plot_fungihuman,
  width = 6, height = 8
)

ggsave(
  filename = here("results", "figures", "phobius_helix_length_fungihuman.png"),
  plot = phobius_composite_plot_fungihuman,
  width = 6, height = 8, bg = "white"
)


# Figure 3 - DeepTMHMM and Phobius comparison scatter

# read string of here/results/proteins/SC_deeptmhmm/predicted_topologies.3line
deeptmhmm_3line <-
  here("results", "deepTMHMM", "S_Cerevisiae", "predicted_topologies.3line") %>%
  read_file() %>%
  # remove irrelevant transcript id
  str_remove_all("-t26_1") %>%
  strsplit(split = ">") %>%
  .[[1]]

deeptmhmm_df <- data.frame(
  seqid = character(),
  DeepTMHMM_type = character(),
  DeepTMHMM_length = numeric()
)

# extract SP and TM regions, within first 60 AA
for (i in 2:length(deeptmhmm_3line)) {
  # split 3line into array of strings
  protein_str <- deeptmhmm_3line[i]
  protein_arr <- str_split(protein_str, "\n")[[1]]

  # extract protein details
  details <- str_split(protein_arr[1], " ")[[1]]
  seqid <- details[1]
  DeepTMHMM_type <- details[3]

  # skip proteins that are not TM or SP
  if (DeepTMHMM_type == "GLOB") {
    next
  }

  # split topology string into groups of same characters
  topology_string <- protein_arr[3]
  topology_groups <- strsplit(topology_string, "(?<=(.))(?!\\1)", perl = TRUE)[[1]]

  # find end positions and filter for those <= 60 AA
  end_positions <- cumsum(lapply(topology_groups, str_length))
  topology_groups <- topology_groups[end_positions <= 60]

  # populate dataframes with regions
  SP_seqid_list <- c()
  for (j in seq_along(topology_groups)) {
    if (substr(topology_groups[j], 1, 1) == "S") {
      new_row <- data.frame(seqid = seqid, DeepTMHMM_type = "SP", DeepTMHMM_length = str_length(topology_groups[j]))
      deeptmhmm_df <- rbind(deeptmhmm_df, new_row)
      break
    } else if (substr(topology_groups[j], 1, 1) == "M") {
      new_row <- data.frame(seqid = seqid, DeepTMHMM_type = "TM", DeepTMHMM_length = str_length(topology_groups[j]))
      deeptmhmm_df <- rbind(deeptmhmm_df, new_row)
      break
    }
  }
}


# --- Comparing to phobius ---#

# read full length phobius results
full_phobius <- read_csv(here("results", "phobius", "S_Cerevisiae_fullSignal.csv")) %>%
  filter(phobius_end != 0) %>%
  mutate(phobius_length = phobius_end - phobius_start + 1) %>%
  select(seqid, Phobius_type = phobius_type, Phobius_length = phobius_length) %>%
  mutate(method = "Phobius")

labelled_phobius <- full_phobius %>%
  mutate(
    `Experimental label` =
      case_when(
        seqid %in% verified_non_srp ~ "Sec63-dependent",
        seqid %in% screened_non_srp ~ "Sec63-dependent",
        seqid %in% verified_srp ~ "SRP-dependent",
        TRUE ~ "Unverified"
      )
  )
verified_phobius <- labelled_phobius %>%
  filter(`Experimental label` != "Unverified")

# combined table of labelled proteins by classification method
combined_labelled <- full_join(labelled_phobius, deeptmhmm_df, by = "seqid")

# Compare DeepTMHMM predictions with verified proteins
# make contingency table for display
contingency_df_deeplabel <-
  combined_labelled %>%
  filter(`Experimental label` != "Unverified") %>%
  mutate(
    `Experimental label` =
      factor(`Experimental label`,
        levels = c("Sec63-dependent", "SRP-dependent")
      )
  ) %>%
  group_by(`Experimental label`) %>%
  summarise(
    SP = sum(DeepTMHMM_type == "SP", na.rm = TRUE),
    TM = sum(DeepTMHMM_type == "TM", na.rm = TRUE)
  )

contingency_table_deeplabel <- as.table(as.matrix(contingency_df_deeplabel[, 2:3]))

# make table pretty for display
names(dimnames(contingency_table_deeplabel)) <- c("Experimental label", "DeepTMHMM label")
rownames(contingency_table_deeplabel) <- c("Sec63", "SRP")

# run chi-squared independence test and extract p-value
# the null hypothesis is that the two categorical variables are independent
# the p-value rejects this, so there is a high association between the phobius label and the experimental label
test_deeplabel <- chisq.test(contingency_table_deeplabel)
p_value_deeplabel <- test_deeplabel$p.value
p_value_deeplabel

knitr::kable(contingency_table_deeplabel, caption = "Contingency table of SP/TM regions predicted by DeepTMHMM and verified translocon", format = "simple")



# contingency table by method
contingency_df_deepphob <-
  combined_labelled %>%
  group_by(DeepTMHMM_type, Phobius_type) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  pivot_wider(names_from = Phobius_type, values_from = count, values_fill = 0)

contingency_table_deepphob <- as.table(as.matrix(contingency_df_deepphob[1:2, 2:3]))

# make table pretty for display
names(dimnames(contingency_table_deepphob)) <- c("DeepTMHMM", "Phobius")
rownames(contingency_table_deepphob) <- c("SP", "TM")

# test if Phobius and DeepTMHMM agree, where they both make predictions
# run chi-squared independence test and extract p-value
# the null hypothesis is that the two categorical variables are independent
# the p-value rejects this, so there is a high association between the phobius label and the experimental label
test_deepphob <- chisq.test(contingency_table_deepphob)
p_value_deepphob <- test_length$p.value
p_value_deepphob

knitr::kable(contingency_table_deepphob, caption = "Contingency table of SP/TM regions predicted by both DeepTMHMM and Phobius", format = "simple")

# scatter of predicted lengths by method, coloured by label match

deepphob_lengthcor_df <-
  combined_labelled %>%
  dplyr::rename(Phobius = Phobius_type, DeepTMHMM = DeepTMHMM_type) %>%
  group_by(Phobius, DeepTMHMM) %>%
  summarise(
    cor_length = cor(Phobius_length, DeepTMHMM_length,
      use = "pairwise.complete.obs"
    ),
    count = n(),
    .groups = "drop"
  ) %>%
  mutate(
    cor_label = paste0("R = ", round(cor_length, digits = 2)),
    count_label = paste0("n = ", count),
    both_label = paste0(count_label, "\n", cor_label)
  )

deepphob_match_df <-
  combined_labelled %>%
  dplyr::rename(Phobius = Phobius_type, DeepTMHMM = DeepTMHMM_type) %>%
  group_by(Phobius, Phobius_length, DeepTMHMM, DeepTMHMM_length) %>%
  summarise(count = n(), .groups = "drop")

deepphob_match_plot <-
  ggplot(
    data = deepphob_match_df %>%
      drop_na() %>%
      dplyr::mutate(DeepTMHMM = factor(DeepTMHMM,
        levels = c("TM", "SP")
      )),
    aes(x = Phobius_length, y = DeepTMHMM_length)
  ) +
  geom_abline(slope = 1, intercept = 0, colour = "grey60") +
  geom_point(aes(size = count), colour = "darkgreen", alpha = 0.5) +
  geom_text(
    data = deepphob_lengthcor_df %>%
      drop_na() %>%
      dplyr::mutate(DeepTMHMM = factor(DeepTMHMM,
        levels = c("TM", "SP")
      )),
    aes(label = both_label),
    x = 40, y = 60, hjust = 0, vjust = 1, size = 3,
    inherit.aes = FALSE
  ) +
  facet_grid(DeepTMHMM ~ Phobius, labeller = label_both) +
  theme(panel.border = element_rect(fill = NA, colour = "grey90")) +
  tune::coord_obs_pred() +
  labs(
    x = "Phobius predicted length",
    y = "DeepTMHMM predicted length"
  ) +
  scale_x_continuous(breaks = seq(0, 100, 10)) +
  scale_y_continuous(breaks = seq(0, 100, 10))

deepphob_match_plot

# save plot
ggsave(
  filename = here("results", "figures", "Phobius_DeepTMHMM_length_match.pdf"),
  plot = deepphob_match_plot,
  width = 6.5, height = 5.5, dpi = 300
)

ggsave(
  filename = here("results", "figures", "Phobius_DeepTMHMM_length_match.png"),
  plot = deepphob_match_plot,
  width = 6.5, height = 5.5, dpi = 300, bg = "white"
)

# Figure 2 - AA Composition

# --- Amino acid composition analysis --- #

# attach path to protein file names
species_df <- here("data", "proteins", "pub", "proteome_table_fungi.txt") %>%
  read_tsv(comment = "#") %>%
  # Next line makes Nicename a factor in same order as given
  mutate(
    Nicename = as_factor(Nicename),
    Nicename_splitline =
      factor(Nicename,
        levels = Nicename,
        labels = str_replace(Nicename,
          pattern = " ",
          replacement = "\n"
        )
      )
  )
protein_paths <- here("data", "proteins", "pub", species_df$Filename)

# read in protein sequences
proteins <- lapply(protein_paths, readAAStringSet)

# get all S. cerevisiae proteins and subset by phobius prediction
SC_AA <- proteins[[1]]
SC_phobius_df <- labelled_df %>%
  select(seqid, window_length, window_start, window_end, window_type, `Experimental label`) %>%
  mutate(
    window_AA_string =
      mapply(
        function(x, y, z) toString(as.character(SC_AA[[x]][y:z])),
        seqid, window_start, window_end
      )
  )

# load KD hydrophocities and sort Amino Acids by KD hydrophobicity
sorted_AA <-
  read.csv(here("data", "scales.csv"), header = TRUE) %>%
  arrange(Kyte.Doolittle) %>%
  pull(aa)

count_AAs <- function(AA_list, group_name = NA, as.prob = FALSE) {
  output <- AA_list %>%
    AAStringSet() %>%
    letterFrequency(
      letters = AA_ALPHABET[1:20],
      as.prob = as.prob
    ) %>%
    as_tibble() %>%
    summarise_all(sum)

  if (!is.na(group_name)) {
    output <- output %>% mutate(group = group_name)
  }

  return(output)
}

AA_cts_verified_SRP <- SC_phobius_df %>%
  filter(`Experimental label` == "SRP-dependent") %>%
  pull(window_AA_string) %>%
  count_AAs(group_name = "Verified SRP-dependent")

AA_cts_verified_Sec63 <- SC_phobius_df %>%
  filter(`Experimental label` == "Sec63-dependent") %>%
  pull(window_AA_string) %>%
  count_AAs(group_name = "Verified Sec63-dependent")

AA_cts_phobius_TM <- SC_phobius_df %>%
  filter(window_type == "TM") %>%
  pull(window_AA_string) %>%
  count_AAs(group_name = "Phobius TM")

AA_cts_phobius_SP <- SC_phobius_df %>%
  filter(window_type == "SP") %>%
  pull(window_AA_string) %>%
  count_AAs(group_name = "Phobius SP")

AA_prob_df <-
  bind_rows(
    AA_cts_verified_SRP,
    AA_cts_verified_Sec63,
    AA_cts_phobius_TM,
    AA_cts_phobius_SP
  ) %>%
  pivot_longer(cols = -group, names_to = "AA", values_to = "count") %>%
  mutate(AA = factor(AA, levels = sorted_AA)) %>%
  group_by(group) %>%
  mutate(prob = count / sum(count))

# parallel coordinates plot for AA composition
AA_prob_plot <-
  ggplot(AA_prob_df, aes(x = AA, y = prob)) +
  geom_line(aes(group = group, colour = group),
    linewidth = 0.8
  ) +
  scale_colour_manual(
    breaks = c(
      "Verified SRP-dependent",
      "Verified Sec63-dependent",
      "Phobius TM",
      "Phobius SP"
    ),
    values = c(
      "Verified SRP-dependent" = "red",
      "Verified Sec63-dependent" = "blue",
      "Phobius TM" = "indianred",
      "Phobius SP" = "skyblue"
    )
  ) +
  scale_y_continuous(
    breaks = c(0, 0.1, 0.2),
    minor_breaks = c(0.05, 0.15, 0.25),
    expand = c(0, 0)
  ) +
  guides(colour = guide_legend(ncol = 2, nrow = 2)) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 8),
    plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "inches"),
    panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_line(colour = "grey90", linewidth = 0.25)
  ) +
  labs(
    x = "Amino Acid, ordered by KD hydrophobicity",
    y = "Proportion in hydrophobic helix"
  )

ggsave(
  filename = here("results", "figures", "AA_prob_plot.pdf"),
  plot = AA_prob_plot,
  width = 4.5, height = 3, dpi = 300
)

ggsave(
  filename = here("results", "figures", "AA_prob_plot.png"),
  plot = AA_prob_plot,
  width = 4.5, height = 3, dpi = 300, bg = "white"
)


# --- PSIPRED analysis --- #

# setwd for python script
setwd(here("src", "psipred"))

# source psipred extraction script
library(reticulate)
source_python(here("src", "psipred", "PSIPRED.py"))

# run psipred extraction
df <- psipred_df(8, 0)

# plots are done by phobius labels
phobius_srp <- read_lines(here("results", "proteins", "Saccharomyces cerevisiae_TM.txt"))
phobius_non_srp <- read_lines(here("results", "proteins", "Saccharomyces cerevisiae_SP.txt"))

labelled_df <- df %>%
  mutate(
    `Phobius Label` =
      case_when(
        seqid %in% phobius_srp ~ "SRP-dependent",
        seqid %in% phobius_non_srp ~ "Sec63-dependent",
        TRUE ~ "Unlabelled"
      )
  ) %>%
  filter(`Phobius Label` != "Unlabelled")

# plots with no minimum length and high confidence
labelled_df %>%
  ggplot(aes(x = window_length, fill = `Phobius Label`)) +
  geom_histogram(binwidth = 1) +
  labs(
    title = "Distribution of window length by Phobius label, confidence >= 8, length >= 0",
    x = "Window length",
    y = "Count"
  ) +
  theme(legend.position = "bottom")
labelled_df %>%
  ggplot(aes(x = window_length, fill = `Phobius Label`, alpha = 0.8)) +
  geom_histogram(binwidth = 1, position = "identity") +
  labs(
    title = "Distribution of structure prediction counts by Phobius label, confidence >= 8, length >= 0",
    x = "Window length",
    y = "Count"
  ) +
  theme(legend.position = "bottom")

# plot using count data - assumes all helix predicted amino acids
# are part of the same helix
df <- psipred_df(8, 0, count = TRUE)

labelled_df <- df %>%
  mutate(
    `Phobius Label` =
      case_when(
        seqid %in% phobius_srp ~ "SRP-dependent",
        seqid %in% phobius_non_srp ~ "Sec63-dependent",
        TRUE ~ "Unlabelled"
      )
  ) %>%
  filter(`Phobius Label` != "Unlabelled")

# recovers the classic distribution from previous analysis
labelled_df %>%
  ggplot(aes(x = window_length, fill = `Phobius Label`, alpha = 0.8)) +
  geom_histogram(binwidth = 1, position = "identity") +
  labs(
    title = "Distribution of structure prediction counts by Phobius label, confidence >= 8, length >= 0",
    x = "Window length",
    y = "Count"
  ) +
  theme(legend.position = "bottom")

# by labels, faceted by confidence and length
confidence_lower <- c(7, 8, 9)
length_lower <- c(0, 3, 5, 8)

# all pairs of confidence and length
pairs <- expand.grid(confidence_lower, length_lower)

# generate data
full_df <- data.frame()

for (i in 1:nrow(pairs)) {
  confidence <- pairs$Var1[i]
  length <- pairs$Var2[i]
  df <- psipred_df(confidence, length) %>%
    mutate(
      conf_lower = confidence,
      length_lower = length
    )

  full_df <- bind_rows(full_df, df)
}

labelled_full_df <- full_df %>%
  mutate(
    `Phobius Label` =
      case_when(
        seqid %in% phobius_srp ~ "SRP-dependent",
        seqid %in% phobius_non_srp ~ "Sec63-dependent",
        TRUE ~ "Unlabelled"
      )
  ) %>%
  filter(`Phobius Label` != "Unlabelled")

# interestingly, confidence = 7 has the valley in Sec63-dependent h-region
# length frequency like the phobius distributions
labelled_full_df %>%
  ggplot(aes(x = window_length, fill = `Phobius Label`, alpha = 0.8)) +
  geom_histogram(binwidth = 1, position = "identity") +
  labs(
    title = "Distribution of window length by Phobius label, confidence >= {7, 8, 9}, length >= {0, 3, 5, 8}",
    x = "Window length",
    y = "Count"
  ) +
  facet_grid(conf_lower ~ length_lower) +
  theme(legend.position = "bottom")
