library(tidyverse)
library(here)
library(treeio)
library(ggtree)

# plot helix length axis
lower <- 5
upper <- 34
helix_delim <- seq(lower, upper, 5)
helix_minor <- seq(lower, upper, 5)
helix_limits <- c(lower, upper)
scale_x_helix_length <-
  scale_x_continuous("Predicted helix length (AA)",
    breaks = helix_delim,
    limits = helix_limits,
    minor_breaks = helix_minor,
    expand = expansion(mult = 0, add = 0.6)
  )


read_phobius <- function(protein_AA_path) {
  # extract file name from path, replace .fasta with _out
  file_name <- gsub(".fasta", "", basename(protein_AA_path))

  # create output path
  out_path <- paste0(here("results", "phobius", file_name), ".csv")

  # check if output path exists, if it does, exit function
  if (file.exists(out_path)) {
    return(read_csv(out_path))
  } else {
    cat("No output file found")
  }
}

phobius_cladogram_plot <- function(species_table, tree_string) {
  # attach path to protein file names
  species_df <- here("data", "proteins", "pub", species_table) %>%
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

  # run phobius
  phobius_results <- lapply(protein_paths, read_phobius)

  for (i in 1:length(phobius_results)) {
    phobius_results[[i]] <- phobius_results[[i]] %>%
      filter(phobius_end != 0) %>%
      mutate(window_length = phobius_end - phobius_start + 1) %>%
      mutate(species = species_df$Nicename_splitline[i])
  }

  # join and reset row names
  phobius_df <- do.call(rbind, phobius_results)
  rownames(phobius_df) <- NULL

  # get GO analysis values
  species_file_names <- lapply(species_df$Filename, function(x) gsub(".fasta", "", x))

  GO_df <- data.frame(
    species = character(),
    prediction = character(),
    GO_term = character(),
    p_value = numeric()
  )
  for (species_file in species_file_names) {
    for (prediction in c("SP", "TM")) {
      prediction_file <- paste(species_file, prediction, sep = "_")
      path <- here("results", "GO", prediction_file, "goEnrichmentResult.tsv")

      # get lowest p-value GO term
      # after removing meaningless "cellular component"
      go_df <- read_tsv(path)
      go_df <- go_df %>%
        filter(Name != "cellular component") %>%
        filter(`P-value` == min(`P-value`))

      GO_df <- rbind(
        GO_df,
        data.frame(
          species = species_file,
          prediction = prediction,
          GO_term = go_df$Name,
          p_value = go_df$`P-value`
        )
      )
    }
  }

  # heights of histogram modes for GO label positions
  sub_figure_heights <- phobius_df %>%
    filter(window_length == 12) %>%
    group_by(species) %>%
    summarise(height = n())

  # create a dataframe with species, prediction, GO term and p-value
  GO_summary_df <- species_df %>%
    mutate(Filename = gsub(".fasta", "", Filename)) %>%
    left_join(GO_df, by = c("Filename" = "species")) %>%
    select(species = Nicename_splitline, prediction, GO_term, p_value) %>%
    mutate(pred_substr = paste(prediction, ": ", GO_term, sep = "")) %>%
    group_by(species) %>%
    summarise(GO_term = paste(pred_substr, collapse = "\n")) %>%
    # shorten longest GO term for display
    mutate(GO_term = stringr::str_remove(GO_term, pattern = "external ")) %>%
    left_join(sub_figure_heights, by = c("species" = "species"))

  phobius_plot <-
    ggplot(phobius_df, aes(x = window_length, fill = phobius_type)) +
    geom_histogram(binwidth = 1, center = 0) +
    geom_vline(xintercept = 13.5, linetype = "dashed") +
    geom_label(
      data = GO_summary_df,
      aes(x = 28, y = height %/% 1.4, label = GO_term),
      size = 2, inherit.aes = FALSE, show.legend = FALSE
    ) +
    facet_wrap(~species,
      scales = "free_y", ncol = 1,
      strip.position = "left"
    ) +
    scale_y_continuous("Number of proteins", position = "right") +
    scale_x_helix_length +
    scale_fill_manual("Phobius prediction",
      values = c("SP" = "skyblue3", "TM" = "indianred")
    ) +
    theme(
      legend.position = "bottom",
      strip.text.y.left = element_text(face = "italic", angle = 0),
      strip.placement = "outside"
    )

  # Make fungal species tree / cladogram to inform phobius plot
  # first define the tree in newick format, read in to treeio format
  fungaltree_data <- tree_string %>%
    textConnection() %>%
    read.newick()

  # Plot the tree using ggtree.
  # ladderize = FALSE preserves input tip order.
  fungaltree_plot <-
    ggtree(fungaltree_data,
      ladderize = FALSE
    ) +
    # geom_tiplab here would print the tip labels, useful for checking.
    # geom_tiplab() +
    scale_y_reverse()

  fungaltree_plot

  # make composite plot, including moving the y-axis title
  phobius_composite_plot <-
    plot_grid(
      fungaltree_plot +
        theme(plot.margin = margin(t = 0, r = 0, b = 0.55, l = 0, unit = "in")),
      phobius_plot,
      nrow = 1,
      rel_widths = c(0.2, 1)
    ) +
    theme(plot.background = element_rect(fill = "white", colour = NA))

  return(phobius_composite_plot)
}


truncate_proteins_fasta <- function(file, dirin, dirout, maxpos = 60) {
    # truncate all strings in an amino acid string set
    # unless they are too short
    myAAstrings <-
        readAAStringSet(
            paste(dirin, file, sep = "/")
        ) 
    long_enough <- ( width(myAAstrings) > maxpos )
    truncAAstrings <- subseq(myAAstrings[long_enough], 
                             end = maxpos)
    writeXStringSet(truncAAstrings, 
                    paste(dirout, file, sep = "/")
                    )
}
