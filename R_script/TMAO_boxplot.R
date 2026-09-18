# This script is used to generate a boxplot that examines TMAO serum metabolite levels from the iPOP study
install.packages(c("dplyr", "ggplot2"))
library(dplyr)
library(ggplot2)

# Input and output files 
met_file <- file.path("input_data", "Metabolomics.csv")
meta_file <- file.path("input_data", "SubjectInfo.csv")
gut_file <- file.path("input_data", "gut_16s_abundance.txt")
output_dir <- file.path("output_data", "TMAO_log10_16S_overlap_results")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(2026)

# Import the input files
met <- read.csv(met_file, check.names = FALSE, stringsAsFactors = FALSE)
meta <- read.csv(meta_file, stringsAsFactors = FALSE, na.strings = c("", "NA"))
gut <- read.delim(gut_file, check.names = FALSE, stringsAsFactors = FALSE)
names(met)[1] <- "MetabolomicsSampleID"

# Identify subjects in the gut 16S Visit 01 cohort
gut_subject <- sub("-.*$", "", gut$SampleID)
gut_visit <- sub("-.*$", "", sub("^[^-]+-", "", gut$SampleID))
visit01_subjects <- unique(gut_subject[gut_visit == "01"])

# Prepare the overlapping metabolomics samples
time_name <- c(E11 = "Baseline", E12 = "2 min", E13 = "15 min",
               E14 = "30 min", E15 = "1 hour")
tmao_column <- "Trimethylamine N-oxide (TMAO)"

dat <- data.frame(
  SubjectID = sub("-E[0-9]+$", "", met$MetabolomicsSampleID),
  Code = sub("^.*-", "", met$MetabolomicsSampleID),
  TMAO_TIC = met[[tmao_column]]
) |>
  filter(Code %in% names(time_name), SubjectID %in% visit01_subjects) |>
  left_join(meta[, c("SubjectID", "IR_IS_classification")], by = "SubjectID") |>
  mutate(
    Group = factor(IR_IS_classification, levels = c("IR", "IS")),
    Time = factor(unname(time_name[Code]), levels = unname(time_name)),
    log10_TMAO = log10(TMAO_TIC)
  ) |>
  filter(!is.na(Group))

if (anyDuplicated(dat[c("SubjectID", "Code")])) stop("Duplicate samples found.")
if (anyNA(dat$log10_TMAO) || any(dat$TMAO_TIC <= 0)) stop("Invalid TMAO values.")

# Run paired Wilcoxon tests:Each post-exercise time is compared with E11 baseline within IR or IS.
post_codes <- c("E12", "E13", "E14", "E15")
tests <- bind_rows(lapply(levels(dat$Group), function(g) {
  bind_rows(lapply(post_codes, function(code) {
    base <- dat |> filter(Group == g, Code == "E11") |>
      select(SubjectID, baseline = log10_TMAO)
    post <- dat |> filter(Group == g, Code == code) |>
      select(SubjectID, post = log10_TMAO)
    paired <- inner_join(base, post, by = "SubjectID")
    fit <- wilcox.test(paired$post, paired$baseline, paired = TRUE,
                       exact = FALSE, correct = TRUE)
    data.frame(Group = g, Code = code, Time = unname(time_name[code]),
               n = nrow(paired), p_value = fit$p.value)
  }))
}))
tests$q_value_BH <- p.adjust(tests$p_value, method = "BH")

# Prepare condition, sample-size, and p-value labels
labels <- bind_rows(lapply(levels(dat$Group), function(g) {
  group_data <- filter(dat, Group == g)
  group_test <- filter(tests, Group == g)
  bind_rows(
    data.frame(Group = g, Time = "Baseline",
               Label = sprintf("Baseline\nn = %d", sum(group_data$Code == "E11"))),
    data.frame(Group = g, Time = group_test$Time,
               Label = sprintf("%s\nn = %d\np = %.3f",
                               group_test$Time, group_test$n, group_test$p_value))
  )
}))

plot_data <- dat |>
  mutate(TimeText = as.character(Time)) |>
  left_join(labels, by = c("Group", "TimeText" = "Time")) |>
  mutate(Label = factor(Label, levels = unique(labels$Label)))

# Generate the box-and-whisker plot
colors <- c("Baseline" = "#686868", "2 min" = "#4C78A8",
            "15 min" = "#59A14F", "30 min" = "#F28E2B",
            "1 hour" = "#B55A9D")
group_labels <- c(IR = "IR (n = 12)", IS = "IS (n = 6)")

p <- ggplot(plot_data, aes(Label, log10_TMAO, fill = TimeText)) +
  geom_boxplot(outlier.shape = NA, width = .58, alpha = .53,
               linewidth = .6, color = "#333333") +
  geom_point(aes(color = TimeText),
             position = position_jitter(width = .10, seed = 2026),
             size = 2, alpha = .72) +
  facet_wrap(~Group, scales = "free_x", labeller = as_labeller(group_labels)) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) +
  scale_y_continuous(expand = expansion(mult = c(.03, .07))) +
  labs(x = NULL, y = expression(log[10] * "(TMAO TIC)")) +
  theme_minimal(base_size = 12) +
  theme(
    strip.text = element_text(face = "bold", size = 14),
    axis.text.x = element_text(size = 9, lineheight = 1.08),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position = "none",
    plot.margin = margin(10, 18, 10, 18)
  )

# Save the graph and result tables
ggsave(file.path(output_dir, "TMAO_log10_boxplots_16S_overlap_clean.png"),
       p, width = 14, height = 7, dpi = 300, bg = "white")
write.csv(dat, file.path(output_dir, "TMAO_log10_subject_data.csv"), row.names = FALSE)
write.csv(tests, file.path(output_dir, "TMAO_log10_paired_tests.csv"), row.names = FALSE)