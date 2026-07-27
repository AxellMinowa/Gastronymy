##############################################################################
##                                                                          ##
##                         G A S T R O N O M E N C L A T U R E              ##
##            Supplementary Code S1 — Full Reproducible Analysis Script     ##
##                                                                          ##
##  Garraffoni et al.                                                       ##
##  Analysis of nomenclatural and etymological patterns in                  ##
##  Gastrotricha species names                                              ##
##                                                                          ##
##  Author       : Axell Minowa                                             ##
##  R version    : 4.5.1 (2025-06-13 ucrt)                                  ##
##  Last revised : 2026-07-25                                               ##
##                                                                          ##
##############################################################################
#

# CODE SECTION <-> MANUSCRIPT OUTPUT INDEX
# -----------------------------------------
#   Section 0  - Setup (packages, data, global plot parameters)
#   Section 1  - Data preprocessing
#   Section 2  - Basic descriptive summaries                  -> Table 1 (overview), Table S1 (environment)
#   Section 3  - Etymology category analysis                  -> Figure 2A
#   Section 4  - Epithet repetition & uniqueness tests         -> Figure 3
#   Section 5  - Name-length distributions & temporal trends   -> Figure Box (S)
#   Section 6  - Cumulative species description curve          -> Figure 1
#   Section 7  - Treemap of species per genus                  -> Figure S1
#   Section 8  - Authorship summary                            -> Figure 4
#   Section 9  - Etymology language origins                    -> Figure S2 (pie)
#   Section 10 - Honorific epithets & gender representation     -> Figure 6
#   Section 11 - GAM models: temporal etymology trends          -> Figure 2B, Table 2
#   Section 12 - Author collaboration network (build + stats)   -> Figure 5, Figures S3-S5, Table S2-S3
#   Section 13 - Nationality x species-name-language analysis   -> Figures S6-S8, Table S4
#   Section 14 - Nationality x etymology analysis                -> Figures S9-S16, Table S5
#   Section 15 - Consolidated figure export (main-text figures)
#   Section 16 - Consolidated table export
#   Section 17 - Session information
#
# REQUIRED INPUT FILES (working directory root)
# -----------------------------------------------
#   Final_dataset_2.csv  - one row per described species (genus, specific_epithet,
#                          Order, Family, year_of_publication, authors,
#                          number_of_authors, Etymologies_group,
#                          Etymologies_subgroup, Language, Environment, ...)
#   genus.csv            - genus-level reference table (";"-delimited, Windows-1252)
#   data.csv             - author & eponym metadata (surname, author_nationality,
#                          author_gender, eponyms, eponym_nationality, eponym_gender)
#   country_colors.csv   - country -> colour -> continent lookup for map/legend colours
#
# OUTPUT STRUCTURE
# -----------------
#   Figures/   - all PDF + PNG figures (main text: Figure1-6, FigureBox;
#                supplementary: supp_Figure1-16)
#   Tables/    - all CSV tables (summaries, network descriptors, CA coordinates)
#
################################################################################


# ==============================================================================
# 0. SETUP
# ==============================================================================

## 0.1 Package management -------------------------------------------------------

library(tidyverse)
library(readr)
library(tidyr)
library(dplyr)
library(purrr)
library(stringr)
library(scales)

library(tidygraph)
library(igraph)
library(ggraph)

library(patchwork)
library(gridExtra)
library(sjPlot)
library(paletteer)
library(ggsci)
library(treemapify)

library(mgcv)
library(mgcViz)
library(tidygam)
library(DHARMa)
library(performance)
library(FactoMineR)
library(factoextra)

set.seed(42)  # reproducibility for GAM REML fits, Louvain, and layout algorithms


## 0.2 Load data -----------------------------------------------------------------

df    <- read.csv("Final_dataset_2.csv")
genus <- read.csv("genus.csv", sep = ";", fileEncoding = "Windows-1252")
data  <- read.csv("data.csv")

country_continent <- read.csv("country_colors.csv", stringsAsFactors = FALSE)

# Author metadata (used throughout Sections 10, 12-14)
data_authors <- data %>%
  select(surname, author_nationality, author_gender)

# Eponym metadata (honorific epithets; Section 10)
data_eponyms <- data %>%
  select(eponyms, eponym_nationality, eponym_gender) %>%
  filter(!is.na(eponyms) & eponyms != "")

## 0.4 Global plot parameters ----------------------------------------------------

theme_custom <- function() {
  theme_bw() +
    theme(
      legend.position   = c(0.7, 0.2),
      legend.background = element_blank(),
      legend.title      = element_blank(),
      legend.text       = element_text(size = 11),
      strip.text        = element_text(size = 10, face = "bold"),
      axis.title        = element_text(size = 12),
      axis.text.x       = element_text(size = 11),
      axis.text.y       = element_text(size = 11),
      axis.ticks.length = unit(-0.75, "mm"),
      axis.ticks        = element_line(linewidth = 0.3),
      panel.grid        = element_blank(),
      plot.caption      = element_text(size = 10, color = "gray50"),
      plot.title        = element_text(face = "bold", size = 12)
    )
}
theme_set(theme_minimal())

aqua <- paletteer::paletteer_d("tvthemes::Aquamarine")[6]

country_colors <- setNames(country_continent$color, country_continent$country)

continent_order <- country_continent %>%
  arrange(continent, country) %>%
  pull(country)

country_continent$country <- factor(country_continent$country, levels = continent_order)

etym_labels <- c(
  M = "Morphological", G = "Geographical", P = "Personal",
  E = "Ecological",     C = "Cultural",     O = "Others"
)

year_range  <- range(df$year_of_publication, na.rm = TRUE)
year_breaks <- seq(floor(year_range[1] / 10) * 10, ceiling(year_range[2] / 10) * 10, by = 10)


# ==============================================================================
# 1. DATA PREPROCESSING
# ==============================================================================

df <- df %>%
  mutate(
    binomial          = paste(genus, specific_epithet),
    genus_length      = nchar(genus),
    epithet_length    = nchar(specific_epithet),
    binomial_length   = nchar(genus) + nchar(specific_epithet),
    number_of_authors = as.integer(number_of_authors),
    year              = as.numeric(year_of_publication),
    authors           = str_trim(authors),
    first_letter      = toupper(substr(specific_epithet, 1, 1)),
    # Rare/ambiguous etymology codes are pooled into "Others" (O)
    Etymologies_group = case_when(
      Etymologies_group %in% c("I", "F", "E / G") ~ "O",
      TRUE ~ Etymologies_group
    )
  ) %>%
  # Records with unresolved etymology codes are excluded from all downstream analyses
  filter(!Etymologies_group %in% c("?", "E / G"))


# ==============================================================================
# 2. BASIC DESCRIPTIVE SUMMARIES
# ==============================================================================

overview_summary <- tibble(
  Metric = c("Orders", "Families", "Genera", "Species", "Earliest year", "Latest year"),
  Value  = c(
    n_distinct(df$Order), n_distinct(df$Family), n_distinct(df$genus), nrow(df),
    min(df$year_of_publication, na.rm = TRUE), max(df$year_of_publication, na.rm = TRUE)
  )
)
cat("\n== BASIC SUMMARY STATISTICS ==\n"); print(overview_summary)

species_per_family <- df %>% count(Family, name = "species_count", sort = TRUE)

species_per_genus <- df %>%
  count(genus, name = "species_count", sort = TRUE) %>%
  mutate(percentage = round(100 * species_count / sum(species_count), 2)) %>%
  arrange(desc(percentage))
cat("\nTop 10 most species-rich genera:\n"); print(head(species_per_genus, 10))

freq_species_per_family <- species_per_family %>% count(species_count, name = "n_families") %>% arrange(species_count)
freq_species_per_genus  <- species_per_genus  %>% count(species_count, name = "n_genus")    %>% arrange(species_count)

# Environment representation (genus- and species-level); reported in Table S1
df_env_genus <- df %>%
  filter(!is.na(Environment), Environment != "") %>%
  count(genus, Environment) %>%
  pivot_wider(names_from = Environment, values_from = n, values_fill = 0) %>%
  mutate(total_env = rowSums(across(where(is.numeric)))) %>%
  arrange(desc(total_env))

df_env_species <- df %>%
  filter(!is.na(Environment), Environment != "") %>%
  distinct(species_name, Environment) %>%
  count(Environment, name = "n_species") %>%
  arrange(desc(n_species))


# ==============================================================================
# 3. ETYMOLOGY CATEGORY ANALYSIS  ->  Figure 2A
# ==============================================================================

etymologies_summary <- df %>%
  count(Etymologies_group, sort = TRUE) %>%
  mutate(label = etym_labels[Etymologies_group], pct = round(n / sum(n), 3))
cat("\nEtymology group distribution:\n"); print(etymologies_summary)

ety_group <- ggplot(etymologies_summary, aes(x = reorder(label, n), y = n, fill = label)) +
  geom_col(color = "black", alpha = 0.8) +
  coord_flip() +
  geom_text(aes(y = n + max(n) * 0.05, 
                label = paste0(n, "\n(", scales::percent(pct), ")")), 
            hjust = 0.1, size = 4, lineheight = 0.9) +
  scale_y_continuous(expand = expansion(mult = c(0.025, 0.2))) +
  labs(x = NULL, y = "Number of species") +
  theme_custom() +
  theme(legend.position = "none")

df_etym <- df %>%
  filter(!is.na(Etymologies_group), !is.na(Etymologies_subgroup)) %>%
  count(Etymologies_group, Etymologies_subgroup) %>%
  group_by(Etymologies_group) %>%
  mutate(total = sum(n)) %>%
  ungroup()

ety_subgr <- ggplot(df_etym, aes(x = reorder(Etymologies_group, total), y = n, fill = Etymologies_subgroup)) +
  geom_col(color = "black") +
  coord_flip() +
  theme_minimal() +
  theme(panel.grid = element_blank(), legend.position = c(0.8, 0.5), legend.title = element_blank()) +
  guides(fill = guide_legend(ncol = 2)) +
  labs(x = NULL, y = "Number of species")


# ==============================================================================
# 4. EPITHET REPETITION & UNIQUENESS TESTS  ->  Figure 3
# ==============================================================================

epithet_counts <- df %>% count(specific_epithet, name = "n")

epithet_leaderboard <- epithet_counts %>% count(n, name = "n_epithets") %>% arrange(desc(n))
cat("\nEpithet repetition summary:\n"); print(head(epithet_leaderboard, 10))

species_repetition <- epithet_counts %>%
  count(n, wt = n, name = "n_species") %>%
  mutate(percentage = 100 * n_species / sum(n_species)) %>%
  arrange(n)

epithet_repetition <- ggplot(epithet_leaderboard, aes(x = factor(n), y = n_epithets)) +
  geom_col(fill = aqua, color = "black") +
  labs(x = "Number of times an epithet is used", y = "Number of epithets") +
  theme_custom()

non_unique_epithets <- df %>%
  count(specific_epithet, Etymologies_group) %>%
  group_by(specific_epithet) %>%
  mutate(total_n = sum(n)) %>%
  ungroup() %>%
  filter(total_n > 1)
cat("\nMost frequently used epithets:\n")
print(non_unique_epithets %>% arrange(desc(total_n)) %>%
        select(specific_epithet, total_n, Etymologies_group) %>% distinct() %>% head(15))

epithet_freq <- ggplot(non_unique_epithets, aes(x = reorder(specific_epithet, total_n), y = n, fill = Etymologies_group)) +
  geom_col(color = "black") +
  coord_flip() +
  labs(x = "Epithet", y = "Count", fill = "Etymology") +
  theme_custom()

## 4.1 Epithet-uniqueness hypothesis tests --------------------------------------
# H: species-rich genera re-use epithets more often (non-unique epithets)

top5_genera <- species_per_genus %>% slice_max(species_count, n = 5) %>% pull(genus)

df_strat <- df %>%
  mutate(in_top5 = genus %in% top5_genera) %>%
  left_join(epithet_counts, by = "specific_epithet") %>%
  mutate(unique_epithet = n == 1) %>%
  left_join(species_per_genus %>% select(genus, species_count), by = "genus")

epithet_uniqueness_chisq <- chisq.test(table(df_strat$in_top5, df_strat$unique_epithet))
cat("\nChi-square test - epithet uniqueness vs. top-5 genus membership:\n")
print(epithet_uniqueness_chisq)

epithet_uniqueness_glm <- glm(unique_epithet ~ species_count, data = df_strat, family = binomial)
cat("\nGLM - epithet uniqueness ~ genus species count:\n")
print(summary(epithet_uniqueness_glm))

df_strat %>%
  group_by(in_top5) %>%
  summarise(pct_unique = mean(unique_epithet) * 100, n = n())

# Cramér's V for the chi-square test
tbl <- table(df_strat$in_top5, df_strat$unique_epithet)
sqrt(chisq.test(tbl)$statistic / sum(tbl))

# ==============================================================================
# 5. NAME-LENGTH DISTRIBUTIONS & TEMPORAL TRENDS  ->  Figure Box (Supplementary)
# ==============================================================================

plot_barchart <- function(data_in, column, title, xlab) {
  vals <- na.omit(data_in[[column]])
  stat <- tibble(mean = mean(vals), median = median(vals), sd = sd(vals))
  df_count <- tibble(x = vals) %>% count(x)
  ggplot(df_count, aes(x = x, y = n)) +
    geom_col(fill = aqua, color = "black") +
    geom_vline(xintercept = stat$mean,   linetype = "dashed", color = "red") +
    geom_vline(xintercept = stat$median, linetype = "solid",  color = "blue") +
    labs(title = title, x = xlab, y = "Count") +
    theme_custom()
}

bin_leng <- plot_barchart(df, "binomial_length", "Binomial Length", "Number of characters")
epi_leng <- plot_barchart(df, "epithet_length",  "Epithet Length",  "Number of characters")

avg_length_by_year <- df %>%
  filter(!is.na(year)) %>%
  group_by(year) %>%
  summarise(mean_length = mean(epithet_length), .groups = "drop")
overall_avg <- mean(df$epithet_length, na.rm = TRUE)

average_over_time <- ggplot(avg_length_by_year, aes(x = year, y = mean_length)) +
  geom_line(color = aqua, linewidth = 1) +
  geom_hline(yintercept = overall_avg, linetype = "dashed") +
  geom_smooth(method = "loess", se = FALSE, color = "gray50", linetype = "dotted") +
  labs(x = "Year", y = "Average epithet length") +
  theme_custom()

alpha_letter_freq <- df %>% count(first_letter)
alph_capital <- ggplot(alpha_letter_freq, aes(x = first_letter, y = n)) +
  geom_col(fill = aqua, color = "black") +
  labs(x = "Initial letter", y = "Count") +
  theme_custom()

mean_true   <- mean(df$epithet_length, na.rm = TRUE)
median_true <- median(df$epithet_length, na.rm = TRUE)

df_div <- df %>%
  group_by(genus) %>%
  summarise(
    mean_epithet_length = mean(epithet_length),
    species_count       = n(),
    earliest_year       = min(year_of_publication),
    .groups = "drop"
  )

diversity <- ggplot(df_div, aes(x = species_count, y = mean_epithet_length)) +
  geom_point(aes(color = earliest_year), size = 3, alpha = 0.7) +
  geom_hline(yintercept = mean_true,   linetype = "dashed", color = "red") +
  geom_hline(yintercept = median_true, linetype = "solid",  color = "blue") +
  scale_color_viridis_c(option = "plasma") +
  labs(x = "Number of species per genus", y = "Mean epithet length", color = "Earliest year") +
  theme_custom() +
  theme(legend.position = "right")


# ==============================================================================
# 6. CUMULATIVE SPECIES DESCRIPTION CURVE  ->  Figure 1
# ==============================================================================

df_cum <- df %>%
  filter(!is.na(year_of_publication)) %>%
  mutate(year = as.numeric(year_of_publication)) %>%
  count(year, name = "n_species") %>%
  arrange(year) %>%
  mutate(
    cumulative     = cumsum(n_species),
    total_species  = sum(n_species),
    cumulative_pct = cumulative / total_species
  )

milestones <- sapply(c(0.25, 0.5, 0.75), function(p) df_cum$year[min(which(df_cum$cumulative_pct >= p))])

cumm <- ggplot(df_cum, aes(x = year, y = cumulative)) +
  geom_area(fill = aqua, alpha = 0.5) +
  geom_line(color = aqua, linewidth = 1.2) +
  geom_vline(xintercept = milestones, linetype = "dashed") +
  annotate("text", x = milestones, y = df_cum$cumulative[match(milestones, df_cum$year)],
           label = paste0(c("25%", "50%", "75%"), "\n", milestones), vjust = -0.5, size = 3.5) +
  labs(x = "Year of Publication", y = "Cumulative Species Count") +
  theme_minimal()


# ==============================================================================
# 7. TREEMAP OF SPECIES PER GENUS  ->  Figure S1
# ==============================================================================

df_treemap <- df %>%
  group_by(Order, genus) %>%
  summarise(species_count = n(), earliest_year = min(year, na.rm = TRUE), .groups = "drop")

treemap <- ggplot(df_treemap, aes(area = species_count, fill = earliest_year, label = genus, subgroup = Order)) +
  geom_treemap() +
  geom_treemap_subgroup_border(color = "black") +
  geom_treemap_text(colour = "black", place = "center") +
  scale_fill_viridis_c(option = "plasma") +
  labs(fill = "Earliest year")

# ==============================================================================
# 8. AUTHORSHIP SUMMARY  ->  Figure 4
# ==============================================================================

authorship_summary <- df %>%
  summarise(
    mean_authors   = round(mean(number_of_authors), 2),
    median_authors = median(number_of_authors),
    max_authors    = max(number_of_authors),
    pct_solo       = mean(number_of_authors == 1),
    pct_multi      = mean(number_of_authors > 1)
  )
cat("\nAuthorship summary:\n"); print(authorship_summary)

df_clean <- df %>% distinct(species_name, year_of_publication, number_of_authors, .keep_all = TRUE)

number_of_authors_plot <- ggplot(df_clean, aes(x = year_of_publication, fill = factor(number_of_authors))) +
  geom_histogram(binwidth = 1, color = "black") +
  scale_x_continuous(limits = year_range, breaks = year_breaks) +
  labs(x = "Year", y = "Number of Publications", fill = "Number of Authors") +
  theme_minimal() +
  theme(panel.grid = element_blank())


# ==============================================================================
# 9. ETYMOLOGY LANGUAGE ORIGINS  ->  Figure S2
# ==============================================================================

lang_summary <- df %>% count(Language, sort = TRUE)
cat("\nLanguage distribution:\n"); print(lang_summary)

language_summary_grouped <- df %>%
  count(Language) %>%
  mutate(origin = case_when(
    Language %in% c("Latin", "Greek", "Latin+Greek") ~ "Latin or Greek",
    str_detect(Language, "Latinisation")              ~ "Latinisation",
    TRUE                                              ~ "Other"
  )) %>%
  group_by(origin) %>%
  summarise(n = sum(n), .groups = "drop") %>%
  mutate(pct = scales::percent(n / sum(n), accuracy = 0.1))

pie <- ggplot(language_summary_grouped, aes(x = 2, y = n, fill = origin)) +
  geom_col(color = "white") +
  coord_polar(theta = "y") +
  xlim(0.5, 2.5) +
  geom_text(aes(label = paste0(origin, "\n", pct)), position = position_stack(vjust = 0.5), size = 4) +
  labs(title = "Etymology Language Origins") +
  theme_void()


# ==============================================================================
# 10. HONORIFIC (PERSONAL) EPITHETS & GENDER REPRESENTATION  ->  Figure 6
# ==============================================================================

df_person <- df %>%
  filter(Etymologies_group == "P") %>%
  left_join(data_eponyms %>% select(eponyms, eponym_gender), by = c("specific_epithet" = "eponyms")) %>%
  rename(gender = eponym_gender) %>%
  filter(!is.na(gender), !is.na(year))

df_barplot <- df_person %>%
  filter(specific_epithet != "brahmsi") %>%
  count(year, gender) %>%
  mutate(n = ifelse(gender == "Female", -n, n))

gender <- ggplot(df_barplot, aes(x = year, y = n, fill = gender)) +
  geom_col(width = 0.8, color = "black") +
  geom_hline(yintercept = 0, color = "black") +
  labs(x = "Year", y = "Number of honorific epithets") +
  theme_minimal() +
  theme(panel.grid = element_blank())

eponyms_summary <- df_person %>% count(gender, sort = TRUE)
authors_summary <- data_authors %>% count(author_gender)
cat("\nGender distribution in eponyms:\n"); print(eponyms_summary)
cat("\nGender distribution in authors:\n"); print(authors_summary)

# First female-authored species description
female_authors <- df %>%
  separate_rows(authors, sep = ",\\s*") %>%
  mutate(authors = str_trim(authors)) %>%
  left_join(data_authors %>% select(surname, author_gender), by = c("authors" = "surname")) %>%
  filter(author_gender == "Female") %>%
  arrange(year) %>%
  select(species_name, year, authors)
cat("\nFirst female-authored species:\n"); print(head(female_authors, 5))

# First female honorific eponym
female_eponym <- df %>%
  filter(Etymologies_group == "P") %>%
  left_join(data_eponyms %>% select(eponyms, eponym_gender), by = c("specific_epithet" = "eponyms")) %>%
  filter(eponym_gender == "Female") %>%
  arrange(year) %>%
  select(species_name, year, specific_epithet, eponym_gender)
cat("\nFirst female eponym:\n"); print(head(female_eponym, 5))

# Two-sided and directional test: is the female-eponym rate lower than the female-authorship rate?
total_eponyms <- 107 + 25  # 132
total_authors <- 52 + 115  # 167

prop_test_two_sided <- prop.test(x = c(25, 52), n = c(total_eponyms, total_authors), alternative = "two.sided")
prop_test_directional <- prop.test(x = c(25, 52), n = c(total_eponyms, total_authors), alternative = "less")
cat("\nProportion test (two-sided) - female eponyms vs. female authors:\n"); print(prop_test_two_sided)
cat("\nProportion test (one-sided: female eponyms < female authors):\n"); print(prop_test_directional)

female <- tibble(
  group      = c("Authors", "Eponyms"),
  female_pct = c(52 / 167, 25 / 132) * 100,
  total      = c(167, 132)
) %>%
  ggplot(aes(x = group, y = female_pct)) +
  geom_col(width = 0.5, alpha = 0.8, color = "black", fill = aqua) +
  geom_text(aes(label = sprintf("%.1f%%\n(n=%d)", female_pct, total)), vjust = -0.5, size = 4, color = "grey40") +
  scale_y_continuous(limits = c(0, 45), labels = function(x) paste0(x, "%")) +
  labs(x = NULL, y = "% Female") +
  theme_minimal() +
  theme(legend.position = "none")


# ==============================================================================
# 11. GAM MODELS: TEMPORAL ETYMOLOGY TRENDS  ->  Figure 2B, Table 2
# =============================================================================

## 11.1 H0_1: temporal stability of etymology-category proportions -------------
# (quasibinomial GAM following Mammola et al. 2023)

db_year_plot <- df %>%
  filter(!is.na(year_of_publication), !is.na(Etymologies_group)) %>%
  count(year_of_publication, Etymologies_group) %>%
  complete(year_of_publication, Etymologies_group, fill = list(n = 0)) %>%
  group_by(year_of_publication) %>%
  mutate(Tot = sum(n), Value = n) %>%
  ungroup() %>%
  rename(Year = year_of_publication, Type = Etymologies_group) %>%
  filter(Year < 2020) %>%                                   # excludes incomplete final year
  mutate(Type = factor(Type, levels = names(etym_labels)))

r1_H0_1 <- gam(
  cbind(Value, Tot - Value) ~ Type + s(Year, by = Type, k = 10),
  family = binomial(link = "logit"), data = db_year_plot
)
cat("\nOverdispersion check (binomial GAM, H0_1):\n")
print(performance::check_overdispersion(r1_H0_1))         # expected: overdispersed

r2_H0_1 <- gam(
  cbind(Value, Tot - Value) ~ Type + s(Year, k = 10) + s(Year, by = Type, k = 10),
  family = quasibinomial(link = "logit"), data = db_year_plot
)
cat("\nFinal model R^2 (H0_1, quasibinomial):\n"); print(performance::r2(r2_H0_1))
cat("\nFinal model summary (H0_1):\n"); print(summary(r2_H0_1))

newdata_H0_1 <- expand.grid(
  Year = seq(min(db_year_plot$Year), max(db_year_plot$Year), length.out = 200),
  Type = unique(db_year_plot$Type)
)
pred_H0_1 <- predict(r2_H0_1, newdata = newdata_H0_1, type = "link", se.fit = TRUE)

df_pred_H0_1 <- newdata_H0_1 %>%
  mutate(
    fit  = plogis(pred_H0_1$fit),
    se   = pred_H0_1$se.fit,
    ymin = plogis(pred_H0_1$fit - 1.96 * se),
    ymax = plogis(pred_H0_1$fit + 1.96 * se)
  )

trends <- ggplot(df_pred_H0_1, aes(x = Year, y = fit, color = Type, fill = Type)) +
  geom_line(size = 1.1) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha = 0.2, color = NA) +
  labs(
    y = "Estimated proportion",
    #title = "Temporal trends in etymology groups"
  ) +
  theme_minimal()

trends_observed <- ggplot(db_year_plot, aes(x = Year, y = Value / Tot, color = Type)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "gam", formula = y ~ s(x),
              method.args = list(family = quasibinomial(link = "logit")), se = TRUE) +
  labs(x = "Year", y = "Observed proportion") +
  theme_minimal() +
  scale_color_manual(values = paletteer::paletteer_d("ggsci::default_jco", 6))

## 11.2 H0_2 (final): per-category temporal trend, phylum-wide -----------------
# Tests whether each etymology category's usage probability changes over time,
# modelled independently per category with no Order-level structure.

df_taxa <- df %>%
  filter(!is.na(Etymologies_group), !is.na(year_of_publication)) %>%
  mutate(year_scaled = scale(year_of_publication)[, 1])

categories <- unique(df_taxa$Etymologies_group)

models_H0_2 <- map(categories, function(code) {
  df_temp <- df_taxa %>% mutate(target = as.numeric(Etymologies_group == code))
  gam(target ~ s(year_scaled, k = 8), family = binomial(link = "logit"),
      data = df_temp, method = "REML")
})
names(models_H0_2) <- categories

pvals_H0_2  <- map_dbl(models_H0_2, ~ summary(.x)$s.table[1, "p-value"])
edf_H0_2    <- map_dbl(models_H0_2, ~ summary(.x)$s.table[1, "edf"])
devexp_H0_2 <- map_dbl(models_H0_2, ~ summary(.x)$dev.expl)

cat("\nH0_2 - per-category temporal smooth p-values:\n");           print(pvals_H0_2)
cat("\nH0_2 - per-category effective degrees of freedom (edf):\n"); print(edf_H0_2)
cat("\nH0_2 - per-category deviance explained:\n");                 print(devexp_H0_2)

# Directional check: does each category's predicted probability rise or fall
# between the earliest and latest years in the dataset?
year_endpoints <- range(df_taxa$year_of_publication, na.rm = TRUE)
year_scaled_endpoints <- (year_endpoints - mean(df_taxa$year_of_publication, na.rm = TRUE)) /
  sd(df_taxa$year_of_publication, na.rm = TRUE)

direction_H0_2 <- map_dfr(names(models_H0_2), function(code) {
  m <- models_H0_2[[code]]
  pred <- predict(m, newdata = data.frame(year_scaled = year_scaled_endpoints), type = "response")
  tibble(
    category   = code,
    prob_early = pred[1],
    prob_late  = pred[2],
    direction  = ifelse(prob_late > prob_early, "increasing", "decreasing"),
    change     = prob_late - prob_early
  )
})
cat("\nH0_2 - direction of change (earliest -> latest year):\n"); print(direction_H0_2)

# Diagnostic smooth plots (one panel per etymology category)
pdf("Figures/supp_FigureS0_GAM_smooths_H0_2.pdf", width = 12, height = 8)
par(mfrow = c(2, 3))
walk(names(models_H0_2), function(code) {
  plot(models_H0_2[[code]], shade = TRUE, main = code, ylab = "logit(probability)")
})
par(mfrow = c(1, 1))
dev.off()


# ==============================================================================
# 12. AUTHOR COLLABORATION NETWORK  ->  Figure 5, Figures S3-S5, Tables S2-S3
# ==============================================================================

## 12.1 Build node and edge tables (single build, reused throughout) -----------

species_authors <- df %>%
  mutate(species_id = row_number(), author_list = str_split(authors, ",\\s*")) %>%
  select(species_id, author_list)

team_size_df <- species_authors %>%
  mutate(team_size = lengths(author_list)) %>%
  count(team_size, name = "n_species")

node_df <- species_authors %>%
  unnest(author_list) %>%
  rename(author = author_list) %>%
  mutate(author = str_trim(author)) %>%
  filter(author != "") %>%
  count(author, name = "n_species") %>%
  left_join(data_authors %>% select(surname, author_nationality, author_gender),
            by = c("author" = "surname")) %>%
  mutate(
    author_nationality = ifelse(is.na(author_nationality) | author_nationality == "", "Unknown", author_nationality),
    author_nationality = factor(author_nationality, levels = continent_order),
    author_gender      = ifelse(is.na(author_gender), "Unknown", author_gender)
  )

edge_df <- species_authors %>%
  filter(lengths(author_list) >= 2) %>%
  mutate(author_list = lapply(author_list, function(x) sort(str_trim(x)))) %>%
  unnest(author_list) %>%
  rename(author = author_list) %>%
  filter(author != "") %>%
  group_by(species_id) %>%
  summarise(pairs = list(combn(author, 2, simplify = FALSE)), .groups = "drop") %>%
  unnest(pairs) %>%
  mutate(from = sapply(pairs, `[[`, 1), to = sapply(pairs, `[[`, 2)) %>%
  select(from, to) %>%
  count(from, to, name = "n_collaborations")

g <- graph_from_data_frame(d = edge_df, directed = FALSE, vertices = node_df)

## 12.2 Global topology, connectivity & bibliometric descriptors ---------------

n_nodes   <- vcount(g)
n_edges   <- ecount(g)
density   <- round(edge_density(g), 6)

components_g <- components(g)
g_giant      <- induced_subgraph(g, which(components_g$membership == which.max(components_g$csize)))
avg_path     <- round(mean_distance(g_giant, directed = FALSE), 4)
diam         <- diameter(g_giant, directed = FALSE)
transitiv    <- round(transitivity(g, type = "global"), 4)

n_comp       <- components_g$no
comp_sizes   <- sort(components_g$csize, decreasing = TRUE)
giant_size   <- comp_sizes[1]
giant_pct    <- round(100 * giant_size / n_nodes, 1)
solo_authors <- sum(comp_sizes == 1)

cat("\n== CO-AUTHORSHIP NETWORK: GLOBAL TOPOLOGY ==\n")
cat(sprintf("Nodes (authors)        : %d\n", n_nodes))
cat(sprintf("Edges (collab. pairs)  : %d\n", n_edges))
cat(sprintf("Network density        : %.6f\n", density))
cat(sprintf("Diameter (giant comp.) : %d\n", diam))
cat(sprintf("Avg. path length       : %.4f\n", avg_path))
cat(sprintf("Global clustering      : %.4f\n", transitiv))
cat(sprintf("Connected components   : %d\n", n_comp))
cat(sprintf("Giant component        : %d nodes (%.1f%%)\n", giant_size, giant_pct))
cat(sprintf("Isolated (solo) authors: %d\n", solo_authors))

## 12.3 Node-level centrality ---------------------------------------------------

deg       <- degree(g)
strength_v <- strength(g, weights = E(g)$n_collaborations)
between   <- betweenness(g, normalized = TRUE)
close     <- closeness(g, normalized = TRUE)
eigen_c   <- eigen_centrality(g, weights = E(g)$n_collaborations)$vector
pg_rank   <- page_rank(g, weights = E(g)$n_collaborations)$vector
clust_loc <- transitivity(g, type = "local", isolates = "zero")

## 12.4 Degree distribution & scale-free test -----------------------------------

deg_freq    <- as.data.frame(table(degree = deg)) %>% mutate(degree = as.integer(as.character(degree)))
deg_nonzero <- deg[deg > 0]
pl_fit      <- fit_power_law(deg_nonzero, implementation = "plfit")
cat(sprintf("\nPower-law exponent (alpha): %.4f | x_min: %.1f | KS: %.4f\n",
            pl_fit$alpha, pl_fit$xmin, pl_fit$KS.stat))

## 12.5 Community structure ------------------------------------------------------

comm_louvain   <- cluster_louvain(g, weights = E(g)$n_collaborations)
comm_walktrap  <- cluster_walktrap(g, weights = E(g)$n_collaborations)
comm_labelprop <- cluster_label_prop(g, weights = E(g)$n_collaborations)

mod_louvain <- round(modularity(comm_louvain), 4)
mod_walk    <- round(modularity(comm_walktrap), 4)
mod_lp      <- round(modularity(comm_labelprop), 4)

cat("\n== COMMUNITY STRUCTURE ==\n")
cat(sprintf("Louvain    : %d communities | Q = %.4f\n", length(comm_louvain), mod_louvain))
cat(sprintf("Walktrap   : %d communities | Q = %.4f\n", length(comm_walktrap), mod_walk))
cat(sprintf("Label prop : %d communities | Q = %.4f\n", length(comm_labelprop), mod_lp))

node_stats <- tibble(
  author           = V(g)$name,
  n_species        = V(g)$n_species,
  degree           = deg,
  strength         = strength_v,
  betweenness_norm = round(between, 6),
  closeness_norm   = round(close, 6),
  eigenvector      = round(eigen_c, 6),
  pagerank         = round(pg_rank, 6),
  local_clustering = round(clust_loc, 4),
  community_louvain = membership(comm_louvain)[V(g)$name]
) %>% arrange(desc(n_species))

cat("\nTop 15 authors by betweenness (bridges):\n")
print(node_stats %>% arrange(desc(betweenness_norm)) %>%
        select(author, n_species, degree, betweenness_norm) %>% head(15))

## 12.6 Assortativity -------------------------------------------------------------

assort_degree  <- round(assortativity_degree(g, directed = FALSE), 4)
assort_species <- round(assortativity(g, V(g)$n_species, directed = FALSE), 4)
cat(sprintf("\nDegree assortativity: %.4f | Species-count assortativity: %.4f\n",
            assort_degree, assort_species))

## 12.7 Bibliometric descriptors --------------------------------------------------

n_species_total <- nrow(df)
n_solo          <- sum(lengths(species_authors$author_list) == 1)
solo_rate       <- round(100 * n_solo / n_species_total, 1)
mean_team       <- round(mean(lengths(species_authors$author_list)), 2)
max_team        <- max(lengths(species_authors$author_list))

cat(sprintf("\nSolo-authored species: %d (%.1f%%) | Mean team size: %.2f | Max team size: %d\n",
            n_solo, solo_rate, mean_team, max_team))

## 12.8 Visualisations ------------------------------------------------------------

# Figure 5 (main text): gender-coded main collaboration network
p_network <- ggraph(g, layout = "fr") +
  geom_edge_link(aes(width = n_collaborations, alpha = n_collaborations), colour = "brown1") +
  scale_edge_width(range = c(1, 4), name = "Collaborations") +
  scale_edge_alpha(range = c(0.2, 0.9), guide = "none") +
  geom_node_point(aes(size = n_species, shape = author_gender), fill = aqua, alpha = 1) +
  scale_shape_manual(values = c("Female" = 21, "Male" = 22, "Unknown" = 23), name = "Gender") +
  scale_size(range = c(5, 15), name = "Species described") +
  guides(colour = "none") +
  geom_node_label(aes(label = ifelse(n_species >= quantile(n_species, 0.75), name, NA), size = n_species),
                  repel = TRUE, size = 3.5, label.size = 0.15, fill = "white", alpha = 0.85, show.legend = FALSE) +
  theme_graph(base_family = "sans") +
  theme(legend.position = "none", panel.border = element_rect(color = "black", fill = NA))

# Figure S2 (supplementary): nationality-coded network
supp_p_network <- ggraph(g, layout = "fr") +
  geom_edge_link(aes(width = n_collaborations, alpha = n_collaborations), colour = "red") +
  scale_edge_width(range = c(1, 4), name = "Collaborations") +
  scale_edge_alpha(range = c(0.2, 0.9), guide = "none") +
  geom_node_point(aes(size = n_species, colour = author_nationality, shape = author_gender), alpha = 1) +
  scale_size(range = c(2, 14), name = "Species described") +
  scale_colour_manual(
    values = country_colors, name = "Nationality", breaks = continent_order,
    labels = function(x) {
      continent_of <- setNames(country_continent$continent, as.character(country_continent$country))
      cont <- continent_of[x]
      prev_cont <- c(NA, head(cont, -1))
      ifelse(is.na(prev_cont) | cont != prev_cont, paste0("[", cont, "]\n", x), x)
    }, drop = FALSE
  ) +
  geom_node_label(aes(label = ifelse(n_species >= quantile(n_species, 0.05), name, NA), size = n_species),
                  repel = TRUE, size = 3, label.size = 0.15, fill = "white", alpha = 0.85, show.legend = FALSE) +
  theme_graph(base_family = "sans") +
  theme(legend.position = "right", panel.border = element_rect(color = "black", fill = NA))

# Figure S3: network coloured by Louvain community
tg <- as_tbl_graph(g) %>% activate(nodes) %>% mutate(community = as.factor(membership(comm_louvain)[name]))

p_net_community <- ggraph(tg, layout = "fr") +
  geom_edge_link(aes(width = n_collaborations, alpha = n_collaborations), colour = "#4E9AF1") +
  scale_edge_width(range = c(0.2, 4), name = "Collaborations") +
  scale_edge_alpha(range = c(0.15, 0.85), guide = "none") +
  geom_node_point(aes(size = n_species, colour = community), alpha = 0.85) +
  scale_size(range = c(1.5, 14), name = "Species described") +
  scale_colour_discrete(guide = "none") +
  geom_node_label(aes(label = ifelse(n_species >= quantile(n_species, 0.80), name, NA)),
                  repel = TRUE, size = 2.8, label.size = 0.12, fill = "white", alpha = 0.85, show.legend = FALSE) +
  labs(
    title    = "Co-authorship network of species descriptions",
    subtitle = "Node size = species count | Edge width = collaborations | Colour = Louvain community",
    caption  = sprintf("n = %d authors | %d collaboration pairs | %d communities | Q = %.3f",
                       n_nodes, n_edges, length(comm_louvain), mod_louvain)
  ) +
  theme_graph(base_family = "sans") +
  theme(plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 9, colour = "grey40"),
        plot.caption = element_text(size = 8, colour = "grey60"),
        legend.position = "right")

# Figure S4 components: degree distribution, team size, top-20 betweenness/degree
p_deg <- ggplot(deg_freq, aes(x = degree, y = Freq)) +
  geom_point(colour = "#4E9AF1", size = 2, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, colour = "#F4845F", linewidth = 0.8, linetype = "dashed") +
  scale_x_log10(labels = label_comma()) +
  scale_y_log10(labels = label_comma()) +
  annotation_logticks(sides = "bl", size = 0.3) +
  labs(title = "Degree distribution (log-log)",
       subtitle = sprintf("Power-law alpha = %.2f | KS = %.4f", pl_fit$alpha, pl_fit$KS.stat),
       x = "Degree (log)", y = "Frequency (log)") +
  theme_minimal(base_size = 11) + theme(plot.title = element_text(face = "bold", size = 11))

p_team <- ggplot(team_size_df, aes(x = factor(team_size), y = n_species)) +
  geom_col(fill = "#4E9AF1", alpha = 0.85, width = 0.7) +
  geom_text(aes(label = n_species), vjust = -0.4, size = 3, colour = "grey30") +
  labs(title = "Team-size distribution", x = "Number of authors per species", y = "Species count") +
  theme_minimal(base_size = 11) + theme(plot.title = element_text(face = "bold", size = 11))

p_between <- node_stats %>%
  arrange(desc(betweenness_norm)) %>% head(20) %>%
  mutate(author = reorder(author, betweenness_norm)) %>%
  ggplot(aes(x = author, y = betweenness_norm)) +
  geom_col(aes(fill = n_species), width = 0.7) +
  scale_fill_gradient(low = "#B5D4F4", high = "#185FA5", name = "Species") +
  coord_flip() +
  labs(title = "Top 20 authors - betweenness centrality", subtitle = "Bridge authors connecting sub-communities",
       x = NULL, y = "Betweenness (normalised)") +
  theme_minimal(base_size = 11) + theme(plot.title = element_text(face = "bold", size = 11))

p_degree <- node_stats %>%
  arrange(desc(degree)) %>% head(20) %>%
  mutate(author = reorder(author, degree)) %>%
  ggplot(aes(x = author, y = degree)) +
  geom_col(aes(fill = n_species), width = 0.7) +
  scale_fill_gradient(low = "#F5C4B3", high = "#993C1D", name = "Species") +
  coord_flip() +
  labs(title = "Top 20 authors - degree centrality", subtitle = "Authors with most unique collaborators",
       x = NULL, y = "Degree") +
  theme_minimal(base_size = 11) + theme(plot.title = element_text(face = "bold", size = 11))

p_stats <- (p_deg | p_team) / (p_between | p_degree) +
  plot_annotation(title = "Co-authorship network - descriptive statistics",
                  theme = theme(plot.title = element_text(face = "bold", size = 13)))

## 12.9 Save network figures -------------------------------------------------------

ggsave("Figures/supp_Figure3_network_by_community.pdf", p_net_community, width = 16, height = 11, device = cairo_pdf)
ggsave("Figures/supp_Figure3_network_by_community.png", p_net_community, width = 16, height = 11, dpi = 300)
ggsave("Figures/supp_Figure4_network_statistics.pdf",   p_stats,         width = 16, height = 12, device = cairo_pdf)
ggsave("Figures/supp_Figure4_network_statistics.png",   p_stats,         width = 16, height = 12, dpi = 300)

## 12.10 Export network tables -----------------------------------------------------

write_csv(node_stats, "Tables/network_nodes.csv")
write_csv(edge_df,    "Tables/network_edges.csv")

network_global_stats <- tibble(
  metric = c("n_nodes", "n_edges", "density", "diameter", "avg_path_length", "global_clustering",
             "n_components", "giant_component_size", "giant_component_pct", "isolated_authors",
             "degree_assortativity", "species_assortativity", "power_law_alpha", "power_law_KS",
             "modularity_louvain", "n_communities_louvain", "modularity_walktrap", "n_communities_walktrap",
             "n_species_total", "solo_rate_pct", "mean_team_size", "max_team_size"),
  value = c(n_nodes, n_edges, density, diam, avg_path, transitiv, n_comp, giant_size, giant_pct, solo_authors,
            assort_degree, assort_species, round(pl_fit$alpha, 4), round(pl_fit$KS.stat, 4),
            mod_louvain, length(comm_louvain), mod_walk, length(comm_walktrap),
            n_species_total, solo_rate, mean_team, max_team)
)
write_csv(network_global_stats, "Tables/network_global_stats.csv")


# ==============================================================================
# 13. NATIONALITY x SPECIES-NAME-LANGUAGE ANALYSIS  ->  Figures S5-S7, Table S4
# ==============================================================================

authors_meta <- data_authors %>% rename(nationality = author_nationality, gender = author_gender)

species_long <- df %>%
  select(species_name, Language, authors) %>%
  filter(!is.na(Language) & Language != "") %>%
  mutate(author_list = str_split(authors, ",\\s*")) %>%
  unnest(author_list) %>%
  rename(surname = author_list) %>%
  mutate(surname = str_trim(surname))

species_nat <- species_long %>%
  left_join(authors_meta %>% select(surname, nationality), by = "surname") %>%
  filter(!is.na(nationality))

MIN_N_LANG <- 5
nat_counts  <- species_nat %>% count(nationality) %>% filter(n >= MIN_N_LANG) %>% pull(nationality)
lang_counts <- species_nat %>% count(Language)    %>% filter(n >= MIN_N_LANG) %>% pull(Language)

species_dedup <- species_nat %>%
  filter(nationality %in% nat_counts, Language %in% lang_counts) %>%
  mutate(nationality = factor(nationality), Language = factor(Language)) %>%
  distinct(species_name, nationality, Language)

ctab_lang <- table(species_dedup$nationality, species_dedup$Language)
chi_lang  <- chisq.test(ctab_lang, simulate.p.value = TRUE, B = 2000)
cramers_v_lang <- sqrt(chi_lang$statistic / (sum(ctab_lang) * (min(dim(ctab_lang)) - 1)))
cat("\n== NATIONALITY x LANGUAGE ==\n"); print(chi_lang)
cat(sprintf("Cramer's V: %.4f\n", cramers_v_lang))

prop_tab_lang <- as.data.frame(prop.table(ctab_lang, margin = 1)) %>%
  rename(nationality = Var1, Language = Var2, proportion = Freq) %>%
  left_join(as.data.frame(ctab_lang) %>% rename(nationality = Var1, Language = Var2, count = Freq),
            by = c("nationality", "Language"))

nat_order_lang  <- names(sort(rowSums(ctab_lang), decreasing = TRUE))
lang_order_lang <- names(sort(colSums(ctab_lang), decreasing = TRUE))
prop_tab_lang <- prop_tab_lang %>%
  mutate(nationality = factor(nationality, levels = rev(nat_order_lang)),
         Language    = factor(Language, levels = lang_order_lang))

p_heat_lang <- ggplot(prop_tab_lang, aes(x = Language, y = nationality, fill = proportion)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = ifelse(count > 0, paste0(count, "\n(", percent(proportion, accuracy = 1), ")"), "")),
            size = 2.6, colour = "white", fontface = "bold", lineheight = 1.1) +
  scale_fill_gradient(low = "#B5D4F4", high = "#042C53", name = "Proportion\nwithin nationality",
                      labels = percent_format(accuracy = 1)) +
  labs(title = "Language preference by author nationality",
       subtitle = "Cell = species count | Colour = row proportion (within nationality)",
       x = "Species name language", y = "Author nationality") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 9), axis.text.y = element_text(size = 9),
        plot.title = element_text(face = "bold", size = 13), plot.subtitle = element_text(colour = "grey40", size = 9),
        legend.position = "right", panel.grid = element_blank())

p_bubble_lang <- ggplot(prop_tab_lang %>% filter(count > 0),
                        aes(x = Language, y = nationality, size = count, colour = proportion)) +
  geom_point(alpha = 0.85) +
  scale_size(range = c(2, 14), name = "Species count") +
  scale_colour_gradient(low = "#B5D4F4", high = "#042C53", name = "Proportion", labels = percent_format(accuracy = 1)) +
  labs(title = "Nationality x language bubble chart",
       subtitle = "Bubble size = species count | Colour = proportion within nationality",
       x = "Language", y = NULL) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 9), axis.text.y = element_text(size = 9),
        plot.title = element_text(face = "bold", size = 13), plot.subtitle = element_text(colour = "grey40", size = 9),
        panel.grid.major = element_line(colour = "grey92"))

ca_lang <- CA(ctab_lang, graph = FALSE)
cat("\nCA explained inertia (nationality x language):\n"); print(round(ca_lang$eig, 4))

p_ca_lang <- fviz_ca_biplot(
  ca_lang, repel = TRUE, col.row = "#185FA5", col.col = "#D85A30",
  shape.row = 17, shape.col = 16, map = "symmetric",
  title = "Correspondence analysis - nationality x language", ggtheme = theme_minimal(base_size = 11)
) +
  labs(subtitle = paste0("Dim 1: ", round(ca_lang$eig[1, 2], 1), "% inertia | Dim 2: ", round(ca_lang$eig[2, 2], 1), "% inertia\n",
                         "Triangles = nationalities | Circles = languages")) +
  theme(plot.title = element_text(face = "bold", size = 13), plot.subtitle = element_text(colour = "grey40", size = 9))

ggsave("Figures/supp_Figure5_natlang_heatmap.pdf", p_heat_lang,   width = 12, height = 7, device = cairo_pdf)
ggsave("Figures/supp_Figure5_natlang_heatmap.png", p_heat_lang,   width = 12, height = 7, dpi = 300)
ggsave("Figures/supp_Figure6_natlang_bubble.pdf",  p_bubble_lang, width = 12, height = 7, device = cairo_pdf)
ggsave("Figures/supp_Figure6_natlang_bubble.png",  p_bubble_lang, width = 12, height = 7, dpi = 300)
ggsave("Figures/supp_Figure7_natlang_CA_biplot.pdf", p_ca_lang,   width = 10, height = 8, device = cairo_pdf)
ggsave("Figures/supp_Figure7_natlang_CA_biplot.png", p_ca_lang,   width = 10, height = 8, dpi = 300)

write_csv(prop_tab_lang, "Tables/natlang_proportions.csv")
write_csv(as.data.frame(ca_lang$row$coord) %>% tibble::rownames_to_column("nationality"),
          "Tables/CA_lang_nationality_coords.csv")
write_csv(as.data.frame(ca_lang$col$coord) %>% tibble::rownames_to_column("Language"),
          "Tables/CA_lang_language_coords.csv")


# ==============================================================================
# 14. NATIONALITY x ETYMOLOGY ANALYSIS  ->  Figures S8-S16, Table S5
# ==============================================================================

MIN_N_ETYM <- 5

make_etym_long <- function(target_col) {
  df %>%
    select(species_name, authors, val = all_of(target_col)) %>%
    filter(!is.na(val), val != "") %>%
    mutate(author_list = str_split(authors, ",\\s*")) %>%
    unnest(author_list) %>%
    rename(surname = author_list) %>%
    mutate(surname = str_trim(surname)) %>%
    left_join(authors_meta %>% select(surname, nationality), by = "surname") %>%
    filter(!is.na(nationality)) %>%
    distinct(species_name, nationality, val)
}

long_group    <- make_etym_long("Etymologies_group")
long_subgroup <- make_etym_long("Etymologies_subgroup")

run_etym_analysis <- function(long_df, label, fig_offset, min_n = MIN_N_ETYM) {
  
  cat(sprintf("\n== NATIONALITY x ETYMOLOGY %s ==\n", label))
  
  keep_nat <- long_df %>% count(nationality) %>% filter(n >= min_n) %>% pull(nationality)
  keep_val <- long_df %>% count(val)         %>% filter(n >= min_n) %>% pull(val)
  
  filt <- long_df %>%
    filter(nationality %in% keep_nat, val %in% keep_val) %>%
    mutate(nationality = factor(nationality), val = factor(val))
  
  if (nlevels(filt$nationality) < 2 || nlevels(filt$val) < 2) {
    cat("  Not enough categories to analyse at this MIN_N_ETYM threshold.\n")
    return(invisible(NULL))
  }
  
  ctab <- table(filt$nationality, filt$val)
  chi  <- chisq.test(ctab, simulate.p.value = TRUE, B = 2000)
  v    <- sqrt(chi$statistic / (sum(ctab) * (min(dim(ctab)) - 1)))
  cat("Contingency table:\n"); print(ctab)
  cat("Chi-square test:\n"); print(chi)
  cat(sprintf("Cramer's V: %.4f\n", v))
  
  prop_long <- as.data.frame(prop.table(ctab, margin = 1)) %>%
    rename(nationality = Var1, val = Var2, proportion = Freq) %>%
    left_join(as.data.frame(ctab) %>% rename(nationality = Var1, val = Var2, count = Freq),
              by = c("nationality", "val"))
  
  nat_ord <- names(sort(rowSums(ctab), decreasing = TRUE))
  val_ord <- names(sort(colSums(ctab), decreasing = TRUE))
  prop_long <- prop_long %>%
    mutate(nationality = factor(nationality, levels = rev(nat_ord)), val = factor(val, levels = val_ord))
  
  p_heat <- ggplot(prop_long, aes(x = val, y = nationality, fill = proportion)) +
    geom_tile(colour = "white", linewidth = 0.6) +
    geom_text(aes(label = ifelse(count > 0, paste0(count, "\n(", percent(proportion, accuracy = 1), ")"), "")),
              size = 2.5, colour = "white", fontface = "bold", lineheight = 1.1) +
    scale_fill_gradient(low = "#B5D4F4", high = "#042C53", name = "Proportion\nwithin nationality",
                        labels = percent_format(accuracy = 1)) +
    labs(title = sprintf("Etymology %s preference by nationality", label),
         subtitle = "Cell = species count | Colour = row proportion (within nationality)",
         x = sprintf("Etymology %s", label), y = "Nationality") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 9), axis.text.y = element_text(size = 9),
          plot.title = element_text(face = "bold", size = 13), plot.subtitle = element_text(colour = "grey40", size = 9),
          legend.position = "right", panel.grid = element_blank())
  
  p_bubble <- ggplot(prop_long %>% filter(count > 0), aes(x = val, y = nationality, size = count, colour = proportion)) +
    geom_point(alpha = 0.85) +
    scale_size(range = c(2, 14), name = "Species count") +
    scale_colour_gradient(low = "#B5D4F4", high = "#042C53", name = "Proportion", labels = percent_format(accuracy = 1)) +
    labs(title = sprintf("Nationality x etymology %s bubble chart", label),
         subtitle = "Bubble size = species count | Colour = proportion within nationality",
         x = sprintf("Etymology %s", label), y = NULL) +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 9), axis.text.y = element_text(size = 9),
          plot.title = element_text(face = "bold", size = 13), plot.subtitle = element_text(colour = "grey40", size = 9),
          panel.grid.major = element_line(colour = "grey92"))
  
  ca_res <- CA(ctab, graph = FALSE)
  cat(sprintf("CA explained inertia (%s):\n", label)); print(round(ca_res$eig, 4))
  
  p_ca <- fviz_ca_biplot(
    ca_res, repel = TRUE, col.row = "#185FA5", col.col = "#D85A30",
    shape.row = 17, shape.col = 16, map = "symmetric",
    title = sprintf("CA biplot - nationality x etymology %s", label), ggtheme = theme_minimal(base_size = 11)
  ) +
    labs(subtitle = paste0("Dim 1: ", round(ca_res$eig[1, 2], 1), "% inertia | Dim 2: ", round(ca_res$eig[2, 2], 1), "% inertia\n",
                           "Triangles = nationalities | Circles = ", label, " categories")) +
    theme(plot.title = element_text(face = "bold", size = 13), plot.subtitle = element_text(colour = "grey40", size = 9))
  
  p_facet <- NULL
  if (label == "Subgroup") {
    prop_long_f <- prop_long %>%
      mutate(parent_group = str_extract(as.character(val), "^[A-Z]")) %>%
      filter(!is.na(parent_group))
    
    p_facet <- ggplot(prop_long_f %>% filter(count > 0), aes(x = val, y = proportion, fill = nationality)) +
      geom_col(position = "dodge", width = 0.7, alpha = 0.85) +
      facet_wrap(~ parent_group, scales = "free_x", nrow = 1) +
      scale_fill_brewer(palette = "Set2", name = "Nationality") +
      scale_y_continuous(labels = percent_format(accuracy = 1)) +
      labs(title = "Etymology subgroup proportions by nationality, faceted by parent group",
           subtitle = "Y-axis = proportion of species within that nationality",
           x = "Etymology subgroup", y = "Proportion") +
      theme_minimal(base_size = 10) +
      theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 8), strip.text = element_text(face = "bold", size = 10),
            plot.title = element_text(face = "bold", size = 12), plot.subtitle = element_text(colour = "grey40", size = 9),
            legend.position = "bottom")
  }
  
  slug <- tolower(label)
  ggsave(sprintf("Figures/supp_Figure%d_etym_%s_heatmap.pdf", fig_offset,     slug), p_heat,   width = 13, height = 7, device = cairo_pdf)
  ggsave(sprintf("Figures/supp_Figure%d_etym_%s_heatmap.png", fig_offset,     slug), p_heat,   width = 13, height = 7, dpi = 300)
  ggsave(sprintf("Figures/supp_Figure%d_etym_%s_bubble.pdf",  fig_offset + 1, slug), p_bubble, width = 13, height = 7, device = cairo_pdf)
  ggsave(sprintf("Figures/supp_Figure%d_etym_%s_bubble.png",  fig_offset + 1, slug), p_bubble, width = 13, height = 7, dpi = 300)
  ggsave(sprintf("Figures/supp_Figure%d_etym_%s_CA.pdf",      fig_offset + 2, slug), p_ca,     width = 10, height = 8, device = cairo_pdf)
  ggsave(sprintf("Figures/supp_Figure%d_etym_%s_CA.png",      fig_offset + 2, slug), p_ca,     width = 10, height = 8, dpi = 300)
  if (!is.null(p_facet)) {
    ggsave(sprintf("Figures/supp_Figure%d_etym_%s_faceted.pdf", fig_offset + 3, slug), p_facet, width = 16, height = 6, device = cairo_pdf)
    ggsave(sprintf("Figures/supp_Figure%d_etym_%s_faceted.png", fig_offset + 3, slug), p_facet, width = 16, height = 6, dpi = 300)
  }
  
  write_csv(prop_long, sprintf("Tables/etym_%s_proportions.csv", slug))
  write_csv(as.data.frame(ca_res$row$coord) %>% tibble::rownames_to_column("nationality"),
            sprintf("Tables/CA_%s_nationality_coords.csv", slug))
  write_csv(as.data.frame(ca_res$col$coord) %>% tibble::rownames_to_column(label),
            sprintf("Tables/CA_%s_%s_coords.csv", slug, slug))
  
  invisible(list(ctab = ctab, ca = ca_res, chi = chi, cramers_v = v))
}

results_group    <- run_etym_analysis(long_group,    label = "Group",    fig_offset = 8)
results_subgroup <- run_etym_analysis(long_subgroup, label = "Subgroup", fig_offset = 12)

if (!is.null(results_group) && !is.null(results_subgroup)) {
  etym_analysis_summary <- tibble(
    analysis        = c("Nationality x Etymology Group", "Nationality x Etymology Subgroup"),
    chi_sq          = c(results_group$chi$statistic, results_subgroup$chi$statistic),
    df              = c(results_group$chi$parameter, results_subgroup$chi$parameter),
    p_value         = c(results_group$chi$p.value, results_subgroup$chi$p.value),
    cramers_v       = round(c(results_group$cramers_v, results_subgroup$cramers_v), 4),
    CA_dim1_inertia = round(c(results_group$ca$eig[1, 2], results_subgroup$ca$eig[1, 2]), 2),
    CA_dim2_inertia = round(c(results_group$ca$eig[2, 2], results_subgroup$ca$eig[2, 2]), 2)
  )
  cat("\n== COMBINED ETYMOLOGY x NATIONALITY SUMMARY ==\n"); print(etym_analysis_summary)
  write_csv(etym_analysis_summary, "Tables/etym_analysis_summary.csv")
}


# ==============================================================================
# 15. CONSOLIDATED FIGURE EXPORT (MAIN-TEXT FIGURES)
# ==============================================================================

pdf("Figures/Figure1.pdf", width = 14, height = 10); grid.arrange(cumm, ncol = 1); dev.off()
png("Figures/Figure1.png", width = 1400, height = 1000, res = 150); grid.arrange(cumm, ncol = 1); dev.off()

pdf("Figures/Figure2.pdf", width = 14, height = 10); grid.arrange(ety_group, trends, ncol = 2); dev.off()
png("Figures/Figure2.png", width = 1400, height = 1000, res = 150); grid.arrange(ety_group, trends, ncol = 2); dev.off()

pdf("Figures/Figure3.pdf", width = 14, height = 10); grid.arrange(epithet_freq, epithet_repetition, ncol = 2); dev.off()

pdf("Figures/Figure4.pdf", width = 15, height = 10); grid.arrange(number_of_authors_plot, ncol = 1); dev.off()

pdf("Figures/Figure5.pdf", width = 15, height = 10); grid.arrange(p_network, ncol = 1); dev.off()

pdf("Figures/Figure6.pdf", width = 15, height = 10); grid.arrange(gender, female, ncol = 1); dev.off()

pdf("Figures/FigureBox.pdf", width = 15, height = 10)
lay_char <- rbind(c(1, 2), c(3, 3), c(4, 5))
grid.arrange(bin_leng, epi_leng, average_over_time, alph_capital, diversity, layout_matrix = lay_char)
dev.off()

pdf("Figures/supp_Figure1_treemap.pdf", width = 15, height = 10); grid.arrange(treemap, ncol = 1); dev.off()

pdf("Figures/supp_Figure2_network_by_nationality.pdf", width = 15, height = 10); grid.arrange(supp_p_network, ncol = 1); dev.off()

ggsave("Figures/supp_FigureS_language_origins_pie.pdf", pie, width = 8, height = 6, device = cairo_pdf)
ggsave("Figures/supp_FigureS_language_origins_pie.png", pie, width = 8, height = 6, dpi = 300)


# ==============================================================================
# 16. CONSOLIDATED TABLE EXPORT
# ==============================================================================

write_csv(overview_summary,        "Tables/overview_summary.csv")
write_csv(species_per_family,      "Tables/species_per_family.csv")
write_csv(species_per_genus,       "Tables/species_per_genus.csv")
write_csv(df_env_genus,            "Tables/environment_per_genus.csv")
write_csv(df_env_species,          "Tables/environment_per_species.csv")
write_csv(etymologies_summary,     "Tables/etymology_summary.csv")
write_csv(epithet_leaderboard,     "Tables/epithet_repetition_leaderboard.csv")
write_csv(lang_summary,            "Tables/language_summary.csv")
write_csv(authorship_summary,      "Tables/authorship_summary.csv")
write_csv(eponyms_summary,         "Tables/eponyms_gender_summary.csv")
write_csv(direction_H0_2,          "Tables/H0_2_direction_of_change.csv")

cat("\nAll figures saved to 'Figures/'; all tables saved to 'Tables/'.\n")


# ==============================================================================
# 17. SESSION INFORMATION
# ==============================================================================

cat("\n", rep("=", 60), "\n", sep = "")
cat("  ANALYSIS COMPLETE\n")
cat(rep("=", 60), "\n", sep = "")
sessionInfo()

