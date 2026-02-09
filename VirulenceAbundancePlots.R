library(tidyverse)
library(ggsci) 
library(scales)
csv_file   <- "virulence_abundance.csv"
out_png    <- "virulence_abundance_plot.png"
sample_dat <- read_csv(csv_file, show_col_types = FALSE) %>%
  filter(Location %in% c("AtlOnly","BoOnly","InOnly")) %>%
  filter(!grepl("escherichia.coli|e\\.\\s*coli", Organism, ignore.case = TRUE)) %>%
  mutate(Sample_Category = paste(Distance_Category, Sample, sep = "_")) %>%
  group_by(Location, Sample_Category, Organism) %>%
  summarize(Total_Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
  group_by(Location, Sample_Category) %>%
  mutate(Relative_Abundance = Total_Abundance / sum(Total_Abundance) * 100) %>%
  ungroup()
sample_counts <- sample_dat %>%
  group_by(Location) %>%
  summarize(n = n_distinct(Sample_Category), .groups = "drop")
facet_labels <- sample_counts %>%
  mutate(
    label = case_when(
      Location == "AtlOnly" ~ paste0("Atlanta, US (n=", n, ")"),
      Location == "BoOnly"  ~ paste0("La Paz, Bolivia (n=", n, ")"),
      Location == "InOnly"  ~ paste0("Kanpur, India (n=", n, ")"),
      TRUE ~ as.character(Location)
    )
  ) %>%
  { set_names(.$label, .$Location) }
levels_order <- unique(sample_dat$Sample_Category)
sample_dat <- sample_dat %>%
  mutate(
    Location = factor(Location, levels = c("AtlOnly","BoOnly","InOnly")),
    Sample_Category = factor(Sample_Category, levels = levels_order),
    Sample_Label = Sample_Category %>%
      str_replace("less10m", "<10m") %>%
      str_replace("10to100m", "10-100m") %>%
      str_replace("100to1000m", "100-1000m") %>%
      str_replace("greater1000m", ">1000m") %>%
      str_replace_all("_", " ")
  )
x_label_map <- sample_dat %>%
  distinct(Sample_Category, Sample_Label) %>%
  { set_names(.$Sample_Label, .$Sample_Category) }
orgs <- unique(sample_dat$Organism)
n_colors <- length(orgs)
base_npg <- pal_npg("nrc", alpha = 0.6)(10)
org_colors <- colorRampPalette(base_npg)(n_colors)
pirate_cols <- set_names(org_colors, orgs)
p <- ggplot(sample_dat,
            aes(x = Sample_Category, y = Relative_Abundance, fill = Organism)) +
  geom_col(width = 0.6) +  # Narrower bars for more spacing
  facet_wrap(~ Location, scales = "free_x", nrow = 1,
             labeller = labeller(Location = facet_labels)
  ) +
  scale_x_discrete(labels = x_label_map, expand = expansion(add = 0.8)) +
  scale_fill_manual(values = pirate_cols) +
  labs(
    title = "Virulence Factor Relative Abundance by Location",
    x     = NULL,
    y     = "Relative Abundance (%)",
    fill  = ""
  ) +
  theme_minimal(base_family = "Calibri", base_size = 14) +
  theme(
    axis.text.x      = element_text(angle = 60, hjust = 1, vjust = 1, size = 9),
    axis.title.y     = element_text(size = 16),
    strip.text       = element_text(face = "bold", size = 16),
    strip.clip       = "off",  
    legend.text      = element_text(face = "italic", size = 14),
    plot.title       = element_text(hjust = 0.5, size = 18, face = "bold"),
    panel.background = element_rect(fill = "white"),
    panel.border     = element_rect(color = "black", fill = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.spacing    = unit(2, "lines"),
    plot.margin      = margin(t = 25, r = 15, b = 10, l = 10)  
  )
print(p)
ggsave(out_png, p, width = 20, height = 8, dpi = 300)
