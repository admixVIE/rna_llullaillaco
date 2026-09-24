basdir="data/"
'%ni%' <- Negate('%in%')
options("scipen"=100)
library(dplyr)
library(ggplot2)
library(forcats)
library(tidyr)
library(readr)
library(patchwork)
library(ggnewscale)
library(scales)
library(cowplot)


load(file=paste(basdir,"/mumstats.Robject",sep=""))
load(file=paste(basdir,"/mumstat_bestDNA.Robject",sep=""))
load(file=paste(basdir,"/mumstat_best.Robject",sep=""))

### plot of stats
mump<-mumstats[c(1:22),]
mump[is.na(mump)]<-0
mump[5,16]<-round(mump[5,16])
xnam<-c("Production reads (M)","Uniquely mapped (%)", "Assigned in mapped (%)","Spliced in assigned (%)", "rRNA in assigned (%)", "tRNAs","miRNAs")
hli<-c("darkorange","darkcyan","darkcyan")
names(hli)<-c("La Doncella_Arm","El Nino_Anus","El Nino_Anus_1_a")

mumpl<-data.frame(samples=rownames(mump),input=as.numeric(mump[,1])/1000000,umap=as.numeric(gsub("%","",mump[,3])),assign=as.numeric(mump[,16]),spliced=as.numeric(mump[,17]),rRNA=as.numeric(mump[,18]),tRNA=as.numeric(mump[,14]),miRNA=as.numeric(mump[,13]))

mumd<-data.frame(names=c("samples","input","umap","assign","spliced","rRNA","tRNA","miRNA"),
                 values=c("DNA",
                          median(as.numeric(mumstatsd[,1]))/100000,
                          median(as.numeric(gsub("%","",mumstatsd[,3]))),
                          median(as.numeric(mumstatsd[,14])),
                          median(as.numeric(mumstatsd[,15])),
                          median(as.numeric(mumstatsd[,16])),
                          0,0))

legn  <- c("El Nino anal swab" = "darkcyan", "La Doncella skin" = "darkorange", "DNA library" = "darkred")



# Metrics (unmap removed)
metrics <- c("input","umap","assign","spliced","rRNA","tRNA","miRNA")

# x labels
xlab_map <- xnam
if (is.null(names(xlab_map))) names(xlab_map) <- metrics

# Long format
long <- mumpl %>%
  pivot_longer(all_of(metrics), names_to = "metric", values_to = "value") %>%
  mutate(metric = factor(metric, levels = metrics))

# DNA medians (assign/spliced/rRNA only)
mumd_lines <- mumd %>%
  filter(names %in% c("assign","spliced","rRNA")) %>%
  transmute(metric = factor(names, levels = metrics),
            value  = readr::parse_number(values))
# Scale to percent if needed
if (max(mumd_lines$value, na.rm = TRUE) <= 1 &&
    max(long$value[long$metric %in% c("assign","spliced","rRNA")], na.rm = TRUE) > 1) {
  mumd_lines$value <- 100 * mumd_lines$value
}

base_theme <- theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "white", color = NA),
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

# Panel 1: Input (0–110M)
p_input <- ggplot(filter(long, metric == "input"),
                  aes(x = metric, y = value)) +
  geom_violin(fill = "#b3cde3", color = "grey35", linewidth = 0.3,
              adjust = 1.0, scale = "width", draw_quantiles = 0.5) +
  geom_jitter(color = "grey10", width = 0.15, height = 0,
              size = 1.2, alpha = 0.9, show.legend = FALSE) +
  scale_y_continuous(limits = c(0, 110)) +
  scale_x_discrete(labels = function(x) xlab_map[x]) +
  labs(x = NULL, y = NULL) + base_theme

# Panel 2: umap (0–1)
p_umap <- ggplot(filter(long, metric == "umap"),
                 aes(x = metric, y = value)) +
  geom_violin(fill = "#b3cde3", color = "grey35", linewidth = 0.3,
              adjust = 1.0, scale = "width", draw_quantiles = 0.5) +
  geom_jitter(color = "grey10", width = 0.15, height = 0,
              size = 1.2, alpha = 0.9, show.legend = FALSE) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_discrete(labels = function(x) xlab_map[x]) +
  labs(x = NULL, y = NULL) + base_theme +
  theme(
    plot.margin = margin(t = 5.5, r = 5.5, b = 5.5, l = 15, unit = "pt")
  )

# Legend data for Percent panel
legend_names <- setdiff(names(legn), "DNA library")
legend_points <- tibble::tibble(
  metric = factor(rep("assign", length(legend_names)), levels = metrics),
  value  = rep(49, length(legend_names)),
  label  = legend_names
)

# Panel 3: Percent (assign/spliced/rRNA), 0–50, legend only for DNA library
p_percent <- ggplot(filter(long, metric %in% c("assign","spliced","rRNA")),
                    aes(x = metric, y = value)) +
  geom_violin(fill = "#b3cde3", color = "grey35", linewidth = 0.3,
              adjust = 1.0, scale = "width", draw_quantiles = 0.5) +
  geom_jitter(color = "grey10", width = 0.15, height = 0,
              size = 1.2, alpha = 0.9, show.legend = FALSE) +
  geom_crossbar(data = transform(mumd_lines, label = "DNA library"),
                aes(x = metric, y = value, ymin = value, ymax = value, color = label),
                inherit.aes = FALSE, width = 0.9, linewidth = 0.2) +
  scale_color_manual(values = c("DNA library" = "darkred"),
                     breaks = "DNA library", name = NULL) +
  guides(color = guide_legend(override.aes = list(size = 2.5))) +
  scale_y_continuous(limits = c(0, 50)) +
  scale_x_discrete(labels = function(x) xlab_map[x]) +
  labs(x = NULL, y = NULL) +
  base_theme +
  theme(
    legend.position = c(0.98, 0.98),
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = scales::alpha("white", 0.85), color = NA),
    legend.key.size = grid::unit(3, "mm"),
    legend.text = element_text(size = 8),
    plot.margin = margin(t = 5.5, r = 5.5, b = 5.5, l = 15, unit = "pt")
  )

# Panel 4: Small RNAs (tRNA/miRNA), 0–6200
p_small <- ggplot(filter(long, metric %in% c("tRNA","miRNA")),
                  aes(x = metric, y = value)) +
  geom_violin(fill = "#b3cde3", color = "grey35", linewidth = 0.3,
              adjust = 1.0, scale = "width", draw_quantiles = 0.5) +
  geom_jitter(color = "grey10", width = 0.15, height = 0,
              size = 1.2, alpha = 0.9, show.legend = FALSE) +
  scale_y_continuous(limits = c(0, 6200)) +
  scale_x_discrete(labels = function(x) xlab_map[x]) +
  labs(x = NULL, y = NULL) + base_theme


## next panel: proportions for best sample
names(cats$Splicing)[1]<-"Unspliced"
names(cats$rRNA)[1]<-"Other genes"
names(catsd$Splicing)[1]<-"Unspliced"
names(catsd$rRNA)[1]<-"Other genes"
cat_levels <- c("Mapping", "Assignment", "Splicing", "rRNA")

# Tidy to proportions
tidy_props <- function(lst) {
  do.call(rbind, lapply(cat_levels, function(cat) {
    v <- lst[[cat]]
    if (is.null(v)) return(NULL)
    v <- unlist(v, use.names = TRUE)
    if (length(v) == 0) return(NULL)
    if (is.null(names(v))) names(v) <- paste0("C", seq_along(v))
    s <- sum(v, na.rm = TRUE)
    p <- if (s > 0) v / s else rep(0, length(v))
    data.frame(cat = cat, component = names(p), prop = as.numeric(p), stringsAsFactors = FALSE)
  })) |>
    transform(cat = factor(cat, levels = cat_levels))
}

df1 <- tidy_props(cats)
df2 <- tidy_props(catsd)

# Per-category palette function
cat_pal_fun <- function(cat) switch(cat,
                                    "Mapping"    = viridisLite::viridis,
                                    "Assignment" = viridisLite::magma,
                                    "Splicing"   = viridisLite::plasma,
                                    "rRNA"       = viridisLite::inferno
)

# Consistent palettes across cats and catsd (union of components)
pal_map <- lapply(cat_levels, function(cat) {
  comps <- sort(unique(c(df1$component[df1$cat == cat], df2$component[df2$cat == cat])))
  setNames(cat_pal_fun(cat)(length(comps)), comps)
})
names(pal_map) <- cat_levels

base_theme <- theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.background = element_rect(fill = "white", color = NA),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.background = element_rect(fill = "grey90", color = NA),
    strip.text = element_text(face = "bold"),
    legend.box = "vertical"
  )

build_plot <- function(df, pal_map) {
  p <- ggplot() +
    facet_wrap(~ cat, nrow = 1) +
    scale_y_continuous(limits = c(0, 1), labels = percent_format(accuracy = 1)) +
    labs(x = NULL, y = "Proportion") +
    base_theme
  
  for (i in seq_along(cat_levels)) {
    cat <- cat_levels[i]
    df_cat <- filter(df, cat == !!cat)
    if (i > 1) p <- p + ggnewscale::new_scale_fill()
    p <- p +
      geom_col(data = df_cat, aes(x = 1, y = prop, fill = component),
               width = 0.8, color = "grey20") +
      scale_fill_manual(values = pal_map[[cat]], name = cat,
                        guide = guide_legend(order = i))
  }
  p
}

p_left  <- build_plot(df1, pal_map)
p_right <- build_plot(df2, pal_map)

# Legend from left plot
leg <- cowplot::get_legend(p_left + theme(legend.position = "right"))

# Remove legends in panels
p_left_noleg  <- p_left  + theme(legend.position = "none")
p_right_noleg <- p_right + theme(legend.position = "none")



# 1) Build per-category legend grobs (using the same palettes pal_map)
make_cat_legend <- function(cat, pal) {
  comps <- names(pal)
  # dummy plot just to generate the legend
  p <- ggplot(data.frame(component = factor(comps, levels = comps))) +
    geom_bar(aes(x = 1, fill = component)) +
    scale_fill_manual(values = pal, name = cat) +
    theme_minimal(base_size = 10) +
    theme(
      legend.position = "right",
      legend.title = element_text(face = "bold"),
      legend.box.margin = margin(0, 0, 0, 0),
      legend.margin = margin(0, 0, 0, 0)
    )
  cowplot::get_legend(p)
}

leg_map  <- make_cat_legend("Mapping",    pal_map[["Mapping"]])
leg_ass  <- make_cat_legend("Assignment", pal_map[["Assignment"]])
leg_spl  <- make_cat_legend("Splicing",   pal_map[["Splicing"]])
leg_rrna <- make_cat_legend("rRNA",       pal_map[["rRNA"]])

# Legend row: categories next to each other
legend_row <- cowplot::plot_grid(
  cowplot::ggdraw(leg_map),
  cowplot::ggdraw(leg_ass),
  cowplot::ggdraw(leg_spl),
  cowplot::ggdraw(leg_rrna),
  ncol = 4, rel_widths = c(1, 1, 1, 1)
)

### fragment misincorporation patterns for RNA
md <- read.delim(paste(basdir,"/misincorporation.txt",sep=""), sep="\t",comment.char="#",
  check.names = FALSE
)

subs <- grep("^[ACGT]>[ACGT]$", names(md), value = TRUE)

damage <- md |>
  mutate(Pos = as.integer(Pos)) |>
  group_by(End, Pos) |>
  summarise(
    across(
      all_of(c("A", "C", "G", "T", subs)),
      \(x) sum(as.numeric(x), na.rm = TRUE)
    ),
    .groups = "drop"
  ) |>
  pivot_longer(
    all_of(subs),
    names_to = "substitution",
    values_to = "mismatches"
  ) |>
  mutate(
    source = substr(substitution, 1, 1),
    denominator = case_when(
      source == "A" ~ A,
      source == "C" ~ C,
      source == "G" ~ G,
      source == "T" ~ T
    ),
    frequency = if_else(
      denominator > 0,
      mismatches / denominator,
      NA_real_
    ),
    distance = abs(Pos),
    end = factor(
      End,
      levels = c("5p", "3p"),
      labels = c("5' end", "3' end")    ),
    substitution = gsub(">", "\u2192", substitution, fixed = TRUE)
  )

damage25 <- damage |>
  mutate(
    sub = gsub("\u2192", ">", substitution, fixed = TRUE),
    end = factor(
      End,
      levels = c("5p", "3p"),
      labels = c("5' end", "3' end")    ),
    # Places position 1 at the outside of each plot
    plot_position = if_else(End == "3p", -distance, distance)
  ) |>
  filter(distance >= 1, distance <= 25)

p_damage <- ggplot(
  damage25,
  aes(plot_position, frequency, group = sub)
) +
  # All other substitutions
  geom_line(
    data = filter(
      damage25,
      !sub %in% c("C>T", "G>A")
    ),
    colour = "grey70",
    linewidth = 0.35,
    na.rm = TRUE
  ) +
  # C>T
  geom_line(
    data = filter(damage25, sub == "C>T"),
    colour = "#D73027",
    linewidth = 1.1,
    na.rm = TRUE
  ) +
  # G>A
  geom_line(
    data = filter(damage25, sub == "G>A"),
    colour = "#2166AC",
    linewidth = 1.1,
    na.rm = TRUE
  ) +
  facet_wrap(
    ~end,
    nrow = 1,
    scales = "free_x"
  ) +
  scale_x_continuous(
    breaks = c(
      -25, -20, -15, -10, -5, -1,
      1, 5, 10, 15, 20, 25
    ),
    labels = function(x) abs(x),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  scale_y_continuous(
    breaks = seq(0, 0.20, 0.05),
    labels = label_percent(accuracy = 1),
    expand = expansion(mult = c(0, 0)),
    sec.axis = dup_axis(
      name = ""
    )
  ) +
  coord_cartesian(ylim = c(0, 0.20)) +
  labs(
    x = "Distance from read end (nt)",
    y = "Mismatch frequency"
  ) +
  theme_classic(base_size = 9) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "none",
    panel.spacing.x = unit(8, "pt"),
    plot.margin = margin(t = 5.5, r = 5.5, b = 5.5, l = 15, unit = "pt")
  )


## damage patterns for DNA!
mdD <- read.delim(paste(basdir,"/DNA_misincorporation.txt",sep=""), sep="\t",comment.char="#",
                 check.names = FALSE)

subsD <- grep("^[ACGT]>[ACGT]$", names(mdD), value = TRUE)

damageD <- mdD |>
  mutate(Pos = as.integer(Pos)) |>
  group_by(End, Pos) |>
  summarise(
    across(
      all_of(c("A", "C", "G", "T", subs)),
      \(x) sum(as.numeric(x), na.rm = TRUE)
    ),
    .groups = "drop"
  ) |>
  pivot_longer(
    all_of(subsD),
    names_to = "substitution",
    values_to = "mismatches"
  ) |>
  mutate(
    source = substr(substitution, 1, 1),
    denominator = case_when(
      source == "A" ~ A,
      source == "C" ~ C,
      source == "G" ~ G,
      source == "T" ~ T
    ),
    frequency = if_else(
      denominator > 0,
      mismatches / denominator,
      NA_real_
    ),
    distance = abs(Pos),
    end = factor(
      End,
      levels = c("5p", "3p"),
      labels = c("5' end", "3' end")    ),
    substitution = gsub(">", "\u2192", substitution, fixed = TRUE)
  )

damage25D <- damageD |>
  mutate(
    sub = gsub("\u2192", ">", substitution, fixed = TRUE),
    end = factor(
      End,
      levels = c("5p", "3p"),
      labels = c("5' end", "3' end")    ),
    # Places position 1 at the outside of each plot
    plot_position = if_else(End == "3p", -distance, distance)
  ) |>
  filter(distance >= 1, distance <= 25)

p_damageD <- ggplot(
  damage25D,
  aes(plot_position, frequency, group = sub)
) +
  # All other substitutions
  geom_line(
    data = filter(
      damage25D,
      !sub %in% c("C>T", "G>A")
    ),
    colour = "grey70",
    linewidth = 0.35,
    na.rm = TRUE
  ) +
  # C>T
  geom_line(
    data = filter(damage25D, sub == "C>T"),
    colour = "#D73027",
    linewidth = 1.1,
    na.rm = TRUE
  ) +
  # G>A
  geom_line(
    data = filter(damage25D, sub == "G>A"),
    colour = "#2166AC",
    linewidth = 1.1,
    na.rm = TRUE
  ) +
  facet_wrap(
    ~end,
    nrow = 1,
    scales = "free_x"
  ) +
  scale_x_continuous(
    breaks = c(
      -25, -20, -15, -10, -5, -1,
      1, 5, 10, 15, 20, 25
    ),
    labels = function(x) abs(x),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  scale_y_continuous(
    breaks = seq(0, 0.20, 0.05),
    labels = label_percent(accuracy = 1),
    expand = expansion(mult = c(0, 0)),
    sec.axis = dup_axis(
      name = "Mismatch frequency"
    )
  ) +
  coord_cartesian(ylim = c(0, 0.20)) +
  labs(
    x = "Distance from read end (nt)",
    y = ""
  ) +
  theme_classic(base_size = 9) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "none",
    panel.spacing.x = unit(8, "pt")
  )



# 2) First row: A–D
first_row <- plot_grid(
  p_input, p_umap, p_percent, p_small,
  ncol = 4, rel_widths = c(1, 1, 2, 1),
  labels = c("A", "B", "C", "D"),
  label_size = 12, label_fontface = "bold",
  label_x = 0.02, label_y = 0.98, hjust = 0, vjust = 1
)

# 3) Second row: E–F (only barplots)
second_row <- plot_grid(
  p_left_noleg, p_right_noleg,
  ncol = 2, rel_widths = c(1, 1),
  labels = c("E", "F"),
  label_size = 12, label_fontface = "bold",
  label_x = 0.02, label_y = 0.98, hjust = 0, vjust = 1
)

# Third row: G–H (RNA and DNA damage patterns)
third_row <- plot_grid(
  p_damage, p_damageD,
  ncol = 2,
  rel_widths = c(1, 1),
  labels = c("G", "H"),
  label_size = 12,
  label_fontface = "bold",
  label_x = 0.02,
  label_y = 0.98,
  hjust = 0,
  vjust = 1
)

# Combine plot rows and shared legend
final_fig <- plot_grid(
  first_row,
  second_row,
  legend_row,
  third_row,
  ncol = 1,
  rel_heights = c(1, 1,  0.55, 0.7),
  align = "v"
)


pdf("/files/Fig1_stats.pdf", width = 7.5, height = 8.5)
final_fig
dev.off()

