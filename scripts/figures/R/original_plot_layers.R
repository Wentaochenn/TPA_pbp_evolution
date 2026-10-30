# Plot-layer excerpts from the supplied manuscript script.
# Each function evaluates the original plotting expression in a prepared data
# environment. Any later aesthetic amendment is explicit in revised_figures.R.

y.theme.strip <- theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
y.theme.strip.partial <- theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
x.theme.strip <- theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
x.theme.axis.rotate <- theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
x.theme.axis.rotate.angle <- theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
y.facet.text.rotate <- theme(strip.text.y = element_text(angle = 0, face = "bold"))
legend.strip <- theme(legend.position = "none")
legend.size <- theme(legend.key.size = unit(0.55, "line"))
standard.textsize <- 9
text.size.within <- 2.4
theme.text.size <- theme(text = element_text(size = standard.textsize))
original_country_count <- function(envir) {
 evalq({
p.country.count.hbar <- country_order_df %>%
    mutate(Geo_Country = factor(Geo_Country, levels = rev(country_levels))) %>%
    ggplot(aes(y = Geo_Country, x = Count, fill = Continent)) +
    geom_col(width = 0.58, alpha = 0.86) +
    geom_text(aes(label = Count), hjust = -0.18, size = text.size.within) +
    facet_grid(Continent ~ ., scales = "free_y", space = "free_y", labeller = as_labeller(CONTINENT_DISPLAY_LABELS)) +
    scale_fill_manual(values = setNames(Continent.cols$Continent.col, Continent.cols$Continent),
                       labels = CONTINENT_DISPLAY_LABELS[Continent.cols$Continent], na.value = "grey70") +
    scale_y_discrete(labels = function(x) gsub("_", " ", x)) +
    coord_cartesian(xlim = c(0, max(country_order_df$Count, na.rm = TRUE) * 1.22)) +
    labs(title = "Country-level sample counts", y = "Country", x = "Sample count", fill = "Continent") +
    theme_bw(base_size = standard.textsize) + theme.text.size + legend.size + y.facet.text.rotate + legend.strip + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
p.country.count.hbar
 },envir=envir)
}

original_country_time <- function(envir) {
 evalq({
p.country.date.bubbleplot <- meta %>%
    filter(is.finite(year)) %>%
    group_by(Geo_Country, Continent, year) %>%
    summarise(Count.per.country.year = n(), .groups = "drop") %>%
    mutate(Geo_Country = factor(Geo_Country, levels = rev(country_levels))) %>%
    ggplot(aes(as.numeric(year), Geo_Country, colour = Continent, group = Geo_Country)) +
    geom_line(alpha = 0.12, linewidth = 0.18) +
    geom_point(alpha = 0.72, aes(size = Count.per.country.year)) +
    facet_grid(Continent ~ ., scales = "free_y", space = "free_y", labeller = as_labeller(CONTINENT_DISPLAY_LABELS)) +
    scale_color_manual(values = setNames(Continent.cols$Continent.col, Continent.cols$Continent),
                        labels = CONTINENT_DISPLAY_LABELS[Continent.cols$Continent], na.value = "grey70") +
    scale_y_discrete(labels = function(x) gsub("_", " ", x)) +
    scale_size_area(max_size = 5.5, breaks = c(1, 5, 10, 25, 50, 100)) +
    labs(title = "Sampling through time", y = "Country", x = "Sample year", size = "Count") +
    theme_bw(base_size = standard.textsize) + theme.text.size + legend.size + y.facet.text.rotate + y.theme.strip + theme(plot.title = element_text(hjust = 0.5, face = "bold")) + theme(plot.title = element_text(hjust = 0.5, face = "bold")) + guides(colour = "none") + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
p.country.date.bubbleplot
 },envir=envir)
}

original_country_lineage <- function(envir) {
 evalq({
p.country.lineage.proportion <- meta %>%
    group_by(Continent, Geo_Country, TPA_Lineage) %>%
    summarise(count = n(), .groups = "drop") %>%
    mutate(Geo_Country = factor(Geo_Country, levels = rev(country_levels))) %>%
    ggplot(aes(y = Geo_Country, x = count, fill = TPA_Lineage)) +
    geom_col(width = 0.65, position = "fill") +
    facet_grid(Continent ~ ., scales = "free_y", space = "free_y", labeller = as_labeller(CONTINENT_DISPLAY_LABELS)) +
    scale_fill_manual(values = setNames(TPA_Lineage.cols$Lineage.col, TPA_Lineage.cols$Lineage), name = "Lineage") +
    scale_x_continuous(expand = c(0.01, 0.01), breaks = seq(0, 1, 0.5), labels = scales::percent) +
    labs(title = "Lineage proportion", y = "Country", x = "Lineage proportion") +
    theme_bw(base_size = standard.textsize) + theme.text.size + legend.size + y.facet.text.rotate + y.theme.strip +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
p.country.lineage.proportion
 },envir=envir)
}

original_map <- function(envir) {
 evalq({
p.map <- ggplot() +
    geom_polygon(data = world, aes(x = long, y = lat, group = group), fill = "grey94", color = "white", linewidth = 0.12) +
    geom_point(data = map_df, aes(x = long, y = lat, size = Count), shape = 21, fill = "#ec7014", color = "grey30", alpha = 0.72, stroke = 0.20) +
    ggrepel::geom_text_repel(
        data = map_label_df,
        aes(x = long, y = lat, label = label_country),
        size = 2.25,
        color = "grey15",
        min.segment.length = 0,
        segment.size = 0.15,
        segment.linetype = "dashed",
        segment.color = "grey45",
        box.padding = 0.18,
        point.padding = 0.08,
        force = 1.0,
        max.overlaps = Inf
    ) +
    scale_size_area(max_size = 9, name = "Samples") +
    coord_quickmap(xlim = c(-170, 180), ylim = c(-55, 80), expand = FALSE) +
    labs(title = "Geographic distribution of sampled genomes") +
    theme_void(base_size = standard.textsize) + theme(legend.position = "right", plot.title = element_text(face = "bold", hjust = 0.5))
p.map
 },envir=envir)
}

original_variant_count <- function(envir) {
 evalq({
p.resistance.counts <- ggplot(top_mut, aes(y = axis_label, x = n_carrier, fill = gene_category)) +
    geom_col(width = 0.65, alpha = 0.86) +
    geom_text(aes(label = n_carrier), hjust = -0.10, size = 1.75, color = "grey25") +
    facet_grid(gene ~ ., scales = "free_y", space = "free_y") +
    scale_fill_manual(values = c("PBP genes" = "#B2182B",
                                  "Peptidoglycan-pathway gene (lytB)" = "grey55"),
                      name = NULL) +
    coord_cartesian(xlim = c(0, max(top_mut$n_carrier, na.rm = TRUE) * 1.20)) +
    labs(title = "Candidate PBP-associated gene variants",
         subtitle = "All detected PBP-associated variants are shown (matches Panel F).",
         y = NULL, x = "Mutant genomes") +
    theme_bw(base_size = standard.textsize) + theme.text.size +
    theme(axis.text.y = element_text(size = 5.5), legend.position = "bottom",
          legend.text = element_text(size = 6), legend.key.size = unit(0.35, "cm"),
          strip.text.y = element_text(angle = 0, face = "bold", size = 6.5),
          plot.subtitle = element_text(size = 7.2, color = "grey35"))
p.resistance.counts
 },envir=envir)
}

original_variant_time <- function(envir) {
 evalq({
p.pbp.timeline <- ggplot(pbp_long, aes(x = year, y = axis_label, size = n, color = focal_status)) +
    geom_point(alpha = 0.40, na.rm = TRUE) +
    facet_grid(gene ~ ., scales = "free_y", space = "free_y") +
    scale_color_manual(values = c("Focal PBP variants" = "#B2182B", "Other PBP variants" = "grey35"), name = NULL) +
    scale_size_area(max_size = 6.5, name = "Mutant genomes") +
    labs(title = "PBP mutation occurrence through time",
         y = "PBP mutation (overall frequency)", x = "Sample year") +
    theme_bw(base_size = standard.textsize) + theme.text.size +
    theme(axis.text.y = element_text(size = 5.5), legend.position = "bottom",
          strip.text.y = element_text(angle = 0, face = "bold", size = 7))
p.pbp.timeline
 },envir=envir)
}

original_country_correlation <- function(envir) {
 evalq({
p2b <- ggplot(lineage_corr_long, aes(x = SS14_pct, y = variant_pct)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey72", linewidth = 0.30) +
    geom_smooth(aes(color = fig2_group), method = "lm", se = FALSE, linewidth = 0.45, linetype = "solid", show.legend = FALSE) +
    geom_point(aes(size = n, fill = fig2_group), shape = 21, color = "grey25", stroke = 0.25, alpha = 0.82) +
    ggrepel::geom_text_repel(aes(label = gsub("_", " ", Geo_Country)), size = 1.8, color = "grey30", max.overlaps = 10, min.segment.length = 0, segment.size = 0.15, segment.linetype = "dashed", segment.alpha = 0.55, show.legend = FALSE) +
    geom_text(data = lineage_corr_stats, aes(x = 5, y = Inf, label = label),
              inherit.aes = FALSE, hjust = 0, vjust = 1.35, size = 2.8, fontface = "bold", color = "grey20", show.legend = FALSE) +
    facet_wrap(~ variant, nrow = 1, scales = "free_y") +
    scale_fill_manual(values = setNames(Fig2.cols$fig2_col, Fig2.cols$fig2_group), name = "Variant colour group") +
    scale_color_manual(values = setNames(Fig2.cols$fig2_col, Fig2.cols$fig2_group), guide = "none") +
    scale_size_area(max_size = 4, name = "Country n") +
    labs(title = "Variant frequency versus SS14 lineage frequency",
         subtitle = "Each point is a country; the panel shows whether focal-variant frequency tracks the country-level SS14-lineage fraction.",
         x = "SS14 lineage frequency by country (%)",
         y = "Variant frequency by country (%)") +
    theme_bw(base_size = standard.textsize) + theme.text.size + legend.size +
    theme(strip.text = element_text(face = "bold"), legend.position = "bottom")
p2b
 },envir=envir)
}
