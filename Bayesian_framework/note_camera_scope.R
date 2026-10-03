# ---- Camera corroboration: pooled across sites, or matched to the site? ----
# PURPOSE   The deterministic analysis and the Bayesian one flag a link as
#           camera-corroborated by different rules. This note quantifies the
#           difference and writes itself to html.
# WHAT      Every number below is computed here from the prediction CSV and the
#           posteriors, so the note cannot drift from the data.
# OUTPUT    Bayesian_framework/note_camera_scope.html

library(tidyverse)

OUT_HTML <- "Bayesian_framework/note_camera_scope.html"

d <- read_csv("results/predictions/serra_marin_loo_prediction_results.csv",
              show_col_types = FALSE) %>%
  mutate(pair = paste(higher_level, lower_level, sep = "___"))

cam <- filter(d, method == "rpi")

# the two definitions
pooled_pairs <- cam %>% filter(ground_truth == 1) %>% pull(pair) %>% unique()
site_pairs   <- cam %>% filter(ground_truth == 1) %>%
  transmute(key = paste(pair, focal_site)) %>% pull() %>% unique()
watched      <- cam %>% transmute(key = paste(pair, focal_site)) %>% pull() %>% unique()

# the deterministic category of every link, from the Bayesian outputs (same rule)
x <- read_csv("results/bayesian_empirical/posteriors_bayesian_example.csv",
              show_col_types = FALSE) %>%
  filter(run == "B-uniform") %>%
  transmute(focal_site, det_category,
            pair = paste(pollinator, plant, sep = "___")) %>%
  mutate(key       = paste(pair, focal_site),
         ev_pooled = pair %in% pooled_pairs,
         ev_site   = key %in% site_pairs,
         watched   = key %in% watched)

by_cat <- x %>%
  group_by(category = det_category) %>%
  summarise(links = n(),
            pooled = sum(ev_pooled), pct_pooled = 100 * mean(ev_pooled),
            site = sum(ev_site), pct_site = 100 * mean(ev_site),
            in_grid = sum(watched),
            pct_watched = 100 * sum(ev_site) / pmax(sum(watched), 1),
            .groups = "drop") %>%
  arrange(desc(links))

by_site <- x %>%
  group_by(site = focal_site) %>%
  summarise(links = n(), watched = sum(watched), pct_watched = 100 * mean(watched),
            pooled = sum(ev_pooled), site_matched = sum(ev_site), .groups = "drop")

RHO <- cor(by_cat$pct_pooled, by_cat$pct_site, method = "spearman")
TOT_POOLED <- sum(x$ev_pooled); TOT_SITE <- sum(x$ev_site); TOT <- nrow(x)
NEVER_WATCHED <- sum(!x$watched)

rich <- x %>% filter(focal_site == "Rocky coastal 1",
                     det_category %in% c("phantom", "possibly forbidden")) %>%
  group_by(category = det_category) %>%
  summarise(links = n(), pooled = sum(ev_pooled), site = sum(ev_site),
            watched = sum(watched), .groups = "drop")

# --- html helpers -----------------------------------------------------------
pc  <- function(x, d = 1) sprintf(paste0("%.", d, "f%%"), x)
tbl <- function(df, align = NULL, bold_rows = integer()) {
  if (is.null(align)) align <- ifelse(sapply(df, is.numeric), "right", "left")
  head <- paste0("<th style='text-align:", align, "'>", names(df), "</th>", collapse = "")
  rows <- vapply(seq_len(nrow(df)), function(i) {
    cells <- paste0("<td style='text-align:", align, "'>",
                    vapply(df[i, ], function(v) as.character(v), ""), "</td>", collapse = "")
    sprintf("<tr%s>%s</tr>", if (i %in% bold_rows) " class='hl'" else "", cells)
  }, "")
  paste0("<table><thead><tr>", head, "</tr></thead><tbody>",
         paste(rows, collapse = ""), "</tbody></table>")
}

cat_tbl <- by_cat %>%
  transmute(Category = category, Links = links,
            `Pooled across sites` = sprintf("%d (%s)", pooled, pc(pct_pooled)),
            `Matched to the site` = sprintf("%d (%s)", site, pc(pct_site)),
            `Of links watched there` = sprintf("%d / %d (%s)", site, in_grid, pc(pct_watched)))

site_tbl <- by_site %>%
  transmute(Site = site, Links = links,
            `Watched by cameras` = sprintf("%d (%s)", watched, pc(pct_watched)),
            `Pooled flags` = pooled, `Site-matched flags` = site_matched)

rank_tbl <- tibble(
  `By pooled rate` = by_cat$category[order(-by_cat$pct_pooled)],
  `By site-matched rate` = by_cat$category[order(-by_cat$pct_site)])

rich_tbl <- rich %>%
  transmute(Category = category, Links = links, `Pooled` = pooled,
            `Site-matched` = site, `Watched there` = watched)

html <- sprintf('
<!-- written by Bayesian_framework/note_camera_scope.R -->
<title>Camera corroboration: pooled or site-matched?</title>
<style>
 :root { --accent: #7b7bb8; --accent-dark: #5b5b93; --accent-pale: #eeeef8; }
 body { font-family: "Segoe UI", Helvetica, Arial, sans-serif; color: #24242e;
        max-width: 54em; margin: 2.5em auto; padding: 0 1.4em; line-height: 1.55; }
 h1 { font-size: 1.6em; color: var(--accent-dark); margin-bottom: 0.1em; }
 h2 { font-size: 1.15em; color: var(--accent-dark); margin-top: 1.8em; }
 .date { color: #777; font-size: 0.9em; margin-top: 0; }
 table { border-collapse: collapse; margin: 1em 0; font-size: 0.93em; }
 th, td { padding: 0.35em 0.8em; border-bottom: 1px solid #e4e4ef; }
 th { background: var(--accent-pale); font-weight: 600; }
 tr.hl td { font-weight: 600; }
 code { background: #f4f4f9; padding: 0.1em 0.3em; border-radius: 3px; font-size: 0.92em; }
 .box { background: var(--accent-pale); border-left: 5px solid var(--accent);
        padding: 0.8em 1.2em; margin: 1.4em 0; }
 .method { background: #f7f7fb; border-left: 3px solid #c6c6de;
           padding: 0.6em 1em; margin: 1em 0; font-size: 0.95em; }
</style>

<h1>Camera corroboration: pooled across sites, or matched to the site?</h1>
<p class="date">Written by <code>Bayesian_framework/note_camera_scope.R</code> on %s. Every number is computed from the prediction CSV and the posteriors.</p>

<h2>The two rules</h2>
<p>The deterministic analysis (<code>serra_marin_link_classification_clean.R</code>, and the Sankey app that reads it) builds <code>interaction_id</code> as <code>plant___pollinator</code>, with no site in it. A link is flagged <em>have evidence</em> when the cameras recorded that species pair <strong>anywhere in the six sites</strong>.</p>
<p>The Bayesian analysis (<code>bayesian_empirical.R</code>) joins the camera records by <code>(plant, pollinator, focal_site)</code>, so corroboration means the cameras recorded the pair <strong>at that site</strong>, and a pair outside that site\'s camera grid is marked as not watched rather than as not recorded.</p>
<p>Across the %s candidate links, the pooled rule flags <strong>%d</strong> and the site-matched rule <strong>%d</strong>. So <strong>%s of the current flags</strong> come from a camera record at a different site.</p>

<h2>What changes, per category</h2>
%s
<p>The last column is the honest denominator under the site-matched rule: corroboration among the links the cameras could actually have recorded at that site.</p>

<div class="box">
<strong>The two headline numbers.</strong> "Cameras confirmed 8.5%% of forbidden and 16.3%% of phantom links" are pooled values. Matched to the site they become <strong>%s</strong> and <strong>%s</strong>, or <strong>%s</strong> and <strong>%s</strong> among the links the cameras watched there.
</div>

<h2>What does not change</h2>
<p>The ranking of categories by corroboration survives the switch (Spearman %.2f):</p>
%s
<p>Categories asserting local realisation stay far better corroborated than those asserting absence, which is the qualitative claim the paper makes.</p>

<h2>The figure site</h2>
<p>At Rocky coastal 1, the site all the maps and Figure 5 use:</p>
%s
<p>So the map changes little: a handful of cells lose their corroborated tint.</p>

<h2>Which rule is right</h2>
<p>It depends on the claim, and the two claims differ.</p>
<ul>
<li><strong>Site-level categories need site-matched evidence.</strong> Locally unique, locally absent, recurrent and model-elusive all assert something about the focal site. A camera record at another site does not corroborate "this link is realised here". The overstatement is largest for locally absent (%s pooled against %s matched) and model-elusive (%s against %s).</li>
<li><strong>Feasibility claims can use the pooled record.</strong> For phantom and possibly forbidden the question is whether the interaction can occur at all, so a record anywhere does bear on it. That reading is legitimate as long as it is called regional feasibility rather than corroboration of the local label.</li>
</ul>

<div class="method">
<strong>A caveat that sharpens under the site-matched rule.</strong> %d of the %s links were never watched at their own site, because one partner is missing from that site\'s camera data. Under the pooled rule this is invisible. Under the site-matched rule, "no evidence" mixes "watched and not recorded" with "never watched", which is why rates among watched links are the ones to report. Coverage also varies between sites:
%s
</div>

<h2>If we switch</h2>
<ul>
<li>Put the site into <code>interaction_id</code>, or build the cross-method sets per site, in <code>serra_marin_link_classification_clean.R</code> (<code>add_obs_for_obs</code> and <code>add_obs_for_rpi</code>), and regenerate the Sankey payload.</li>
<li>Update the case-study numbers in the manuscript and in <code>CLAUDE.md</code>.</li>
<li>Report corroboration among watched links, and keep one pooled sentence for the feasibility reading of phantom and possibly forbidden.</li>
<li>The Bayesian analysis needs no change: it is already site-matched, and so is its camera test.</li>
</ul>
',
  format(Sys.Date(), "%d %B %Y"),   # not a sprintf template, so no escaping
  format(TOT, big.mark = ","), TOT_POOLED, TOT_SITE,
  pc(100 * (TOT_POOLED - TOT_SITE) / TOT_POOLED, 0),
  tbl(cat_tbl, bold_rows = which(cat_tbl$Category %in% c("possibly forbidden", "phantom"))),
  pc(by_cat$pct_site[by_cat$category == "possibly forbidden"]),
  pc(by_cat$pct_site[by_cat$category == "phantom"]),
  pc(by_cat$pct_watched[by_cat$category == "possibly forbidden"]),
  pc(by_cat$pct_watched[by_cat$category == "phantom"]),
  RHO, tbl(rank_tbl), tbl(rich_tbl),
  pc(by_cat$pct_pooled[by_cat$category == "locally absent"]),
  pc(by_cat$pct_site[by_cat$category == "locally absent"]),
  pc(by_cat$pct_pooled[by_cat$category == "model-elusive"]),
  pc(by_cat$pct_site[by_cat$category == "model-elusive"]),
  NEVER_WATCHED, format(TOT, big.mark = ","), tbl(site_tbl))

writeLines(html, OUT_HTML)
cat(sprintf("written: %s\n", OUT_HTML))
cat(sprintf("  %d links, pooled %d, site-matched %d, never watched %d, Spearman %.2f\n",
            TOT, TOT_POOLED, TOT_SITE, NEVER_WATCHED, RHO))
