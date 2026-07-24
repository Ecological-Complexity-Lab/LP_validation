# ---- Interactive Sankey app: Serra-Martin link classification ----
#
# Builds a single self-contained HTML file with an interactive version of the
# 4-axis Sankey from serra_martin_link_classification_clean.R:
#
#   Within-network evaluation -> Contextual evidence -> Link category
#     -> Additional method (have / no evidence)
#
# Controls
#   - Sampling method toggle (obs / rpi); the "additional method" is the other one
#   - Threshold slider for binarising predicted probabilities
#   - Split mode: only the ambiguous categories (phantom, weak support,
#     possibly forbidden) split at axis 4, or all eight categories
#
# Interaction
#   - Click any node or flow: everything outside the selected subgraph fades
#     but stays visible, the full upstream path is highlighted, and a readout
#     reports n, % of total, and % of the parent (upstream) category.
#
# Threshold range: probabilities are floored at 0.5 by the pre-sigmoid zero
# floor, and above ~0.89 categories start emptying out. The slider is therefore
# restricted to [0.50, 0.88], where all eight link classes stay populated for
# both methods (verified by sweep; see THRESH_MIN/THRESH_MAX below).
#
# Style follows possibly-missing-explorer: same CSS variables, card layout,
# sticky control bar, info modals, tabular-numeric readouts.
#
# Output: results/figures/interactive_sankey_explorer.html

# ---- 1. Libraries ----
required_cran_packages <- c("dplyr", "readr", "tidyr", "purrr", "jsonlite")

missing_cran <- setdiff(required_cran_packages, rownames(installed.packages()))
if (length(missing_cran) > 0) {
  message("Installing missing packages: ", paste(missing_cran, collapse = ", "))
  install.packages(missing_cran)
}

library(dplyr)
library(readr)
library(tidyr)
library(purrr)
library(jsonlite)

# ---- 2. Parameters ----

INPUT_CSV  <- "results/predictions/serra_martin_loo_prediction_results.csv"
OUTPUT_DIR <- "results/interactive_sankey"
OUTPUT_HTML <- file.path(OUTPUT_DIR, "interactive_sankey_explorer.html")

# Threshold grid precomputed into the app. Restricted to the range where every
# link class stays non-empty for both methods (see header note).
THRESH_MIN  <- 0.50
THRESH_MAX  <- 0.88
THRESH_STEP <- 0.01
THRESH_DEFAULT <- 0.70

thresholds <- seq(THRESH_MIN, THRESH_MAX, by = THRESH_STEP)

method_labels <- c(obs = "Direct observation", rpi = "Raspberry Pi camera")

# ---- 3. Aesthetics (mirrors the static script + explorer palette) ----

lav <- c(
  locally_unique     = "#C3BAD5",
  phantom            = "#AFA2C4",
  weak_support       = "#9B8BB4",
  possibly_forbidden = "#8878A4"
)

col_subcats <- c(
  recurrent        = "#F08A5D",
  locally_unique   = "#C3BAD5",
  model_elusive    = "#E2705A",
  weak_support     = "#9B8BB4",
  possibly_missing = "#E8A33D",
  phantom          = "#AFA2C4",
  locally_absent   = "#B0C4DE",
  possibly_forbidden = "#8878A4"
)

col_confusion <- c(
  TP = "#C6D8EA", FP = "#B0C4DE",
  TN = "#BC8F8F", FN = "#D8A7A7"
)

col_validation <- c(
  `Observed elsewhere`     = "#F4A460",
  `Not observed elsewhere` = "#C3BAD5"
)

COL_HAVE_EVIDENCE <- "#66CDAA"   # aquamarine3, as in the static figure

# Display names for the eight link classes
cat_display <- c(
  recurrent          = "Recurrent",
  locally_unique     = "Locally unique",
  possibly_missing   = "Possibly missing",
  phantom            = "Phantom",
  locally_absent     = "Locally absent",
  possibly_forbidden = "Possibly forbidden",
  weak_support       = "Weak support",
  model_elusive      = "Model elusive"
)

# Vertical ordering within each axis (small = top)
rank_confusion <- c(TN = 10, FN = 20, FP = 30, TP = 40)
rank_validation <- c(`Not observed elsewhere` = 10, `Observed elsewhere` = 20)
rank_category <- c(
  possibly_forbidden = 10, weak_support = 20, phantom = 30,
  locally_unique = 40, locally_absent = 50, model_elusive = 60,
  possibly_missing = 70, recurrent = 80
)

# Categories split into have / no evidence when split_all = FALSE
ambiguous_cats <- c("phantom", "weak_support", "possibly_forbidden")

# ---- 4. Data ----

df_all <- read_csv(INPUT_CSV, show_col_types = FALSE) %>%
  rename(location = focal_site) %>%
  mutate(interaction_id = paste(higher_level, lower_level, sep = "___"))

stopifnot(all(c("obs", "rpi") %in% unique(df_all$method)))

cat(sprintf("Loaded %d link-site records across %d sites and %d methods.\n",
            nrow(df_all), n_distinct(df_all$location), n_distinct(df_all$method)))

# Cross-method evidence: interaction IDs observed by the *other* method.
add_obs_ids <- list(
  obs = df_all %>% filter(method == "rpi", ground_truth == 1) %>%
    pull(interaction_id) %>% unique(),
  rpi = df_all %>% filter(method == "obs", ground_truth == 1) %>%
    pull(interaction_id) %>% unique()
)

# ---- 5. Classification ----

# Classifies every link-site record for one method at one threshold, using
# spatial corroboration across the six sites. Identical logic to
# classify_by_location() in the static script, with the binary prediction
# recomputed from the probability at the supplied threshold.
classify_at_threshold <- function(df_method, threshold, evidence_ids) {

  obs_counts <- df_method %>%
    group_by(interaction_id) %>%
    summarise(n_obs_total = sum(ground_truth == 1), .groups = "drop")

  df_method %>%
    left_join(obs_counts, by = "interaction_id") %>%
    mutate(
      original_binary = ground_truth,
      predicted_bin   = as.integer(probability > threshold),

      is_all_zero   = n_obs_total == 0,
      is_unique     = n_obs_total == 1,
      is_shared     = n_obs_total >= 2,
      obs_elsewhere = (n_obs_total - as.integer(original_binary == 1)) >= 1,

      validation = if_else(obs_elsewhere, "Observed elsewhere", "Not observed elsewhere"),

      confusion = case_when(
        original_binary == 1 & predicted_bin == 1 ~ "TP",
        original_binary == 1 & predicted_bin == 0 ~ "FN",
        original_binary == 0 & predicted_bin == 1 ~ "FP",
        original_binary == 0 & predicted_bin == 0 ~ "TN"
      ),

      link_category = case_when(
        is_unique   & original_binary == 1 & predicted_bin == 1 ~ "locally_unique",
        is_unique   & original_binary == 1 & predicted_bin == 0 ~ "weak_support",
        is_shared   & original_binary == 1 & predicted_bin == 1 ~ "recurrent",
        is_shared   & original_binary == 1 & predicted_bin == 0 ~ "model_elusive",
        is_all_zero & original_binary == 0 & predicted_bin == 0 ~ "possibly_forbidden",
        is_all_zero & original_binary == 0 & predicted_bin == 1 ~ "phantom",
        obs_elsewhere & original_binary == 0 & predicted_bin == 0 ~ "locally_absent",
        obs_elsewhere & original_binary == 0 & predicted_bin == 1 ~ "possibly_missing",
        TRUE ~ "unclassified"
      ),

      add_obs = interaction_id %in% evidence_ids
    ) %>%
    filter(link_category != "unclassified", !is.na(confusion))
}

# Aggregates one classified data frame into the 4-column path table the app
# consumes: one row per unique (L1, L2, L3, L4) path with its count.
build_paths <- function(df_cat, split_all) {
  df_cat %>%
    mutate(
      L1 = confusion,
      L2 = validation,
      L3 = link_category,
      L4 = if (split_all) {
        paste0(link_category, if_else(add_obs, "::have", "::none"))
      } else {
        if_else(
          link_category %in% ambiguous_cats,
          paste0(link_category, if_else(add_obs, "::have", "::none")),
          link_category
        )
      }
    ) %>%
    count(L1, L2, L3, L4, name = "value") %>%
    arrange(L1, L2, L3, L4)
}

# ---- 6. Precompute every (method x threshold x split mode) combination ----

cat("Precomputing classifications across", length(thresholds), "thresholds x 2 methods x 2 split modes...\n")

payload_frames <- list()

for (m in c("obs", "rpi")) {
  df_method <- df_all %>% filter(method == m)

  for (thr in thresholds) {
    df_cat <- classify_at_threshold(df_method, thr, add_obs_ids[[m]])

    for (sa in c(FALSE, TRUE)) {
      key <- paste(m, sprintf("%.2f", thr), if (sa) "all" else "amb", sep = "|")
      payload_frames[[key]] <- build_paths(df_cat, split_all = sa)
    }
  }
  cat(sprintf("  %s done\n", m))
}

# Sanity check: no empty link class anywhere in the precomputed grid.
empty_check <- imap_dfr(payload_frames, function(tbl, key) {
  parts <- strsplit(key, "|", fixed = TRUE)[[1]]
  missing_cats <- setdiff(names(cat_display), unique(tbl$L3))
  if (length(missing_cats) == 0) return(NULL)
  tibble(method = parts[1], threshold = parts[2], mode = parts[3],
         missing = paste(missing_cats, collapse = ", "))
})

if (nrow(empty_check) > 0) {
  warning("Some link classes are empty within the threshold range:")
  print(empty_check)
} else {
  cat(sprintf("Check passed: all 8 link classes populated across [%.2f, %.2f].\n",
              THRESH_MIN, THRESH_MAX))
}

# ---- 7. Assemble the JSON payload ----

# Node metadata: axis, display label, colour, vertical rank.
node_meta <- bind_rows(
  tibble(
    id    = names(rank_confusion),
    axis  = 1L,
    label = c(TN = "True negative", FN = "False negative",
              FP = "False positive", TP = "True positive")[names(rank_confusion)],
    short = names(rank_confusion),
    color = unname(col_confusion[names(rank_confusion)]),
    rank  = unname(rank_confusion)
  ),
  tibble(
    id    = names(rank_validation),
    axis  = 2L,
    label = names(rank_validation),
    short = names(rank_validation),
    color = unname(col_validation[names(rank_validation)]),
    rank  = unname(rank_validation)
  ),
  tibble(
    id    = names(rank_category),
    axis  = 3L,
    label = unname(cat_display[names(rank_category)]),
    short = unname(cat_display[names(rank_category)]),
    color = unname(col_subcats[names(rank_category)]),
    rank  = unname(rank_category)
  ),
  # Axis 4: pass-through nodes reuse the axis-3 identity and colour
  tibble(
    id    = names(rank_category),
    axis  = 4L,
    label = unname(cat_display[names(rank_category)]),
    short = unname(cat_display[names(rank_category)]),
    color = unname(col_subcats[names(rank_category)]),
    rank  = unname(rank_category)
  ),
  # Axis 4: split nodes
  expand_grid(cat = names(rank_category), ev = c("none", "have")) %>%
    transmute(
      id    = paste0(cat, "::", ev),
      axis  = 4L,
      label = paste0(cat_display[cat], if_else(ev == "have", " — have evidence", " — no evidence")),
      short = if_else(ev == "have", "Have evidence", "No evidence"),
      color = if_else(ev == "have", COL_HAVE_EVIDENCE,
                      if_else(cat %in% names(lav), lav[cat], col_subcats[cat])),
      rank  = rank_category[cat] + if_else(ev == "have", 2, 1)
    )
) %>%
  mutate(key = paste(axis, id, sep = "@"))

payload <- list(
  meta = list(
    generated   = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    nRecords    = nrow(df_all),
    nSites      = n_distinct(df_all$location),
    sites       = sort(unique(df_all$location)),
    threshMin   = THRESH_MIN,
    threshMax   = THRESH_MAX,
    threshStep  = THRESH_STEP,
    threshDefault = THRESH_DEFAULT,
    thresholds  = sprintf("%.2f", thresholds),
    methodLabels = as.list(method_labels),
    ambiguous   = ambiguous_cats,
    catDisplay  = as.list(cat_display),
    axisLabels  = c("Within-network evaluation", "Contextual evidence",
                    "Link category", "Additional method")
  ),
  nodes = node_meta,
  data  = map(payload_frames, ~ as.list(.x))
)

payload_json <- toJSON(payload, auto_unbox = TRUE, dataframe = "columns", digits = 6)

cat(sprintf("Payload assembled: %.1f KB of JSON.\n", nchar(payload_json) / 1024))

# ---- 8. HTML template ----
# Style deliberately mirrors possibly-missing-explorer: same CSS custom
# properties, card/summary structure, sticky control bar, info modal.

html_template <- '<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Link classification Sankey explorer</title>
<style>
  :root{
    --accent:#E07B39; --accent-d:#B4531A;
    --text:#1f2933; --muted:#5b6472; --grid:#e5e8ec;
    --card:#ffffff; --bg:#f5f6f8; --border:#e2e5ea;
  }
  *{box-sizing:border-box}
  body{margin:0;background:var(--bg);color:var(--text);
    font-family:system-ui,-apple-system,Segoe UI,Roboto,Helvetica,Arial,sans-serif;
    line-height:1.5;-webkit-font-smoothing:antialiased}
  .wrap{max-width:1180px;margin:0 auto;padding:24px 20px 60px}
  h1{font-size:22px;margin:0 0 4px;font-weight:700}
  .sub{color:var(--muted);font-size:14px;margin:0 0 20px;max-width:900px}
  .sub code{background:#eef0f3;padding:1px 5px;border-radius:4px;font-size:13px}
  .card{background:var(--card);border:1px solid var(--border);border-radius:12px;
    padding:16px 18px;margin-bottom:16px;box-shadow:0 1px 2px rgba(20,30,50,.04);position:relative}
  .card h2{font-size:14px;margin:0 0 12px;font-weight:600;letter-spacing:.01em;padding-right:32px}
  .hint{font-size:12px;color:var(--muted);margin:0 0 8px}
  .muted{color:var(--muted);font-weight:400}
  /* controls */
  .ctrlrow{display:grid;grid-template-columns:auto 1fr auto;gap:22px;align-items:end}
  @media(max-width:860px){.ctrlrow{grid-template-columns:1fr;gap:16px}}
  .ctrl label{display:flex;justify-content:space-between;font-size:13px;font-weight:600;margin-bottom:6px;gap:12px}
  .ctrl .val{color:var(--accent-d);font-variant-numeric:tabular-nums}
  .ctrl .desc{font-size:12px;color:var(--muted);margin-top:4px;font-weight:400}
  input[type=range]{width:100%;accent-color:var(--accent);height:22px}
  .seg-wrap{display:inline-flex;border:1px solid var(--border);border-radius:9px;overflow:hidden;background:#fff}
  .seg{border:none;background:#fff;padding:8px 15px;font-size:12.5px;color:var(--muted);cursor:pointer;
    font-family:inherit;transition:background .12s,color .12s}
  .seg+.seg{border-left:1px solid var(--border)}
  .seg:hover{background:#f0f2f5}
  .seg.active{background:var(--accent);color:#fff;font-weight:600}
  .seg.active:hover{background:var(--accent-d)}
  @media(max-width:520px){.seg-wrap{display:flex;width:100%}.seg{flex:1}}
  .topbar{position:sticky;top:0;z-index:30;background:var(--bg);padding:14px 0 8px}
  .topbar h1{margin:0}
  .card.ctrlbar{position:sticky;top:56px;z-index:20}
  .ctrlbar.stuck{box-shadow:0 4px 16px rgba(20,30,50,.12)}
  .reset{font-size:12px;background:#eef0f3;border:1px solid var(--border);
    border-radius:6px;padding:6px 12px;cursor:pointer;color:var(--muted);font-family:inherit}
  .reset:hover{background:#e4e7ec}
  /* sankey */
  #sankey{width:100%;height:auto;display:block;user-select:none}
  #sankey .flow{cursor:pointer;transition:opacity .15s}
  #sankey .node{cursor:pointer;transition:opacity .15s}
  #sankey .axlab{font-size:12px;font-weight:700;fill:var(--text)}
  #sankey .nlab{font-size:11px;fill:#2c3542;pointer-events:none}
  #sankey .nlab .n{font-weight:400;fill:var(--muted)}
  .faded{opacity:.13}
  .legend{display:flex;gap:14px;flex-wrap:wrap;font-size:12px;color:var(--muted);margin-top:12px}
  .legend span{display:inline-flex;align-items:center;gap:6px}
  .sw{width:12px;height:12px;border-radius:3px;display:inline-block}
  /* readout */
  .stat{display:flex;align-items:flex-start;gap:16px;padding:14px 0;border-top:1px solid var(--border)}
  .stat:first-of-type{border-top:none;padding-top:2px}
  .stat .big{font-size:40px;font-weight:750;color:var(--accent-d);
    font-variant-numeric:tabular-nums;line-height:1;min-width:118px}
  .statmeta{flex:1;min-width:0}
  .statlabel{font-size:13px;color:var(--muted);line-height:1.45}
  .statlabel b{color:var(--text)}
  .pathline{font-size:12.5px;color:var(--muted);margin-top:8px;line-height:1.7}
  .chip{display:inline-block;background:#f2f4f7;border:1px solid var(--border);border-radius:6px;
    padding:2px 8px;margin-right:4px;color:var(--text);font-size:12px}
  .chip .dot{display:inline-block;width:8px;height:8px;border-radius:2px;margin-right:6px;
    vertical-align:middle;position:relative;top:-1px}
  .arrow{color:var(--muted);margin-right:4px}
  .empty{font-size:13px;color:var(--muted);padding:6px 0 2px}
  table{border-collapse:collapse;width:100%;font-size:12.5px}
  th,td{padding:6px 9px;text-align:right;border-bottom:1px solid var(--border);
    font-variant-numeric:tabular-nums;white-space:nowrap}
  th{font-weight:600;color:var(--muted);font-size:11.5px}
  td.name,th.name{text-align:left}
  tr.hl td{background:#fbeee2;font-weight:600;color:var(--accent-d)}
  /* collapsible */
  details.coll{padding:0;overflow:hidden}
  details.coll>summary{list-style:none;cursor:pointer;padding:14px 18px;font-size:14px;font-weight:600;
    display:flex;align-items:center;gap:8px;padding-right:46px}
  details.coll>summary::-webkit-details-marker{display:none}
  details.coll>summary::before{content:"\\25B8";color:var(--muted);font-size:12px;
    transition:transform .15s;display:inline-block}
  details.coll[open]>summary::before{transform:rotate(90deg)}
  .coll-body{padding:0 18px 16px}
  /* info modal */
  .info{position:absolute;top:13px;right:15px;z-index:6;display:inline-flex;align-items:center;
    justify-content:center;width:21px;height:21px;border-radius:50%;border:none;background:var(--accent);
    color:#fff;font-size:12px;font-weight:700;cursor:pointer;line-height:1;font-style:italic;
    font-family:Georgia,"Times New Roman",serif;padding:0;box-shadow:0 1px 3px rgba(20,30,50,.22)}
  .info:hover{background:var(--accent-d)}
  .modal-back{position:fixed;inset:0;background:rgba(20,28,44,.46);z-index:100;display:none;
    align-items:center;justify-content:center;padding:24px}
  .modal-back.open{display:flex}
  .modal-wrap{position:relative;max-width:640px;width:100%}
  .modal{background:var(--card);border-radius:14px;width:100%;max-height:84vh;overflow:auto;
    box-shadow:0 24px 60px rgba(15,23,42,.32);padding:20px 24px 24px}
  .modal h3{margin:0 0 12px;font-size:17px;padding-right:36px}
  .mbody{font-size:14px;color:var(--text);line-height:1.55}
  .mbody p{margin:0 0 11px}
  .mbody p:last-child{margin-bottom:0}
  .mbody b{color:var(--text)}
  .mbody code{background:#eef0f3;padding:1px 4px;border-radius:4px;font-size:13px}
  .modal-close{position:absolute;top:14px;right:16px;border:none;background:#eef0f3;border-radius:8px;
    width:30px;height:30px;cursor:pointer;font-size:18px;color:var(--muted);line-height:1;z-index:1}
  .modal-close:hover{background:#e2e5ea;color:var(--text)}
</style>
</head>
<body>
<div class="wrap">

  <div class="topbar">
    <h1>Link classification &mdash; <span id="hdrMethod" style="color:var(--accent-d)">direct observation</span></h1>
  </div>
  <p class="sub">Every link&ndash;site record flows left to right: from the <b>within-network confusion matrix</b>,
  through <b>contextual evidence</b> (was the same species pair recorded at any other site?), into the
  <b>eight link classes</b>, and finally to whether the <b>additional sampling method</b> found evidence for it.
  Move the threshold to rebinarise the predicted probabilities, and click any band or block to isolate its path.</p>

  <div class="card ctrlbar" id="ctrlbar">
    <button class="info" data-info="controls" aria-label="About the controls">i</button>
    <h2>Controls</h2>
    <div class="ctrlrow">
      <div class="ctrl">
        <label><span>Sampling method</span></label>
        <div class="seg-wrap" id="methodSeg">
          <button class="seg active" data-method="obs">Direct observation</button>
          <button class="seg" data-method="rpi">Raspberry&nbsp;Pi camera</button>
        </div>
        <div class="desc" id="methodDesc">corroborated against the camera records</div>
      </div>
      <div class="ctrl">
        <label><span>Classification threshold</span><span class="val" id="vThr">0.70</span></label>
        <input type="range" id="sThr" min="0" max="1" step="1" value="0">
        <div class="desc" id="thrDesc">predicted probability above which a link is called present</div>
      </div>
      <div class="ctrl">
        <label><span>Axis 4 split</span></label>
        <div class="seg-wrap" id="splitSeg">
          <button class="seg active" data-split="amb">Ambiguous only</button>
          <button class="seg" data-split="all">All classes</button>
        </div>
        <div class="desc">which classes split into have / no evidence</div>
      </div>
    </div>
  </div>

  <div class="card">
    <button class="info" data-info="sankey" aria-label="About the diagram">i</button>
    <h2>Flow of link&ndash;site records <span class="muted" id="totalNote"></span></h2>
    <p class="hint" id="selHint">Click any block or band to isolate it. Click the background to clear.</p>
    <svg id="sankey" viewBox="0 0 1120 620" preserveAspectRatio="xMidYMid meet"
         role="img" aria-label="Interactive Sankey of link classification"></svg>
    <div class="legend" id="legend"></div>
  </div>

  <div class="card" id="readoutCard">
    <button class="info" data-info="readout" aria-label="About the readout">i</button>
    <h2>Selection</h2>
    <div id="readout"></div>
  </div>

  <details class="card coll" id="tableCard">
    <summary>Link class counts at this threshold <span class="muted">&mdash; all eight classes</span></summary>
    <div class="coll-body">
      <table id="catTable"><thead><tr>
        <th class="name">Link class</th><th>n</th><th>% of total</th>
        <th>Have evidence</th><th>% corroborated</th>
      </tr></thead><tbody></tbody></table>
    </div>
  </details>

</div>

<div class="modal-back" id="modalBack">
  <div class="modal-wrap">
    <button class="modal-close" id="modalClose">&times;</button>
    <div class="modal"><h3 id="modalTitle"></h3><div class="mbody" id="modalBody"></div></div>
  </div>
</div>

<script id="payload" type="application/json">__PAYLOAD__</script>
<script>
(function(){
"use strict";

const P = JSON.parse(document.getElementById("payload").textContent);
const META = P.meta;

// ---- rebuild node metadata as a lookup keyed by "axis@id" ----
const NODES = {};
(function(){
  const n = P.nodes;
  for (let i = 0; i < n.id.length; i++) {
    NODES[n.key[i]] = {
      id: n.id[i], axis: n.axis[i], label: n.label[i],
      short: n.short[i], color: n.color[i], rank: n.rank[i]
    };
  }
})();

function nodeInfo(axis, id){
  return NODES[axis + "@" + id] || {id:id, axis:axis, label:id, short:id, color:"#ccc", rank:999};
}

// ---- state ----
const state = {
  method: "obs",
  thrIdx: Math.max(0, META.thresholds.indexOf(META.threshDefault.toFixed(2))),
  split: "amb",
  sel: null          // {type:"node", axis, id} | {type:"flow", axis, from, to}
};

// ---- layout constants ----
const W = 1120, H = 620;
const M = {top: 46, right: 250, bottom: 22, left: 18};
const NODE_W = 11;
const GAP = 7;                 // vertical gap between nodes on an axis
const PLOT_H = H - M.top - M.bottom;
const COL_X = [0, 1, 2, 3].map(i => M.left + i * ((W - M.left - M.right - NODE_W) / 3));

const svg = document.getElementById("sankey");
const SVGNS = "http://www.w3.org/2000/svg";

function el(tag, attrs, parent){
  const e = document.createElementNS(SVGNS, tag);
  for (const k in attrs) e.setAttribute(k, attrs[k]);
  if (parent) parent.appendChild(e);
  return e;
}

function currentPaths(){
  const key = state.method + "|" + META.thresholds[state.thrIdx] + "|" + state.split;
  const d = P.data[key];
  if (!d) return [];
  const out = [];
  for (let i = 0; i < d.value.length; i++){
    out.push({L1:d.L1[i], L2:d.L2[i], L3:d.L3[i], L4:d.L4[i], value:d.value[i]});
  }
  return out;
}

// ---- build the layout from the path table ----
function buildLayout(paths){
  const total = paths.reduce((s,p)=>s+p.value, 0);

  // node totals per axis
  const totals = [{}, {}, {}, {}];
  paths.forEach(p => {
    [p.L1,p.L2,p.L3,p.L4].forEach((id, ax) => {
      totals[ax][id] = (totals[ax][id] || 0) + p.value;
    });
  });

  // vertical scale: leave room for the inter-node gaps on the busiest axis
  const maxNodes = Math.max(...totals.map(t => Object.keys(t).length));
  const usable = PLOT_H - (maxNodes - 1) * GAP;
  const scale = usable / total;

  // position nodes, ordered by rank
  const layout = [{}, {}, {}, {}];
  totals.forEach((t, ax) => {
    const ids = Object.keys(t).sort((a,b)=> nodeInfo(ax+1,a).rank - nodeInfo(ax+1,b).rank);
    const nGaps = ids.length - 1;
    const y0 = M.top + (PLOT_H - (total * scale + nGaps * GAP)) / 2;
    let y = y0;
    ids.forEach(id => {
      const h = Math.max(t[id] * scale, 1.2);
      layout[ax][id] = {y0:y, y1:y+h, h:h, value:t[id], x:COL_X[ax]};
      y += h + GAP;
    });
  });

  // aggregate flows between consecutive axes, ordered so ribbons do not cross
  // more than necessary: sort by the rank of the partner node.
  const flows = [];
  for (let ax = 0; ax < 3; ax++){
    const agg = {};
    paths.forEach(p => {
      const from = [p.L1,p.L2,p.L3][ax], to = [p.L2,p.L3,p.L4][ax];
      const k = from + "\\x01" + to;
      agg[k] = (agg[k] || 0) + p.value;
    });
    const rows = Object.keys(agg).map(k => {
      const [from, to] = k.split("\\x01");
      return {axis:ax, from:from, to:to, value:agg[k]};
    });

    // stack on the source side, ordered by target rank
    const srcCursor = {}, tgtCursor = {};
    rows.sort((a,b)=>{
      const ra = nodeInfo(ax+1,a.from).rank - nodeInfo(ax+1,b.from).rank;
      if (ra !== 0) return ra;
      return nodeInfo(ax+2,a.to).rank - nodeInfo(ax+2,b.to).rank;
    }).forEach(r => {
      const s = layout[ax][r.from], t = layout[ax+1][r.to];
      const h = r.value * scale;
      const sy = srcCursor[r.from] === undefined ? s.y0 : srcCursor[r.from];
      r.sy0 = sy; r.sy1 = sy + h; srcCursor[r.from] = sy + h;
      r.x0 = s.x + NODE_W; r.x1 = t.x;
      r._th = h;
    });

    // stack on the target side, ordered by source rank
    rows.slice().sort((a,b)=>{
      const ra = nodeInfo(ax+2,a.to).rank - nodeInfo(ax+2,b.to).rank;
      if (ra !== 0) return ra;
      return nodeInfo(ax+1,a.from).rank - nodeInfo(ax+1,b.from).rank;
    }).forEach(r => {
      const t = layout[ax+1][r.to];
      const ty = tgtCursor[r.to] === undefined ? t.y0 : tgtCursor[r.to];
      r.ty0 = ty; r.ty1 = ty + r._th; tgtCursor[r.to] = ty + r._th;
    });

    flows.push(...rows);
  }

  return {total:total, totals:totals, layout:layout, flows:flows, paths:paths};
}

// ---- selection logic: which paths survive the current selection ----
function pathMatchesSelection(p, sel){
  if (!sel) return true;
  const ids = [p.L1, p.L2, p.L3, p.L4];
  if (sel.type === "node") return ids[sel.axis] === sel.id;
  return ids[sel.axis] === sel.from && ids[sel.axis+1] === sel.to;
}

function activeSets(paths, sel){
  // returns {nodes: Set("ax@id"), flows: Set("ax@from>to")} for the surviving
  // subgraph, which is the full upstream *and* downstream path of the selection
  const nodes = new Set(), flows = new Set();
  paths.filter(p => pathMatchesSelection(p, sel)).forEach(p => {
    const ids = [p.L1,p.L2,p.L3,p.L4];
    ids.forEach((id, ax) => nodes.add(ax + "@" + id));
    for (let ax = 0; ax < 3; ax++) flows.add(ax + "@" + ids[ax] + ">" + ids[ax+1]);
  });
  return {nodes:nodes, flows:flows};
}

// ---- ribbon path ----
function ribbon(r){
  const cx0 = r.x0 + (r.x1 - r.x0) * 0.42, cx1 = r.x0 + (r.x1 - r.x0) * 0.58;
  return "M" + r.x0 + "," + r.sy0 +
         "C" + cx0 + "," + r.sy0 + " " + cx1 + "," + r.ty0 + " " + r.x1 + "," + r.ty0 +
         "L" + r.x1 + "," + r.ty1 +
         "C" + cx1 + "," + r.ty1 + " " + cx0 + "," + r.sy1 + " " + r.x0 + "," + r.sy1 + "Z";
}

// ---- render ----
let LAY = null;

function render(){
  const paths = currentPaths();
  LAY = buildLayout(paths);
  const act = activeSets(paths, state.sel);
  const hasSel = !!state.sel;

  while (svg.firstChild) svg.removeChild(svg.firstChild);

  // background catcher: clears the selection
  const bg = el("rect", {x:0, y:0, width:W, height:H, fill:"transparent"}, svg);
  bg.addEventListener("click", () => { state.sel = null; render(); renderReadout(); });

  // axis labels
  const axLabels = META.axisLabels.slice();
  axLabels[3] = "Additional method: " + META.methodLabels[state.method === "obs" ? "rpi" : "obs"];
  axLabels.forEach((t, ax) => {
    el("text", {x:COL_X[ax], y:M.top - 22, class:"axlab"}, svg).textContent = t;
  });

  // flows first, so nodes sit on top
  const gF = el("g", {}, svg);
  LAY.flows.forEach(r => {
    const info = nodeInfo(r.axis + 2, r.to);
    const on = act.flows.has(r.axis + "@" + r.from + ">" + r.to);
    const p = el("path", {
      d: ribbon(r), fill: info.color, "fill-opacity": 0.62,
      class: "flow" + (hasSel && !on ? " faded" : "")
    }, gF);
    p.addEventListener("click", ev => {
      ev.stopPropagation();
      const same = state.sel && state.sel.type === "flow" &&
                   state.sel.axis === r.axis && state.sel.from === r.from && state.sel.to === r.to;
      state.sel = same ? null : {type:"flow", axis:r.axis, from:r.from, to:r.to, value:r.value};
      render(); renderReadout();
    });
    el("title", {}, p).textContent =
      nodeInfo(r.axis+1, r.from).label + "  →  " + info.label + " : " + r.value;
  });

  // nodes + labels
  const gN = el("g", {}, svg);
  LAY.layout.forEach((col, ax) => {
    Object.keys(col).forEach(id => {
      const nd = col[id], info = nodeInfo(ax+1, id);
      const on = act.nodes.has(ax + "@" + id);
      const cls = "node" + (hasSel && !on ? " faded" : "");

      const rect = el("rect", {
        x:nd.x, y:nd.y0, width:NODE_W, height:nd.h, rx:2,
        fill:info.color, class:cls
      }, gN);
      rect.addEventListener("click", ev => {
        ev.stopPropagation();
        const same = state.sel && state.sel.type === "node" &&
                     state.sel.axis === ax && state.sel.id === id;
        state.sel = same ? null : {type:"node", axis:ax, id:id};
        render(); renderReadout();
      });
      el("title", {}, rect).textContent = info.label + " : " + nd.value;

      // label only where the band is tall enough, or on the last axis
      if (nd.h >= 9 || ax === 3){
        const t = el("text", {
          x: nd.x + NODE_W + 6, y: nd.y0 + nd.h/2 + 3.6,
          class: "nlab" + (hasSel && !on ? " faded" : "")
        }, gN);
        const lab = (ax === 3 && id.indexOf("::") >= 0) ? info.short : info.label;
        t.textContent = lab + " ";
        const s = el("tspan", {class:"n"}, t);
        s.textContent = nd.value + " (" + pct(nd.value, LAY.total) + ")";
      }
    });
  });

  document.getElementById("totalNote").textContent =
    "— " + LAY.total.toLocaleString() + " records, " + META.nSites + " sites";
  renderLegend();
  renderTable();
}

function pct(a, b){ return b ? (100*a/b).toFixed(1) + "%" : "–"; }

function renderLegend(){
  const box = document.getElementById("legend");
  const order = ["recurrent","locally_unique","possibly_missing","phantom",
                 "locally_absent","possibly_forbidden","weak_support","model_elusive"];
  let h = order.map(c => {
    const i = nodeInfo(3, c);
    return `<span><i class="sw" style="background:${i.color}"></i>${i.label}</span>`;
  }).join("");
  h += `<span><i class="sw" style="background:${nodeInfo(4,"phantom::have").color}"></i>` +
       `Have evidence (other method)</span>`;
  box.innerHTML = h;
}

function renderTable(){
  const paths = LAY.paths;
  const per = {};
  paths.forEach(p => {
    const r = per[p.L3] || (per[p.L3] = {n:0, have:0});
    r.n += p.value;
    if (p.L4.indexOf("::have") >= 0) r.have += p.value;
  });
  const order = ["recurrent","locally_unique","possibly_missing","phantom",
                 "locally_absent","possibly_forbidden","weak_support","model_elusive"];
  const selCat = state.sel && state.sel.type === "node" && state.sel.axis === 2 ? state.sel.id : null;
  const tb = document.querySelector("#catTable tbody");
  tb.innerHTML = order.map(c => {
    const r = per[c] || {n:0, have:0};
    const splitHere = state.split === "all" || META.ambiguous.indexOf(c) >= 0;
    return `<tr class="${c === selCat ? "hl" : ""}">` +
      `<td class="name">${nodeInfo(3,c).label}</td>` +
      `<td>${r.n}</td><td>${pct(r.n, LAY.total)}</td>` +
      `<td>${splitHere ? r.have : "–"}</td>` +
      `<td>${splitHere ? pct(r.have, r.n) : "–"}</td></tr>`;
  }).join("");
}

// ---- readout ----
function chip(axis, id){
  const i = nodeInfo(axis, id);
  const lab = (axis === 4 && id.indexOf("::") >= 0) ? i.label : i.label;
  return `<span class="chip"><i class="dot" style="background:${i.color}"></i>${lab}</span>`;
}

function renderReadout(){
  const box = document.getElementById("readout");
  const sel = state.sel;

  if (!sel){
    box.innerHTML = `<div class="empty">Nothing selected. Click a block or a band to see its ` +
      `counts, its share of the total, and its share of the class it came from.</div>`;
    return;
  }

  const paths = LAY.paths;
  const kept = paths.filter(p => pathMatchesSelection(p, sel));
  const n = kept.reduce((s,p)=>s+p.value, 0);

  // parent = the node one axis upstream that feeds this selection
  let parentAxis, parentN, parentLabel, ownLabel, ownAxis, ownId;
  if (sel.type === "node"){
    ownAxis = sel.axis; ownId = sel.id;
    ownLabel = nodeInfo(sel.axis+1, sel.id).label;
    if (sel.axis === 0){
      parentN = LAY.total; parentLabel = "all records"; parentAxis = null;
    } else {
      // sum over every upstream node feeding into it, i.e. the whole prior axis
      // restricted to the paths that reach this node is the node itself, so the
      // meaningful parent is each distinct upstream node; report the largest
      const up = {};
      kept.forEach(p => {
        const id = [p.L1,p.L2,p.L3,p.L4][sel.axis-1];
        up[id] = (up[id] || 0) + p.value;
      });
      const ids = Object.keys(up);
      if (ids.length === 1){
        parentAxis = sel.axis - 1;
        parentN = LAY.totals[parentAxis][ids[0]];
        parentLabel = nodeInfo(parentAxis+1, ids[0]).label;
      } else {
        parentN = LAY.total; parentLabel = "all records"; parentAxis = null;
      }
    }
  } else {
    ownAxis = sel.axis + 1; ownId = sel.to;
    ownLabel = nodeInfo(sel.axis+1, sel.from).label + " → " + nodeInfo(sel.axis+2, sel.to).label;
    parentAxis = sel.axis;
    parentN = LAY.totals[sel.axis][sel.from];
    parentLabel = nodeInfo(sel.axis+1, sel.from).label;
  }

  // representative full path (the modal upstream route)
  const routes = {};
  kept.forEach(p => {
    const k = [p.L1,p.L2,p.L3,p.L4].join("\\x01");
    routes[k] = (routes[k] || 0) + p.value;
  });
  const topRoute = Object.keys(routes).sort((a,b)=>routes[b]-routes[a])[0].split("\\x01");
  const nRoutes = Object.keys(routes).length;

  let h = "";
  h += `<div class="stat"><div class="big">${n}</div><div class="statmeta">` +
       `<div class="statlabel">link–site records in <b>${ownLabel}</b>, which is ` +
       `<b>${pct(n, LAY.total)}</b> of all ${LAY.total} records</div></div></div>`;
  h += `<div class="stat"><div class="big">${pct(n, parentN)}</div><div class="statmeta">` +
       `<div class="statlabel">of <b>${parentLabel}</b> (${parentN} records) flows here` +
       `</div></div></div>`;

  h += `<div class="pathline"><b>Upstream path` +
       (nRoutes > 1 ? ` (largest of ${nRoutes})` : "") + `:</b><br>`;
  for (let ax = 0; ax < 4; ax++){
    if (ax > 0) h += `<span class="arrow">→</span>`;
    h += chip(ax+1, topRoute[ax]);
  }
  h += "</div>";

  box.innerHTML = h;
}

// ---- controls ----
const sThr = document.getElementById("sThr");
sThr.max = META.thresholds.length - 1;
sThr.value = state.thrIdx;

function syncLabels(){
  document.getElementById("vThr").textContent = META.thresholds[state.thrIdx];
  const other = state.method === "obs" ? "rpi" : "obs";
  document.getElementById("hdrMethod").textContent =
    META.methodLabels[state.method].toLowerCase();
  document.getElementById("methodDesc").textContent =
    "corroborated against " + META.methodLabels[other].toLowerCase() + " records";
  document.getElementById("thrDesc").textContent =
    "predicted probability above which a link is called present (range " +
    META.threshMin.toFixed(2) + "–" + META.threshMax.toFixed(2) +
    ", where all eight classes stay populated)";
}

sThr.addEventListener("input", () => {
  state.thrIdx = +sThr.value;
  syncLabels(); render(); renderReadout();
});

document.getElementById("methodSeg").addEventListener("click", ev => {
  const b = ev.target.closest("button[data-method]"); if (!b) return;
  state.method = b.dataset.method; state.sel = null;
  [...ev.currentTarget.children].forEach(c => c.classList.toggle("active", c === b));
  syncLabels(); render(); renderReadout();
});

document.getElementById("splitSeg").addEventListener("click", ev => {
  const b = ev.target.closest("button[data-split]"); if (!b) return;
  state.split = b.dataset.split; state.sel = null;
  [...ev.currentTarget.children].forEach(c => c.classList.toggle("active", c === b));
  render(); renderReadout();
});

// sticky shadow
const ctrlbar = document.getElementById("ctrlbar");
window.addEventListener("scroll", () => {
  ctrlbar.classList.toggle("stuck", ctrlbar.getBoundingClientRect().top <= 57);
});

// ---- info modals ----
const INFO = {
  controls: {
    t: "Controls",
    b: "<p><b>Sampling method</b> selects which method is being classified. The other method " +
       "supplies the independent evidence shown on axis 4, so switching also swaps the corroborating method.</p>" +
       "<p><b>Threshold</b> rebinarises the leave-one-out predicted probabilities: a record is called " +
       "present when its probability exceeds the threshold. Raising it moves records from the " +
       "predicted-present classes into the predicted-absent ones.</p>" +
       "<p>The slider is limited to <code>" + META.threshMin.toFixed(2) + "–" +
       META.threshMax.toFixed(2) + "</code>. Probabilities are floored at 0.5 by the pre-sigmoid " +
       "zero floor, so nothing changes below that; above the upper limit some classes empty out " +
       "entirely and their proportions become undefined.</p>" +
       "<p><b>Axis 4 split</b> controls which classes are broken into have / no evidence. " +
       "By default only phantom, weak support and possibly forbidden are split, since those are the " +
       "classes where cross-method evidence is diagnostic.</p>"
  },
  sankey: {
    t: "Reading the diagram",
    b: "<p>Each record is one species pair assessed at one site. Bands are scaled by the number of records.</p>" +
       "<p><b>Axis 1</b> is the within-network confusion matrix at the current threshold. " +
       "<b>Axis 2</b> asks whether the same pair was recorded at any <i>other</i> site. " +
       "<b>Axis 3</b> is the resulting link class. <b>Axis 4</b> asks whether the independent method found it.</p>" +
       "<p>Clicking a block or band isolates its subgraph: everything else fades but stays visible, " +
       "so the highlighted route can be read against the whole. Click again, or click the background, to clear.</p>"
  },
  readout: {
    t: "The selection readout",
    b: "<p>The first figure is the number of link–site records in the selection and its share of all records " +
       "at this threshold.</p>" +
       "<p>The second is the share of the <b>upstream</b> category that flows into the selection: for a band, " +
       "the source block it leaves; for a block, the single upstream block feeding it, when there is only one. " +
       "Where several upstream categories converge, the share is reported against the full total instead.</p>" +
       "<p>The path chips trace the full route through all four axes. When a selection covers more than one " +
       "route, the largest is shown and the count of routes is noted.</p>"
  }
};

document.querySelectorAll(".info").forEach(b => {
  b.addEventListener("click", ev => {
    ev.preventDefault(); ev.stopPropagation();
    const k = b.dataset.info; if (!INFO[k]) return;
    document.getElementById("modalTitle").textContent = INFO[k].t;
    document.getElementById("modalBody").innerHTML = INFO[k].b;
    document.getElementById("modalBack").classList.add("open");
  });
});
function closeModal(){ document.getElementById("modalBack").classList.remove("open"); }
document.getElementById("modalClose").addEventListener("click", closeModal);
document.getElementById("modalBack").addEventListener("click", ev => {
  if (ev.target.id === "modalBack") closeModal();
});
document.addEventListener("keydown", ev => { if (ev.key === "Escape") closeModal(); });

// ---- go ----
syncLabels();
render();
renderReadout();

})();
</script>
</body>
</html>
'

# ---- 9. Write the file ----

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

html_out <- sub("__PAYLOAD__", payload_json, html_template, fixed = TRUE)
writeLines(html_out, OUTPUT_HTML, useBytes = TRUE)

cat(sprintf("\nWrote %s (%.1f KB)\n", OUTPUT_HTML, file.size(OUTPUT_HTML) / 1024))
cat("Open it in any browser; the file is self-contained, no server needed.\n")
