# ---- Interactive Sankey (Plotly + htmlwidgets) --------------------------------
# Single-trace approach: one Sankey trace is initialised; Plotly.react() swaps
# link data on threshold / method change.  Node definition is embedded as clean
# JSON (no Plotly-computed positions), so the layout is recalculated fresh each
# time instead of inheriting stale coordinates from the previous render.
#
# Run from project root: source("interactive_sankey/make_sankey.R")
# Required: dplyr, readr, plotly, htmlwidgets, jsonlite

library(dplyr)
library(readr)
library(plotly)
library(htmlwidgets)
library(jsonlite)

# ---- Load data ---------------------------------------------------------------
df_raw <- read_csv(
  "results/predictions/serra_martin_loo_prediction_results.csv",
  show_col_types = FALSE
) %>%
  rename(location = focal_site) %>%
  mutate(interaction_id = paste(higher_level, lower_level, sep = "___"))

add_obs_for_obs <- df_raw %>%
  filter(method == "rpi", ground_truth == 1) %>%
  pull(interaction_id) %>% unique()

add_obs_for_rpi <- df_raw %>%
  filter(method == "obs", ground_truth == 1) %>%
  pull(interaction_id) %>% unique()

# ---- Node definitions (25 nodes, 0-indexed) ----------------------------------
NODE_IDX <- c(
  "L1_TP" = 0L, "L1_FP" = 1L, "L1_FN" = 2L, "L1_TN" = 3L,
  "L2_Observed elsewhere"     = 4L,
  "L2_Not observed elsewhere" = 5L,
  "L3_recurrent"              = 6L,  "L3_possibly_missing"      = 7L,
  "L3_model_elusive"          = 8L,  "L3_locally_absent"        = 9L,
  "L3_locally_unique"         = 10L, "L3_phantom"               = 11L,
  "L3_weak support"           = 12L, "L3_likely_forbidden"      = 13L,
  "L4_recurrent"              = 14L, "L4_possibly_missing"      = 15L,
  "L4_model_elusive"          = 16L, "L4_locally_absent"        = 17L,
  "L4_locally_unique"         = 18L,
  "L4_phantom — Have evidence"          = 19L,
  "L4_phantom — No evidence"            = 20L,
  "L4_weak support — Have evidence"     = 21L,
  "L4_weak support — No evidence"       = 22L,
  "L4_likely_forbidden — Have evidence" = 23L,
  "L4_likely_forbidden — No evidence"   = 24L
)

node_labels <- c(
  "TP", "FP", "FN", "TN",
  "Observed elsewhere", "Not observed elsewhere",
  "Recurrent", "Possibly missing", "Model elusive", "Locally absent",
  "Locally unique", "Phantom", "Weak support", "Likely forbidden",
  "Recurrent", "Possibly missing", "Model elusive", "Locally absent",
  "Locally unique",
  "Phantom — Have evidence",          "Phantom — No evidence",
  "Weak support — Have evidence",     "Weak support — No evidence",
  "Likely forbidden — Have evidence", "Likely forbidden — No evidence"
)

node_colors_hex <- c(
  "#CAE1FF", "#B0C4DE", "#EEA9A9", "#BC8F8F",       # L1 confusion
  "#F4A460", "#CDB5CD",                               # L2 evidence
  "#FF7F50", "#FF7256", "#EE7056", "#CD5B45",         # L3 TP/FP group
  "#C3BAD5", "#AFA2C4", "#9B8BB4", "#8878A4",         # L3 FN/TN group
  "#FF7F50", "#FF7256", "#EE7056", "#CD5B45",         # L4 pass-throughs (TP/FP)
  "#C3BAD5",                                           # L4 locally_unique
  "#7FFFD4", "#AFA2C4",                               # L4 phantom have/no
  "#76EEC6", "#9B8BB4",                               # L4 weak support have/no
  "#66CDAA", "#8878A4"                                 # L4 likely_forbidden have/no
)

# Fixed node positions: x pins each node to its axis column; y sets vertical
# order (0 = top, 1 = bottom) to match the static R Sankey.
# arrangement = "fixed" makes Plotly respect these coordinates exactly,
# which also stabilises the layout when many categories are empty at high thresholds.
evenly_y <- function(n) seq(0.02, 0.98, length.out = n)

# Node order within each axis matches node_labels (top → bottom):
# L1: TP, FP, FN, TN         (4 nodes, indices 0–3)
# L2: Observed, Not observed  (2 nodes, indices 4–5)
# L3: Recurrent … Likely forbidden (8 nodes, indices 6–13)
# L4: Recurrent … LF—No evidence  (11 nodes, indices 14–24)
node_x <- c(rep(0.01, 4), rep(0.34, 2), rep(0.67, 8), rep(0.99, 11))
node_y <- c(evenly_y(4), evenly_y(2), evenly_y(8), evenly_y(11))

node_def_r <- list(
  pad       = 15,
  thickness = 20,
  label     = node_labels,
  color     = node_colors_hex,
  x         = node_x,
  y         = node_y,
  line      = list(color = "white", width = 2)  # <-- outline: change width to adjust thickness
)

hex_to_rgba <- function(hex, alpha = 0.45) {
  sprintf("rgba(%d,%d,%d,%.2f)",
          strtoi(substr(hex, 2L, 3L), 16L),
          strtoi(substr(hex, 4L, 5L), 16L),
          strtoi(substr(hex, 6L, 7L), 16L),
          alpha)
}

# ---- Classification ----------------------------------------------------------
classify_sankey <- function(df, method_name, threshold, add_obs_ids) {
  df %>%
    filter(method == method_name) %>%
    group_by(interaction_id) %>%
    mutate(n_obs_total = sum(ground_truth == 1L)) %>%
    ungroup() %>%
    mutate(
      cls           = as.integer(probability >= threshold),
      is_all_zero   = n_obs_total == 0L,
      is_unique     = n_obs_total == 1L,
      is_shared     = n_obs_total >= 2L,
      obs_elsewhere = (n_obs_total - as.integer(ground_truth == 1L)) >= 1L,
      confusion = case_when(
        ground_truth == 1L & cls == 1L ~ "TP",
        ground_truth == 1L & cls == 0L ~ "FN",
        ground_truth == 0L & cls == 1L ~ "FP",
        TRUE                           ~ "TN"
      ),
      validation    = if_else(obs_elsewhere,
                              "Observed elsewhere", "Not observed elsewhere"),
      link_category = case_when(
        is_unique    & ground_truth == 1L & cls == 1L ~ "locally_unique",
        is_unique    & ground_truth == 1L & cls == 0L ~ "weak support",
        is_shared    & ground_truth == 1L & cls == 1L ~ "recurrent",
        is_shared    & ground_truth == 1L & cls == 0L ~ "model_elusive",
        is_all_zero  & ground_truth == 0L & cls == 0L ~ "likely_forbidden",
        is_all_zero  & ground_truth == 0L & cls == 1L ~ "phantom",
        obs_elsewhere & ground_truth == 0L & cls == 0L ~ "locally_absent",
        obs_elsewhere & ground_truth == 0L & cls == 1L ~ "possibly_missing",
        TRUE ~ "unclassified"
      ),
      add_obs = interaction_id %in% add_obs_ids,
      L4 = case_when(
        link_category %in% c("phantom", "weak support", "likely_forbidden") & add_obs ~
          paste0(link_category, " — Have evidence"),
        link_category %in% c("phantom", "weak support", "likely_forbidden") & !add_obs ~
          paste0(link_category, " — No evidence"),
        TRUE ~ link_category
      )
    ) %>%
    filter(link_category != "unclassified")
}

build_links <- function(df_cl) {
  flows <- df_cl %>%
    count(confusion, validation, link_category, L4, name = "value")

  lnks <- bind_rows(
    flows %>% transmute(
      source = NODE_IDX[paste0("L1_", confusion)],
      target = NODE_IDX[paste0("L2_", validation)],
      value
    ),
    flows %>% transmute(
      source = NODE_IDX[paste0("L2_", validation)],
      target = NODE_IDX[paste0("L3_", link_category)],
      value
    ),
    flows %>% transmute(
      source = NODE_IDX[paste0("L3_", link_category)],
      target = NODE_IDX[paste0("L4_", L4)],
      value
    )
  ) %>%
    filter(!is.na(source), !is.na(target)) %>%
    group_by(source, target) %>%
    summarise(value = sum(value), .groups = "drop")

  lnks$color <- hex_to_rgba(node_colors_hex[lnks$source + 1L])
  lnks
}

# ---- Pre-compute all link datasets -------------------------------------------
# Threshold range 0.50 – 1.00 in steps of 0.05 (11 thresholds × 2 methods)
thresholds <- seq(0.50, 1.00, by = 0.05)

cat("Pre-computing link datasets...\n")
link_store <- list()

for (meth in c("obs", "rpi")) {
  add_ids <- if (meth == "obs") add_obs_for_obs else add_obs_for_rpi
  for (thresh in thresholds) {
    df_cl <- classify_sankey(df_raw, meth, thresh, add_ids)
    lnks  <- build_links(df_cl)
    key   <- sprintf("%s_%.2f", meth, thresh)
    link_store[[key]] <- list(
      source = as.integer(lnks$source),
      target = as.integer(lnks$target),
      value  = as.numeric(lnks$value),
      color  = lnks$color
    )
    cat(sprintf("  %s  thresh=%.2f  links=%d\n", meth, thresh, nrow(lnks)))
  }
}

link_json <- jsonlite::toJSON(link_store, auto_unbox = FALSE)
node_json <- jsonlite::toJSON(node_def_r, auto_unbox = TRUE)

# ---- Initial trace (obs, threshold = 0.50) -----------------------------------
init <- link_store[["obs_0.50"]]

fig <- plot_ly(
  type        = "sankey",
  arrangement = "fixed",
  node = node_def_r,
  link = list(
    source = init$source,
    target = init$target,
    value  = init$value,
    color  = init$color
  )
) %>%
  layout(
    title = list(
      text = paste0(
        "Pollinator Interaction Link Classification — Interactive Sankey",
        "<br><sup>Click any node to highlight its upstream path  ·  ",
        "click again to clear  ·  hover for counts</sup>"
      ),
      font = list(size = 13, color = "#444"),
      x = 0.01, xanchor = "left"
    ),
    font          = list(family = "Segoe UI, system-ui, sans-serif", size = 16),
    paper_bgcolor = "#f7f7f5",
    plot_bgcolor  = "#f7f7f5",
    height        = 600,
    margin        = list(l = 8, r = 8, t = 65, b = 8)
  )

# ---- onRender: controls + upstream click highlighting -----------------------
js <- paste0("
function(el, data) {

  // Pre-baked link datasets keyed by 'method_threshold' (e.g. 'obs_0.50')
  var linkData = ", link_json, ";
  // Clean node definition without Plotly-computed positions
  var nodeData = ", node_json, ";

  // Deep-copy the initial layout so react always restores the same settings
  var origLayout = JSON.parse(JSON.stringify(el.layout));

  var thresholds    = ['0.50','0.55','0.60','0.65','0.70','0.75',
                       '0.80','0.85','0.90','0.95','1.00'];
  var currentMethod = 'obs';
  var currentThreshI = 0;          // index 0 -> 0.50
  var currentLinks  = linkData['obs_0.50'];
  var highlightedNode = null;

  // ── Helpers ────────────────────────────────────────────────────────────────
  function setAlpha(rgba, a) {
    // Replace the alpha in 'rgba(r,g,b,OLD)' with NEW
    return rgba.replace(/,[\\d.]+\\)$/, ',' + a + ')');
  }

  function getUpstreamLinkSet(nodeIdx) {
    var src = currentLinks.source;
    var tgt = currentLinks.target;
    var upLinks = new Set();
    var upNodes = new Set([nodeIdx]);

    (function trace(n) {
      for (var i = 0; i < src.length; i++) {
        if (tgt[i] === n && !upLinks.has(i)) {
          upLinks.add(i);
          if (!upNodes.has(src[i])) {
            upNodes.add(src[i]);
            trace(src[i]);
          }
        }
      }
    })(nodeIdx);

    return upLinks;
  }

  function applyHighlight(nodeIdx) {
    var upLinks   = getUpstreamLinkSet(nodeIdx);
    var newColors = currentLinks.color.map(function(c, i) {
      return upLinks.has(i) ? setAlpha(c, 0.85) : 'rgba(210,210,210,0.18)';
    });
    Plotly.restyle(el, {'link.color': [newColors]});
  }

  function clearHighlight() {
    Plotly.restyle(el, {'link.color': [currentLinks.color.slice()]});
    highlightedNode = null;
  }

  // ── Data update (called on method/threshold change) ────────────────────────
  function update() {
    var key = currentMethod + '_' + thresholds[currentThreshI];
    var ld  = linkData[key];
    if (!ld) { console.warn('No data for key:', key); return; }
    currentLinks    = ld;
    highlightedNode = null;

    // Use the clean node definition — Plotly recalculates layout from scratch
    Plotly.react(el, [{
      type: 'sankey',
      arrangement: 'fixed',
      node: nodeData,
      link: { source: ld.source, target: ld.target,
              value:  ld.value,  color:  ld.color }
    }], origLayout);
  }

  // ── Click handler: upstream highlight only ---------------------------------
  el.on('plotly_click', function(ev) {
    if (!ev || !ev.points || !ev.points.length) return;
    var pt = ev.points[0];

    // Sankey link clicks expose a 'source' property; node clicks do not
    if (pt.source !== undefined) { clearHighlight(); return; }

    var nodeIdx = pt.pointNumber;
    if (nodeIdx == null) { clearHighlight(); return; }

    if (highlightedNode === nodeIdx) {
      clearHighlight();            // click same node again -> clear
    } else {
      highlightedNode = nodeIdx;
      applyHighlight(nodeIdx);
    }
  });

  // ── Inject controls above the chart ───────────────────────────────────────
  var ff = 'Segoe UI,system-ui,sans-serif';

  function btn(id, label, active) {
    return '<button id=\"sk-' + id + '\" style=\"' +
      'padding:5px 12px;border-radius:4px;cursor:pointer;font-size:12px;' +
      'font-family:' + ff + ';transition:all .12s;' +
      'border:1.5px solid ' + (active ? '#3b6ea5' : '#ccc') + ';' +
      'background:' + (active ? '#3b6ea5' : '#fff') + ';' +
      'color:' + (active ? '#fff' : '#666') + ';' +
      'font-weight:' + (active ? '600' : '400') + '\">' + label + '</button>';
  }

  var ctrl = document.createElement('div');
  ctrl.style.cssText = [
    'padding:8px 12px', 'background:#fff', 'border-radius:6px',
    'box-shadow:0 1px 4px rgba(0,0,0,.1)', 'margin-bottom:8px',
    'display:flex', 'align-items:center', 'gap:18px', 'flex-wrap:wrap',
    'font-family:' + ff
  ].join(';');

  ctrl.innerHTML =
    '<div style=\"display:flex;align-items:center;gap:5px\">' +
      '<span style=\"font-size:10px;font-weight:700;color:#888;' +
        'text-transform:uppercase;letter-spacing:.06em\">Method</span>' +
      btn('obs', 'Direct observation', true) +
      btn('rpi', 'Camera (rpi)', false) +
    '</div>' +
    '<div style=\"display:flex;align-items:center;gap:8px\">' +
      '<span style=\"font-size:10px;font-weight:700;color:#888;' +
        'text-transform:uppercase;letter-spacing:.06em\">Threshold</span>' +
      '<input type=\"range\" id=\"sk-thr\" min=\"0\" max=\"10\" step=\"1\" value=\"0\"' +
        ' style=\"width:160px;cursor:pointer;accent-color:#3b6ea5\">' +
      '<span id=\"sk-tv\" style=\"font-weight:700;color:#3b6ea5;' +
        'min-width:3.5ch;font-variant-numeric:tabular-nums\">0.50</span>' +
    '</div>';

  el.parentNode.insertBefore(ctrl, el);

  function syncMethodBtns(m) {
    ['obs','rpi'].forEach(function(id) {
      var b = document.getElementById('sk-' + id);
      var on = (id === m);
      b.style.background  = on ? '#3b6ea5' : '#fff';
      b.style.color       = on ? '#fff'    : '#666';
      b.style.borderColor = on ? '#3b6ea5' : '#ccc';
      b.style.fontWeight  = on ? '600'     : '400';
    });
  }

  document.getElementById('sk-obs').onclick = function() {
    currentMethod = 'obs'; syncMethodBtns('obs'); update();
  };
  document.getElementById('sk-rpi').onclick = function() {
    currentMethod = 'rpi'; syncMethodBtns('rpi'); update();
  };

  document.getElementById('sk-thr').oninput = function() {
    currentThreshI = parseInt(this.value);
    document.getElementById('sk-tv').textContent = thresholds[currentThreshI];
    update();
  };
}
")

fig <- htmlwidgets::onRender(fig, js)

# ---- Save and open -----------------------------------------------------------
out_path <- "interactive_sankey/sankey_interactive.html"

htmlwidgets::saveWidget(
  widget        = fig,
  file          = normalizePath(out_path, mustWork = FALSE),
  selfcontained = TRUE,
  title         = "Pollinator Link Classification — Interactive Sankey"
)

cat(sprintf("\nSaved: %s  (%.0f KB)\n", out_path, file.size(out_path) / 1024))
browseURL(normalizePath(out_path))
