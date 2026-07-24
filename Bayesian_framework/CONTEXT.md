# Box 2 — Bayesian reading of the link taxonomy (context for Claude Code)

This folder holds the figure script and the interactive explorer for **Box 2** of the
LP-validation / corroborative-inference manuscript, plus this context file. Read this
first: it captures the mathematical model, the conceptual assumptions, and the design
decisions we made, so you can extend the code without re-deriving everything.

Related existing code in this repo: `code/posterior_accumulation.R` and
`code/posterior_accumulation_no_link_uncertainty.R` (the accumulation analysis on real
data), and `interactive_sankey/` (the same "interactive HTML next to its R" pattern).
The manuscript source is in the Overleaf project (`main_tree.tex`), where this is **Box 2**.
Interactive version is deployed at <https://contextualevidence.ecomplab.com>.

---

## 1. The link taxonomy (8 categories)

Every candidate link is placed by **three binary pieces of evidence**:

| symbol | meaning | 1 / 0 |
|---|---|---|
| `Yhat` (Ŷ) | model **prediction** | predicted present / absent |
| `O_local` | **local** observation | observed locally / not |
| `O_rep` | **contextual** evidence — observed **elsewhere** in within-system replicates | observed elsewhere / not |

`Yhat` and `O_local` form the confusion matrix (TP/FP/TN/FN); `O_rep` splits each cell
in two. The eight categories, with their error-free signatures `(Yhat, O_local, O_rep)`:

| category | Yhat | O_local | O_rep | confusion cell | contextual |
|---|:--:|:--:|:--:|---|---|
| Locally unique      | 1 | 1 | 0 | TP (pred +, local +) | not observed elsewhere |
| Recurrent           | 1 | 1 | 1 | TP | observed elsewhere |
| Phantom             | 1 | 0 | 0 | FP (pred +, local −) | not observed elsewhere |
| Possibly missing    | 1 | 0 | 1 | FP | observed elsewhere |
| Possibly forbidden  | 0 | 0 | 0 | TN (pred −, local −) | not observed elsewhere |
| Locally absent      | 0 | 0 | 1 | TN | observed elsewhere |
| Weak support        | 0 | 1 | 0 | FN (pred −, local +) | not observed elsewhere |
| Model-elusive       | 0 | 1 | 1 | FN | observed elsewhere |

The taxonomy is Fig. 1A of the manuscript. Each category carries a recommended action
(add orthogonal method / increase sampling / evidence is sufficient / within-system
replicates).

---

## 2. Box 2 — the simplest Bayesian case (binary evidence)

The main text assigns each link to one category **deterministically** from the three
evidence bits. Box 2 attaches a **confidence** to that assignment by treating the true
category `C` as unknown and inverting the observed evidence `E = (Ŷ, O_local, O_rep)`
with Bayes:

$$P(C\mid E) = \frac{P(E\mid C)\,P(C)}{P(E)}.$$

Each category `c` corresponds to exactly one error-free evidence pattern `z_c` (its
signature above). Errors enter through the likelihood. With independent error rates
`ε_Y`, `ε_local`, `ε_rep` (one per evidence source):

$$P(E\mid C=c)=\prod_{d\in\{Y,\text{local},\text{rep}\}}
\varepsilon_d^{\,\mathbf 1[E_d\neq z_{c,d}]}\,(1-\varepsilon_d)^{\,\mathbf 1[E_d = z_{c,d}]}.$$

**Assumptions of the simplest case (this is what the figure/explorer implement):**
binary evidence (a link is "observed in ≥1 replicate", yes/no — not a proportion),
**fixed** error rates, **uniform prior** `P(C)=1/8`, and **independent** errors.

Consequences worth knowing:
- With binary evidence + uniform prior the 8 likelihoods sum to 1, so the normaliser is
  1 and `P(C|E) = P(E|C)`.
- Confidence in the *observed* category (evidence = its own signature) is
  `(1−ε_Y)(1−ε_local)(1−ε_rep)`. Under a uniform prior this is **the same for every
  category** — confidence depends only on the error rates, not on which box you're in.
  A **non-uniform prior** breaks that symmetry (rarer categories need stronger evidence).
- **Feasibility.** A link is "truly feasible" if it is realised *somewhere* (locally or
  in replicates), i.e. its true signature has `z_local=1` OR `z_rep=1`. Equivalently
  `P(feasible) = 1 − P(Phantom) − P(Possibly forbidden)` — the two categories where the
  link is absent everywhere.

Relation to Banville et al. 2025: they give the **forward** probability that a link is
feasible/realised given ecological conditions; we place a **posterior on the
categorisation**. Their factors (co-occurrence, abundance, traits, environment,
sampling) govern realisation/detectability, which underlie our error rates.

---

## 3. Quantitative extension — accumulation of confidence

The accumulation figure/panel steps **outside** the binary Box 2: instead of a single
bit `O_rep` with one error rate `ε_rep`, it uses the **proportion of replicates** in
which the link is seen. This is what creates an accumulation to plot.

Per-replicate generative model for the contextual axis:
- A **truly realised-elsewhere** link (`z_rep = 1`) is **realised** in a replicate with
  probability `ρ` (realisation rate) and then **detected** with the same sampling
  success as locally, `1 − ε_local` — *a replicate is just another local sample*, so the
  per-replicate detection error **is** `ε_local` (there is **no separate detection
  parameter**). Net per-replicate detection probability:
  $$p_1 = \rho\,(1-\varepsilon_{\text{local}}).$$
- A **truly absent** link (`z_rep = 0`) yields a **false detection** with probability
  $$p_0 = f \quad(\text{false-detection / misID rate}).$$

Given `k` detections out of `R` replicates, the per-category score is the two binary
within-network factors times a binomial replicate term (the binomial coefficient is
common to all categories and cancels):

$$\text{score}(c)=\pi_c\;\underbrace{[\,1-\varepsilon_Y\ \text{or}\ \varepsilon_Y\,]}_{\text{match/mismatch of }\hat Y}\;
\underbrace{[\,1-\varepsilon_{\text{local}}\ \text{or}\ \varepsilon_{\text{local}}\,]}_{\text{match/mismatch of }O_{\text{local}}}\;
p_c^{\,k}(1-p_c)^{\,R-k},$$

$$P(C=c\mid E)=\frac{\text{score}(c)}{\sum_{c'}\text{score}(c')},\qquad
p_c = p_1\ \text{if}\ z_{c,\text{rep}}=1,\ \text{else}\ p_0.$$

**Two scenarios of consistent evidence** (both plotted vs `R`, the number of replicates):
- **Detected in every replicate:** `k = R`.
- **Empty in every replicate:** `k = 0`.

### Key results / behaviours (verified)
- **Ceiling.** As `R → ∞` the contextual axis is fully resolved, but replicate data
  never touches `ε_Y` or `ε_local`. The chosen category saturates at
  `(1−ε_Y)(1−ε_local)` (dashed line in the plots). With `ε_Y=0.20, ε_local=0.30` this is
  **0.56**. (Same value the binary headline confidence reaches as `ε_rep → 0`.)
- **Detected scenario, a pair on the O_rep axis.** For two categories sharing `(Ŷ,
  O_local)` (e.g. *possibly missing* vs *phantom*) everything cancels except the
  replicate term, so their posterior ratio is `(p_1/p_0)^R`. Hence:
  - `f < ρ(1−ε_local)` → possibly missing rises, phantom decays (expected).
  - `f = ρ(1−ε_local)` → the two curves coincide and stay flat.
  - `f > ρ(1−ε_local)` → **they flip** (phantom rises). Crossover at
    `f* = ρ(1−ε_local)`. E.g. ρ=0.3, ε_local=0.3 ⇒ f* = 0.21. This is correct Bayesian
    behaviour: if false detections outnumber true ones per replicate, a detection is
    evidence for the artifact (phantom).
- **ρ only matters when `f > 0`.** With `f = 0` (no misIDs) the surviving realised
  categories all share the same `p_1^R`, which cancels under all-detections, so ρ has no
  effect there (a single detection is conclusive). A nonzero `f` restores the
  ρ-dependent transient. → We keep `f` in the model for this reason.
- **Plateaus = the bar chart.** Under detections, the four "observed-elsewhere"
  categories rise and the four "not-observed-elsewhere" ones decay; each plateau equals
  the binary posterior in the `ε_rep → 0` limit. So the accumulation is literally the
  bar chart emerging as replicates accumulate. (Empty scenario is the mirror image.)
- **All four confusion-cell pairs behave identically** under a uniform prior (the
  within-network factors cancel in each pair's ratio and the normalisation is the same
  across cells). Breaking that symmetry requires a non-uniform prior or per-cell params.

---

## 4. Files in this folder

- **`box2_accumulation.R`** — ggplot2 figure of the accumulation curves. Parameters at
  the top (`eY, eL, rho, f, Rmax, chosen`). Current figure: two panels (detected / empty
  scenarios), all 8 categories, one distinct colour each (ColorBrewer Dark2), dashed
  ceiling line. Posterior is always normalised over all 8 categories. Greek labels use
  plotmath (render on any locale). Writes `box2_accumulation.{pdf,png}`.
  - Figure parameters: **ρ = 0.15, f = 0.05, ε_Y = 0.20, ε_local = 0.30**,
    `chosen = "Possibly missing"` (sets the observed Ŷ, O_local for the all-8 view).
- **`box2_accumulation.pdf` / `.png`** — rendered figure.
- **`link_taxonomy_explorer.html`** — standalone interactive explorer (no build step,
  no dependencies; all CSS/JS inline, canvas-based). Lets the user pick a category on a
  clickable taxonomy tree, edit the three error rates (+ a draggable phase-space
  heatmap), an editable per-category prior, a live likelihood table, and the accumulation
  panel (ρ, f sliders; detected/empty toggle). This is the source for
  contextualevidence.ecomplab.com. The 8-category accumulation colours match the R
  figure (Dark2).

### Reproduce the figure
```r
# from this folder
source("box2_accumulation.R")   # needs ggplot2; writes box2_accumulation.{pdf,png}
```

---

## 5. Open threads / things to reconcile
- Reconcile this standalone `box2_accumulation.R` with the repo's existing
  `code/posterior_accumulation.R` (real-data version) — shared helper for the posterior
  would avoid drift.
- The interactive currently fixes the *observed* pattern to the chosen category; the
  per-category prior editor is the main lever for asymmetry between categories.
- Figure caption (single-panel, detected-only, possibly missing vs phantom) is drafted
  in the manuscript; it explains the ceiling and that the posterior odds grow as
  `(p_1/p_0)^R = (ρ(1−ε_local)/f)^R`.
