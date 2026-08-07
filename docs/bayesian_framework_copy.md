# Page copy — Link taxonomy uncertainty explorer

Agreed text for `docs/bayesian_framework_interactive.html`. Notation follows the
Supplementary Information (`SI.tex`). This file is the source of truth for the
page's wording: settle the text here first, then implement it in the HTML.

Status: **Title + intro — agreed. Prior info modal — agreed.** Remaining blocks
(panel headers, slider descriptions, notes paragraphs, remaining info modals)
not yet drafted.

---

## Title

> **How confident are we in a *possibly missing* link?**

The italicised category name is dynamic — it updates as the reader clicks leaves
on the taxonomy tree.

---

## Intro

Three labelled paragraphs: what it is → what it assumes → how to use it. The
bold labels are lead-ins, not headings, and run inline with the text.

> **What this is:** an interactive guide to the Bayesian framework developed in
> the Supplementary Information of *Contextual evidence for classifying
> interactions and guiding their discovery*. The taxonomy assigns each candidate
> link to one of eight categories from three pieces of evidence — the model
> prediction Ŷ, the local observation O_l, and the replicate observation O_r.
> That assignment is deterministic, but every source is error-prone, so it
> carries uncertainty. The framework therefore treats the true category C as
> unknown and returns a posterior P(C | E) over the eight.
>
> **Scope:** the guide follows the simplest version of the framework — binary
> evidence, one symmetric error rate per source, a uniform prior, and
> independent errors — with the illustrative rates of the worked example
> (ε_Y = 0.2, ε_l = 0.3, ε_r = 0.1, ρ = 0.15, f = 0.05), chosen to show the
> behaviour of the framework rather than to describe any particular system. The
> Supplementary Information states these assumptions explicitly and shows what
> each would take to relax.
>
> **How to use:** pick a category on the taxonomy tree. Its error-free signature
> z_c becomes the observed evidence E, and the confidence in that label is
> (1−ε_Y)(1−ε_l)(1−ε_r). Drag the sliders to see how each error rate moves it,
> and where the residual mass goes.

### Notes on this block

- **How to use** is deliberately last so that it sits directly above the
  taxonomy card — the instruction to "pick a category on the taxonomy tree" is
  adjacent to the thing it points at.
- **Scope** was chosen over *Assumptions* because the paragraph does three
  things: states the assumptions, flags the rates as illustrative, and points to
  the SI. *Assumptions* would cover only the first.
- No Box references. The main text no longer has Boxes, so the page cites the
  Supplementary Information and Fig. 1A instead.
- `z_c` in the **How to use** paragraph is rendered as inline code in the page
  and updates with the selection (e.g. `z_c = (1, 0, 1)`).

---

## Category prior card — info modal

Placement: an **(i) button on the prior card's summary bar**, matching the
taxonomy card. This text was previously an inline `.explain` box at the foot of
the card; moving it to a modal keeps the card to its sliders and short hint.

Modal title: **What makes a prior non-uniform?**

> The uniform prior places half its mass on the link being present locally,
> since four of the eight categories have z_{c,l} = 1. More generally, the prior
> takes the same product form as the likelihood,
>
> > P(C=c) = Π_{d∈{Y,l,r}} π_d^{z_{c,d}} (1−π_d)^{1−z_{c,d}}
>
> with π_Y = π_l = π_r = 1/2 recovering the uniform prior. However, ecological
> networks are sparse and have skewed degree distributions, so most interactions
> are unlikely a priori while a few species have a high potential to interact.
>
> Two conditions must be met for an interaction: individuals must **meet**,
> which means they need to co-occur in space and time, and once they meet their
> **traits must match**. In setting priors, abundance acts on the first
> condition, since the probability that two species encounter each other rises
> with how common each of them is; generalism, which depends on both abundance
> and traits, acts in the same direction, because a species with many partners
> is more likely to have any particular one.
>
> Connectance sets the floor for π_l and π_r alike, a trait mismatch lowers
> both, and partners that do not overlap here but do elsewhere lower π_l alone.
> We leave π_Y uniform, since the model score is already evidence. The eight
> weights in this panel are free, so any such prior can be entered directly.

### Notes on this block

- The product form is the **general** prior; the uniform prior is the special
  case π_d = 1/2. An earlier draft had this backwards ("with a uniform prior,
  the prior takes the product form … with π = 1/2 recovering the uniform
  prior"), which is circular.
- Generalism earns its place only by saying what it does to the prior — a
  species with many partners is more likely to have any particular one — not by
  stating that generalism depends on abundance and traits.
- π_Y is explicitly disposed of (left uniform, because the model score is
  already evidence); otherwise it appears in the formula and never returns.
- The closing sentence connects π_d back to the eight free category sliders,
  which are what the panel actually exposes. Phrased "in this panel" rather than
  "above" because a modal has no spatial relation to the card.
- Deliberately omitted, available if wanted: trait-based exclusions are labile,
  since within-species trait variation lets pairs judged forbidden from database
  means interact — a caution against setting π too hard on trait mismatch.
