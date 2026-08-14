# Prior tab — copy for review

Edit the text blocks below in place. Once you're happy, I implement them verbatim.
Nothing here is in the page yet.

**The rule we agreed:** only two kinds of text on screen.

1. **How to use and read the plots** → the card hint, directly above the controls.
2. **The mapping of ecology to the framework** → the lede under the pill, before the card.

Everything else — closed-form boundaries, why the regions are rectangles, the
`Σ P(C=c) = 1` identity, the π_l < π_r convention, the 80% invariant — moves to the
(i) button, which already contains all of it word for word.

| block | now | after |
|---|---:|---:|
| lede under the pill | 117 | 68 |
| card hint | 95 | 60 |
| slider descriptions | 110 | 55 |
| notes under the map | 184 | **0** |
| caption under the cut | 55 | 35 |
| **total on screen** | **561** | **~218** |

Axis titles, legends, the heading and the in-figure labels are *not* in that count —
they are block 4, listed separately for review.

---

## 1. Lede under the pill — *the ecology*

Sits above the card, before any control. Job: what π_l and π_r **mean** ecologically.
(The other two pills use their lede the same way — to say what the version is.)

> Ecological knowledge provides beliefs about a species pair (e.g., a plant and a bee) before
> any evidence: how likely they interact **anywhere in the region** (π_r), and **at this
> site** (π_l). Such priors can be informed for instance, by species abundance, degree or traits. In this example the map runs from two specialists (priors closer to 0, indicateing we do not believe they will interact) to two generalists (priors closer to 1).

---

## 2. Card hint — *how to use and read*

Sits under the card heading, directly above the sliders and plots. Job: what to do, and
how to read what happens.

> **Drag the marker.** A block's colour is the category that belief makes most probable,
> so a boundary between blocks is where one reading gives way to another. The plot on the right follows
> all four categories along the dotted diagonal path through the marker: read it for how big a
> shift in belief it takes to change the answer. Evidence and error rates stay fixed. The vertical orange line on the right plot corresponds to the marker on the left plot

---

## 3. Slider labels and descriptions

Labels are the short bold text beside each slider; descriptions sit under it.

| label (current) | keep? |
|---|---|
| `π_r  regional prior` | |
| `π_l  local prior` | |
| `π_Y  model prior` | |
| `R  replicates` | |

**π_r — regional prior**

> Prior that the link is realisable in the replicate systems (e.g., Rises with regional degree).

**π_l — local prior**

> Prior that the link is realised where we sampled (e.g., Rises with local degree).

**π_Y — model prior**

> Prior that a truth-fitted model would predict the link (could be a scaled model output). Scales the four curves
> together; the map does not move.

**R — replicates**

> How many replicates stand behind the single bit O_r, which records only whether the
> link was seen in **at least one** of them — not how many.

---

## 4. Text inside and around the figures

Not counted in the table above, and not yet reviewed. All of it is current wording —
edit anything you want changed, leave the rest.

### Card heading

> **Which source explains the evidence?** — an informative prior on the two ecological truths

### Sweep toggle (sits above the right-hand plot)

> `Vary the pair`  |  `Vary the site`

These were "Vary the pair (SI a)" / "Vary the site (SI b)" until the stand-alone pass.
They now say what they do rather than which published panel they reproduce. An
alternative pairing, if these read as too terse:
`Vary how generalist the pair is` | `Vary how much transfers to this site`

### Map — axis titles

> x: **π_r — realisable in the replicates**
> y: **π_l — realised here**

### Map — labels drawn on the boundary lines

> **π_r\* = 0.347**  (on the vertical line)
> **π_l\* = 0.760**  (on the horizontal line)

Numbers are computed. The `π_d*` notation is only explained in the (i) button — if the
notes under the map go (block 6), these become the sole unexplained symbols on screen.
Options: leave as is, drop to unlabelled lines, or relabel as e.g. `switches here`.

### Map — legend

> *(one swatch per category)* — the leading category among the four the model bit supports

### Cut — axis titles

> x, when varying the pair: **π_r, with π_l/π_r held fixed**
> x, when varying the site: **π_l/π_r, with π_r held fixed**
> y: **Posterior probability P(C | E)**

### Cut — legend

Four category names with their line styles. No extra wording.

---

## 5. Caption under the cut — *generated, adapts to the leaf*

The **bold** parts are computed live. The rest is fixed wording.

> Evidence: predicted, not recorded here, **recorded in at least one of 5 replicates** —
> the signature of *possibly missing*. The leader changes at **0.35** and **0.84**. These
> four categories hold **80.0%** of the posterior between them; the other four are their mirrors.
> In this example, both priors rise with the pair's degree, so the axis runs from two interacting specialists on the left to two interacting generalists on the right. Dotted lines mark where the leading category changes. The same evidence and the same error rates therefore support **three** different readings, depending only on what is ecologically known about the pair.

Variants it has to produce:

- `O_r = 1` → "recorded in at least one of *R* replicates"
- `O_r = 0` → "recorded in none of *R* replicates"
- no switch on the cut → "One category leads across the whole cut."
- one switch → "The leader changes at **X**."

---

## 6. Deleted — both notes under the map

Cut entirely. Every sentence already appears in the (i) button.

~~Every point is a valid prior — the two axes are independent, so the whole square is
available. The shaded band below the diagonal is the ecologically usual case, π_l < π_r —
a convention rather than a constraint. The replicates are *other* places, so a link can be
realised here and nowhere else, which is exactly what **locally unique** is. Regional pools
are broader than local realisations, so π_r exceeds π_l on average, but for a locally
dominant specialist it need not. Drag the marker, or use the sliders. The dotted path is
the cut drawn on the right.~~

~~**The regions are rectangles**, and the two solid lines are why: the prior and the
likelihood both factorise across the axes, so each axis is settled by its own comparison
and neither knows about the other. A boundary sits where the prior odds exactly offset that
axis's likelihood ratio, π_d* = L_d(0)/[L_d(0)+L_d(1)]. The three regions of the cut on the
right are one straight path crossing these two perpendicular lines.~~

Two things from these are load-bearing and have been folded elsewhere rather than lost:

- "Drag the marker" → now opens the card hint (block 2).
- "The dotted path is the cut drawn on the right" → now "the dotted path through the
  marker" in the card hint.

---

## 7. Flag — not a page issue, a paper one

The replicate evidence in this panel is the **single bit** `O_r`, recorded in **at least
one** of R replicates. That is what `SI_prior_l.R` computes and what the page does.

The SI caption for the figure says "recorded in **every one** of R = 5 replicates", which
describes the count `n = R` — a different and much stronger piece of evidence:

```
binary bit, seen in >=1 of 5 :  0.4257 vs 0.2262  -> evidence ratio   1.88
count n = 5, seen in all     :  1.28e-5 vs 3.13e-7 -> evidence ratio  40.84
```

The main text has it right ("the binary replicate evidence after five replicates"). The
figure's published numbers (0.347, 0.844) come from the binary route, so the script is
correct and the SI caption is the thing to fix.

Note the other two tabs *do* use consistent evidence (`n = R` or `n = 0`, count retained).
The prior tab is deliberately the odd one out, matching the published figure — which is
why block 3's R slider and block 5's caption now say so explicitly.

---

## 8. Afterwards

If this lands well, the same treatment probably wants applying to the other two tabs — the
effort panel especially still carries a long note under its plot. Separate pass.
