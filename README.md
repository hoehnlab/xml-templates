# xml-templates

[Dowser](https://github.com/immcantation/dowser)-compatible BEAST2 XML templates for running [TyCHE](https://github.com/hoehnlab/tyche)
and non-type-linked analyses. These templates are
designed for use with Dowser 2.4 or later and BEAST2 2.7.7.

## Requirements

- BEAST2 2.7.7
- Dowser 2.4 or later
- TyCHE (for `TypeLinked/` templates)

## Using these templates

Download the templates you want (or clone the whole repo), and provide the
path to `getTimeTrees()`/`getTimeTreesIterate()`. See
[Build Time Trees Using TyCHE and BEAST](https://dowser.readthedocs.io) for
the full walkthrough, and
[Choosing an xml-templates variant](https://tyche.readthedocs.io) for how
the `StrictClock`/`UCLD`/`TypeLinked` families, `InstantSwitch`/
`ExpectedOccupancy`, `Est`/`FixedTraitClockRates`, `_EmpFreq`, and
`simplified/` naming pieces fit together.

Before you run one, it's worth knowing what a template actually needs from
you and what it supports -- that's what the rest of this page is about.

## Finding out what a template needs: `describeTemplate()`

Every template has some `${VAR}` placeholders Dowser fills in for you
automatically, and some free-form ones you have to supply yourself as named
arguments to `getTimeTreesIterate()`. Rather than opening the XML and
grepping for `${...}`, run this from R (dowser package, not this repo):

```r
describeTemplate("TypeLinked/TypeLinkedExpectedOccupancy_EstTraitClockRates_EmpFreq.xml")
```

```
TypeLinkedExpectedOccupancy_EstTraitClockRates_EmpFreq.xml

Structural elements (built-in Dowser behavior):
  [present] GermlineRootTree - Dowser will NOT add its own germline
            TipDatesRandomWalker operator if you pass
            include_germline_as_tip=TRUE; this template's own GRT operators
            handle it instead.
  [present] ROOTFREQS placeholder - include_germline_as_root=TRUE will take
            effect (root freqs enforced).
  [absent]  nodeTypes parameter (id='nodeTypes') - tree_states=TRUE
            with a starting tree will ERROR on this template.
  ...

Variables (16 found, 16 documented):
  KAPPA_PRIOR_M     [required]  Lognormal prior mean (M) for the HKY
                                 substitution model's transition/transversion
                                 ratio, kappa.
  TRAIT_RATE_MEAN_1 [required]  Prior mean clock rate for whichever cell/
                                 lineage type is coded 0 in your trait...
  ...

All 16 variables are documented in this template's dowser-template-info block.
```

This works on *any* template, documented or not -- the structural checklist
(the `[present]`/`[absent]` section) is built entirely from Dowser's own
substitution code, so it's accurate even for a template nobody's annotated
yet. It's the free-form variable descriptions (the second section) that
depend on the template author having documented them -- see below for how.

**Read this before assuming an argument you passed did anything.** Several
of the structural checks below are silent no-ops, not errors, when a
template doesn't support what you're asking for -- `describeTemplate()` is
how you catch that before a long MCMC run, not after.

## Making your own templates

Start from an existing template, or from an example XML exported from
BEAUti, and replace your alignment/trait data with the appropriate
variables (`${DATA}`, `${DATE}`, etc. -- see the table below for the full
set Dowser recognizes). Add a `dowser-template-info` block documenting any
free-form variable you introduce, then check your work:

```sh
Rscript .base/lint_templates.R your_new_template.xml
```

```r
describeTemplate("your_new_template.xml")
```

### Placeholders Dowser requires or handles specially

The behavior of these placeholders is fixed by Dowser's
code (`write_clone_to_xml()`), not by the template author.

| Placeholder / element | Required? | What happens |
|---|---|---|
| `${DATA}` | Always required | The sequence alignment. A template without this can't be used by Dowser at all. |
| `${CLONE}` | Always required | The clone/lineage ID, used to namespace IDs throughout the file. Also always required. |
| `${DATE}` | Optional | If present, `time=` sample dates get attached to tips. If absent, `time=` is silently ignored. |
| `${TRAIT}` / `${TRAIT_NAME}` | Optional | If present, `trait=`/`trait_data_type=` get attached and reconstructed. If absent, `trait=` is silently ignored. Templates with only `${TRAIT}` (no `${TRAIT_NAME}`) fall back to replacing `tag="location"` -- an older-style convention. |
| `${NODES}` | Optional | Replaces with the total number of nodes in the tree, used e.g. for dimensions of per-node parameters. |
| `${MRCA}` | Optional | If present, `start_date`/`max_start_date`/germline-tip MRCA priors take effect. If absent, those arguments are silently ignored. |
| `${OPERATORS}` | Optional | If present, lets Dowser add a germline tip-dating operator when `include_germline_as_tip=TRUE` (skipped, with a warning, if the template uses `GermlineRootTree` -- see below). |
| `${ROOTFREQS}` | Optional | If present, `include_germline_as_root=TRUE` takes effect (root nucleotide frequencies enforced via the external `rootfreqs` package). **If absent, this argument is silently ignored -- no error, no warning.** |
| `${EMP_EQFREQS}` | Optional | If present, empirical per-clone nucleotide frequencies are computed from your data at runtime. Expected to be absent from non-`_EmpFreq` templates. |
| `${ROOT_TYPE_PROBABILITIES}` | Optional | If present, `root_trait=` pins the root's discrete type. **This is the one case that warns rather than silently doing nothing** if absent. |
| `${NODES_TYPE_INIT}` | Optional | Sets the initial internal-node type index at MCMC start; defaults to `0` if omitted even when present. Only meaningful for type-linked templates. |
| `<init>...</init>` tag | Special behavior | A custom starting tree will replace this tag. |
| `id='nodeTypes'` parameter | Special behavior |  For Type Linked models, this is the parameter that stores the type of internal/ambiguous nodes. Automatically updated if you pass `tree_states=TRUE` with a starting tree. |
| `GermlineRootTree` | Special behavior | Germline handling automatically compatible with GermlineRootTree. |

There are also free-form variables, e.g. `TRAIT_RATE_MEAN_1`, `KAPPA_PRIOR_M`, `RATE_INDICATORS`. These are passed as named arguments (`...`) to
`getTimeTreesIterate()`. Template authors should document any free-form variables with the convention below.

### Documenting free-form variables: `dowser-template-info`

Every template should carry a single comment block documenting its
free-form variables, placed anywhere in the file (beginning or end -- the
parser looks for the marker text, not a fixed position):

```xml
<!-- dowser-template-info
@description | TyCHE type-linked clock using the expected-occupancy model, with each type's clock rate estimated during MCMC.
@var TRAIT_RATE_MEAN_1 | required | Prior mean clock rate for whichever cell type is coded 0 in your traitmap.
@var RATE_INDICATORS | required | "1 0" for irreversible type switching, "1 1" for reversible.
@var NODES_TYPE_INIT | optional, default 0 | Initial type assigned to internal nodes before MCMC starts.
@struct GermlineRootTree | Any template-specific note beyond the generic explanation above.
-->
```

Grammar: `@description | free-text summary` is an optional, one-per-file
overall summary of what the template is/does -- `describeTemplate()` prints
it first, before the structural checklist and variable list, and
`generate_catalog.R` (below) uses it to build a one-line-per-template table
for this README. `@var NAME | required|optional[, default X] | free-text
description`. `@struct ELEMENT | note` adds an extra template-specific note
onto one of the built-in structural rows from `describeTemplate()`'s
checklist -- most templates don't need this at all, since the built-in
checklist already covers the standard behavior. A line that isn't a new
`@description`/`@var`/`@struct` entry is treated as a continuation of the
previous entry's text, so long descriptions can wrap. Like `@var` and
`@struct`, `@description` always needs the pipe -- `@description text` with
no `|` isn't recognized.

Only one `@description` is read per template. If a `.base/` file's
description changes, re-run `make_all_templates_from_base.R` so the
generated `StrictClock/`, `UCLD/`, `TypeLinked/`, and `simplified/` copies
pick up the change too (descriptions aren't placeholder-specific, so they
always carry through unfiltered).

**Rules:**

1. **Never write a variable's name as `${NAME}` in the description text.**
   Dowser's placeholder substitution scans the
   whole file, including comments, and will fill in a literal
   `${NAME}` it finds there. Refer to variables by their bare name in prose
   ("the `EMP_EQFREQS` placeholder", not "the `${EMP_EQFREQS}` placeholder").
2. **Never write `--` anywhere inside the comment.** This is an XML rule, as `--` is illegal inside any XML comment's content. 

Run `lint_templates()` (below) to catch both of these automatically, along
with a documentation block that's out of sync with the file.

## Checking your templates: `lint_templates()`

`.base/lint_templates.R`, in this repo, checks a template's `dowser-template-info`
block against its actual placeholders. Run it from the command line, from the
repo root:

```sh
Rscript .base/lint_templates.R path/to/a/template.xml
# or lint every .xml file under a directory, recursively:
Rscript .base/lint_templates.R .
```

or from an R session, for the functions directly:

```r
source(".base/lint_templates.R")
lint_templates(".")
```

It fails on:

- a free-form `${VAR}` in the file with no matching `@var` entry
- an `@var` entry documenting a placeholder that isn't in the file 
- a `dowser-template-info` block that's opened but never closed with `-->`
- a literal `${NAME}` written inside the comment
- a `--` sequence anywhere in the block 
- `${DATA}` or `${CLONE}` missing entirely

It also prints, as **info, not as a failure**, if some of the optional
structural elements from the table above are absent. This is info is designed
for developing base templates that are maximally compatible and is not complete.

## Keeping the catalog in sync: `generate_catalog.R`

`.base/generate_catalog.R` collects every template's `@description` (see
above) into the table under [Template catalog](#template-catalog) below,
grouped by folder. Run it after adding or changing a `@description`:

```sh
Rscript .base/generate_catalog.R .
```

It rewrites everything between the `BEGIN TEMPLATE CATALOG`/
`END TEMPLATE CATALOG` comments in this README -- edits made by hand inside
that block will be overwritten the next time this runs. A template with no
`@description` yet still gets a row, marked `(no @description yet)`, rather
than being left out silently; that's a quick way to see which templates
still need one. Like `lint_templates.R`, this only reads whatever
`@description` is actually in each file -- if you changed a `.base/` file,
run `make_all_templates_from_base.R` first so the generated copies pick up
the change before cataloging them.

## Contributing to this repo/Developing new base templates

If you're re-running `make_all_templates_from_base.R` after editing a
`.base/` file's `dowser-template-info` block, note that the generator
automatically drops any `@var` entry whose placeholder gets filled in by
that generation step (e.g. `simplified/`'s baked-in defaults). You don't
need to hand-maintain separate doc blocks for every generated variant.

If you're adding a new base template intended to generate `simplified/` and
`_EmpFreq`-stripped variants the way the existing `TypeLinked/`,
`StrictClock/`, and `UCLD/` families do, add it under `.base/` and see
`make_all_templates_from_base.R` for the generation contract.

## Template catalog

<!-- BEGIN TEMPLATE CATALOG (auto-generated by generate_catalog.R; do not edit by hand) -->

<!-- END TEMPLATE CATALOG -->
