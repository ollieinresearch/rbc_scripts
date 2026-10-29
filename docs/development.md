# Working on the code and these pages

## How this site is made

These pages are built by [Sphinx](https://www.sphinx-doc.org), the documentation generator
behind NumPy, SciPy and Dedalus, with the [PyData theme](https://pydata-sphinx-theme.readthedocs.io)
(the look of NumPy, SciPy and pandas). Sphinx reads three kinds of source and writes a
static website, a folder of HTML files that any web server, or a browser opening the files
directly, can show:

1. **Pages** in `docs/`, written in Markdown (with the [MyST](https://myst-parser.readthedocs.io)
   extensions for cross-references, maths and directives).
2. **Docstrings** of the `rbc` package, imported by Sphinx to build the
   {doc}`reference/api`.
3. **Command-line parsers** (`get_parser()` in each module), read to build the
   {doc}`reference/cli`, and **example scripts** in `docs/examples/`, run at build time to
   make the {doc}`gallery/index`.

So the reference pages can never fall out of date with the code: they are the code.

| Extension | What it adds |
|---|---|
| `myst_parser` | Markdown pages |
| `autodoc`, `autosummary`, `numpydoc` | API pages from NumPy-style docstrings |
| `sphinxarg.ext` | command-line pages from the argparse parsers |
| `intersphinx` | links to the NumPy, SciPy, Matplotlib, pandas, h5py, mpi4py and Dedalus docs |
| `viewcode` | `[source]` links to highlighted source code |
| `mathjax` | LaTeX maths |
| `sphinx_design` | grids, cards, tabs, dropdowns, badges |
| `sphinx_copybutton` | a copy button on every code block |
| `sphinx_tippy` | hover previews of links, glossary terms and equations |
| `sphinx_gallery` | the figure gallery |

`docs/conf.py` holds all the settings.

## Building the pages

```bash
pip install -e ".[docs]"                          # once
sphinx-build -b html docs docs/_build/html        # build
python -m http.server -d docs/_build/html         # view at http://localhost:8000
```

or, while writing, `pip install sphinx-autobuild` and
`sphinx-autobuild docs docs/_build/html`, which rebuilds and reloads the browser on every
save. Add `-W --keep-going` to treat warnings (a broken cross-reference, say) as errors.
A full build runs the gallery examples and takes about a minute; delete `docs/_build` to
start clean.

## Writing pages

A page is a Markdown file in `docs/`, listed in the `toctree` of its section's
`index.md`. The essentials:

````markdown
# Page title

Text with a {term}`glossary term`, a link to a {doc}`page <user-guide/storage>`, to a
labelled section {ref}`graceful-stop`, to a function {func}`rbc.cases.slices` or
module {mod}`rbc.collect`, and maths $Nu - 1 \sim Ra^{\beta}$.

(my-label)=
## A section others can link to with {ref}`my-label`

$$
Nu = 1 + \sqrt{Ra\,Pr}\,\langle wT\rangle
$$

```{dropdown} Derivation
:class-container: derivation

Hidden until clicked.
```

```{note}
Also tip, warning, important, seealso.
```

```{code-block} bash
:emphasize-lines: 2

echo one
echo two      # highlighted
```
````

Side-by-side comparisons use `::::{grid}` with `:::{grid-item-card}`; see the source of
{doc}`user-guide/clusters`. New terms go in `docs/glossary.md`.

## Docstrings

Every function has a NumPy-style docstring; these become the API pages:

```python
def slice_c(log_ra, log_pr, a=-0.5):
    """Label of the slice through a point: c = log10 Pr_0.

    Longer explanation, if needed.

    Parameters
    ----------
    log_ra, log_pr : float or array_like
        log10 Ra and log10 Pr of the point(s).
    a : float, optional
        Slice slope.

    Returns
    -------
    numpy.ndarray or float
        c for each point.
    """
```

## Publishing on GitHub Pages

The site is published by GitHub Actions on every push to `main`
(`.github/workflows/docs.yml`): it installs the package with the docs extras, builds the
site, and deploys it. To switch it on once:

1. Push the repository to GitHub. On a free account, GitHub Pages needs the repository to
   be **public**; a site published from a private repository needs a paid plan, and is
   public on the web even then. `data/results.csv` is git-ignored, so no results are
   published.
2. In the repository on GitHub: **Settings → Pages → Build and deployment → Source:
   GitHub Actions**.
3. Set `GITHUB_URL` in `docs/conf.py` to the repository's address.
4. Push. The **Actions** tab shows the build; when it is green the site is at
   `https://<user>.github.io/<repository>/`.

A failed build (for example a broken example) does not replace the live site.

[Read the Docs](https://about.readthedocs.com) is the other common host (free for public
repositories, with advertisements). It would need a `.readthedocs.yaml` file instead of the
workflow; nothing else changes.

## Optional features not switched on

`Edit this page` button
: A link on every page to the page's Markdown source on GitHub, opening GitHub's editor.
  Anyone with write access can fix a typo from the browser, and the change goes through
  git like any other. Turn on with `"use_edit_page_button": True` in `html_theme_options`
  and `html_context = {"github_user": ..., "github_repo": ..., "github_version": "main",
  "doc_path": "docs"}`.

Changelog
: A page listing, per version, what changed for users ("`converged_stat` added", "AVG_TIME
  accepts auto"). Unlike commit messages it is written for readers and grouped by release.
  Useful once others depend on specific versions, e.g. the version used for the paper.

Link check
: `sphinx-build -b linkcheck docs docs/_build/linkcheck` visits every external link and
  reports the broken ones. Run it by hand now and then, or as a scheduled workflow; it
  does not change the site.

(adding-a-figure)=
## Adding a figure

1. Create `rbc/plots/<name>.py` with a module docstring (first line = description), a
   `get_parser()` returning an `argparse.ArgumentParser` with `--results`, `--out`,
   `--png` and `--no-tex`, and `main(argv=None)` that calls `style.use("paper")`, builds
   the figure at `style.figsize(...)` and saves with `style.save`.
2. Load data with {func}`rbc.cases.load_results`, and colour slices with
   {func}`rbc.style.slice_kw`, curves at one Ra or one Pr with
   {func}`rbc.style.category_kw`, and runs ranked within the figure with
   {func}`rbc.style.ranked_colours` (lines through {func}`rbc.style.line_kw`, so the pale
   colours get their outline), so it matches the other figures ({ref}`expense-colours`).
3. Add it to the `jobs` list of {mod}`rbc.plots.all`.
4. Add `docs/examples/plot_<name>.py` (copy an existing one) for the gallery, an entry in
   {doc}`user-guide/figures`, and the command to {doc}`reference/cli`.

## Conventions

- Comments are sparse: explain why, not what. Docstrings are complete.
- Canadian spelling in prose, docstrings and labels: colour, centre, behaviour; but
  normalize, analyze, visualize. Code identifiers follow their libraries (`color=`).
- Every command is `python -m rbc.<module>`, with a `get_parser()`.
- Figures are made at printed size; change `rbc.style.TEXTWIDTH_IN`, not individual
  figure sizes.
- Data paths come from {mod}`rbc.cases` (`DATA_DIR`, overridable with `$RBC_DATA`), never
  hard-coded.
