# Notebook Rendering

The conos tutorial notebooks (`*.ipynb`) are rendered from `*.Rmd` with **Quarto**, not with a bare
`rmarkdown::render()` + ad-hoc converter. Quarto emits real notebook code cells for the R chunks.

GitHub's notebook viewer is sensitive to the notebook `language_info` metadata: with R syntax metadata it
renders the code input inside a CodeMirror container that **collapses indentation/whitespace**. The fix is
to keep an R kernelspec but set `language_info` to plain text, so GitHub wraps each code cell in a `<pre>`
block and preserves layout.

## Always use `render_notebook.sh`

[`render_notebook.sh`](render_notebook.sh) runs `quarto render --execute` **and** applies the
`language_info`/kernelspec fixup in one step, so the fixup is never forgotten:

```sh
doc/render_notebook.sh conos-walkthrough        # novice / common-case tutorial
doc/render_notebook.sh conos-advanced           # advanced: alignment spaces, facets, label transfer
```

Each call writes `<name>.ipynb` (GitHub-viewable) and `<name>.html`.

## Requirements

* **Quarto** (tested with `1.9.38`; the script expects it under `$HOME/.local/quarto-1.9.38/bin`).
* The **R `IRkernel`** (`ir`) and a working `jupyter`, used by Quarto to execute R cells.
* **conos** and **pagoda2** installed (the notebooks `library(conos)` / `library(pagoda2)`), plus
  [`conosPanel`](https://github.com/kharchenkolab/conosPanel) for the example data:
  `install.packages("conosPanel", repos="https://kharchenkolab.github.io/drat/", type="source")`.

> Render against the **same conos/pagoda2 versions you intend to ship**: the notebooks execute live code,
> so the rendered output reflects the installed packages' API (e.g. the `run*` verbs).
