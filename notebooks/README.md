# Jupyter tools

Some tools to style and export `jupyter` notebooks to `HTML` reports and slides.
It uses `pandoc` to render citations from `bibtex` files.
- report: optional code visibility (see this [repository](https://github.com/csaid/polished_notebooks))
- slides: inspired by this [project](https://github.com/datitran/jupyter2slides)

Predefined `tex` macros included in the notebooks. There is probably a way to do better (eg load macros from a `tex` file see [this nb extension](https://github.com/ipython-contrib/jupyter_contrib_nbextensions/tree/master/src/jupyter_contrib_nbextensions/nbextensions/load_tex_macros))

Use `rm_in` and `rm_out` tags on cells where you want to remove input or output.


### Usage

For report:

```bash
nbcv --file nb.ipynb --bib ./tex/biblio.bib
```


For slides:

```bash
nbcv --slides --file pres.ipynb --bib ./tex/biblio.bib
```
