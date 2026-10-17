# Legacy mathematical notes

`document.tex` and `images/` are the original RandLib mathematical reference. They describe the 1.x
API and may contain obsolete implementation or accuracy statements, so they are retained as source
material rather than published as the current manual. Corrected definitions and descriptions of the
v2 algorithms are being migrated into the [current mathematical reference](../reference/index.md).

Build from this directory with a suitably complete LaTeX installation:

```sh
pdflatex document.tex
pdflatex document.tex
```

The generated PDF is intentionally not versioned. Current API documentation is generated from the
2.x public headers with Doxygen.
