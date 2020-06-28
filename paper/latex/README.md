Based on [Style and Template for Preprints](https://www.overleaf.com/latex/templates/style-and-template-for-preprints-arxiv-bio-arxiv/fxsnsrzpnvwc).

Use [arxiv-collector](https://github.com/djsutherland/arxiv-collector) to generate files for arXiv. 

Just run: 

```
arxiv-collector ilmrisk.tex --latexmk ./path/to/latexmk
```

The argument ```--latexmk ./path/to/latexmk``` is needed only if the verson of latexmk in TexLive is 4.63b, arxiv-collector can download the latest version of latexmk and save it in an arbitrary location for use.