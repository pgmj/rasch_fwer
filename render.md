# Preprint PDF (Typst) → index.pdf
quarto render index.qmd --to preprint-typst

# Blinded Word manuscript → index.docx
quarto render index.qmd --to apaquarto-docx --profile blind

# Split the docx into main text and online supplement
python3 split_appendix.py index.docx

The split writes index-main.docx and index-appendix.docx.

For an unblinded Word version, leave out --profile blind. Both renders write to index.*, so render the format you need last, or rename the output before the next render.
