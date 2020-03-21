#!/bin/bash

# Build the paper with analyses, images, tables, references etc.

# Run analyses
# cd analysis
# Rscript analyses_plots.R
# Rscript analyses_stats_final_final.R
# cd ..

# Fix table snippets

## Remove $ around numbers
sed -i 's/\$//g' manuscript/include/*.tex

## Align numerical columns on decimals
sed -i 's/cccc/lS[table-format=3.2]S[table-format=3.2]S[table-format=3.2]/g' manuscript/include/gls_aic_comp.tex

sed -i 's/cccccc/llS[table-format=3.2]S[table-format=3.2]S[table-format=3.2]S[table-format=3.2]/g' manuscript/include/*.tex

## Add captions
sed -i 's/caption{}/caption{Model comparison of Generalised Least Squares models predicting damaged sapling bark area using different spatial autocorrelation structures. Models are ordered by increasing AIC value.}/g' include/gls_aic_comp.tex

sed -i 's/caption{}/caption{Model comparison of logistic generalised linear mixed effects models predicting the likelihood of a sapling being attacked by \\textit{H. abietis}. Models are sorted according to increasing AIC.}/g' include/binom_lmer_aic_comp.tex

sed -i 's/caption{}/caption{Model comparison of general linear mixed effects models predicting the damaged bark area of a sapling, for those saplings which have been initially damaged. Models are sorted according to increasing AIC.}/g' include/lmer_aic_comp.tex

## Put column headers in curly braces
sed -i '10s/.*/{Cor. Struct.} \& {AIC} \& {logLik} \& {r2c} \\\\/' manuscript/include/gls_aic_comp.tex

sed -i '10s/.*/{Fixed eff.} \& {Random eff.} \& {AIC} \& {logLik} \& {r2c} \& {r2m} \\\\/' manuscript/include/binom_lmer_aic_comp.tex

sed -i '10s/.*/{Fixed eff.} \& {Random eff.} \& {AIC} \& {logLik} \& {r2c} \& {r2m} \\\\/' manuscript/include/lmer_aic_comp.tex

## Replace shorthand column headers with longhand
sed -i 's/r2m/R\\textsuperscript{2}\\textsubscript{m}/g' manuscript/include/*.tex

sed -i 's/r2c/R\\textsuperscript{2}\\textsubscript{c}/g' manuscript/include/*.tex

# Crop images
convert analysis/damage_high.jpg -gravity East -chop 120x0 manuscript/img/damage_high_crop.jpg
convert analysis/damage_low.jpg -gravity South -chop 0x60 manuscript/img/damage_low_crop.jpg

# Run LaTeX and BibTeX
cd manuscript
BASE="${1%.*}"
pdflatex weevils.tex
if [ $? -ne 0 ]; then
    echo "Compilation error. Check log."
    exit 1
fi
pdflatex weevils.tex
bibtex weevils
bibtex weevils
pdflatex weevils.tex
pdflatex weevils.tex
pdflatex weevils.tex
latexmk -c weevils.tex
exit 0


# Create a word document version for distribution
pandoc -f latex -t docx weevils.tex -o weevils.docx
