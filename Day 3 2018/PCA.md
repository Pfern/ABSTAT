---
title: "Advanced Biostatistics 2018 - Practical Exercises"
#author: "carina"
#date: '2018-04-05'
output:
  html_document:
    keep_md: true
---

x <- rmarkdown::render("file.RMD", run_pandoc = FALSE, clean = FALSE)
knit_meta <- attr(x, "knit_meta") 
rmarkdown::render(input = 'file.knit.md', knit_meta = knit_meta )

<style type="text/css"> body, td { font-size: 18px; } code.r{ font-size: 18px; } pre { font-size: 16px } </style> 





Consider the iris dataset (included with R) which gives the petal
width, petal length, sepal width, sepal length and species for 150
irises. To view more information about the dataset, enter
help(iris).
Perform a PCA analysis to the iris data.




