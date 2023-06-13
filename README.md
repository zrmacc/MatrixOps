---
title: "README"
author: "Zachary McCaw"
date: "2023-06-13"
output: 
  html_document: 
    keep_md: TRUE
--- 

# Matrix Operations

Zachary R. McCaw <br>
Updated: 2023-06-13




## Linear regression


```r
set.seed(1010)
n <- 1e4
x <- data.matrix(stats::rnorm(n))
y <- stats::rnorm(n)
microbenchmark::microbenchmark(
  stats::lm(y ~ x),
  MatrixOps::FitOLS(y, x)
)
```

```
## Unit: microseconds
##                     expr      min       lq      mean   median       uq      max
##         stats::lm(y ~ x) 2151.353 2333.980 2940.1124 2418.534 2659.827 14376.97
##  MatrixOps::FitOLS(y, x)  124.338  133.044  912.6233  266.896  291.297 69481.43
##  neval cld
##    100   b
##    100  a
```
