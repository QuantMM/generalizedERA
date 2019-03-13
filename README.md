---
title: "README"
author: "Sunmee Kim"
date: "March 12, 2019"
output: html_document
bibliography: Refer.bib
---

# Generalized Extended Redundancy Analysis (gERA) in R

Here we provide some basic R functions to implement generalized extended redundancy analysis (gERA)!

## About gERA

Extended redundancy analysis (ERA) is a statistical method for relating multiple sets of predictors to response variables. This method aims to extract a weighted composite or component from each set of predictors in such a way that it explains the maximum variation of response variables. Thus, ERA performs linear regression and data reduction simultaneously, providing a simpler description of directional relationships among many sets of variables.

Generalized ERA (gERA) is a flexible generalization of ordinary ERA! gERA handles various types of response variables (e.g., continuous, binary, or count) that are assumed to follow an exponential family distribution.

The official reference to ERA is the following paper:

- @Takane2005, 

You also can find plenty of extensions:

- @Hwang2000, @Takane2005, @takane2007regularized, @Hwang2012, @Tan2015, @Hwang2015, @DeSarbo2015, @Lee2016, and @Lee2018

## Authors

-   **Heungsun Hwang** - *Initial work in Matlab* - heungsun.hwang@mcgill.ca
-   **Sunmee Kim** - *Code maintainer and translator to R* - sunmee.kim@mail.mcgill.ca

Come visit our lab, here: [The Quantitative Methods Lab](https://sites.google.com/view/hwanglab/home?authuser=0)

# References
