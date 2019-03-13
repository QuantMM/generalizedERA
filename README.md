Generalized Extended Redundancy Analysis (gERA) in R
====================================================

Here we provide some basic R functions to implement generalized extended redundancy analysis (gERA)!

About gERA
----------

Extended redundancy analysis (ERA) is a statistical method for relating multiple sets of predictors to response variables. This method aims to extract a weighted composite or component from each set of predictors in such a way that it explains the maximum variation of response variables. Thus, ERA performs linear regression and data reduction simultaneously, providing a simpler description of directional relationships among many sets of variables.

Generalized ERA (gERA) is a flexible generalization of ordinary ERA! gERA handles various types of response variables (e.g., continuous, binary, or count) that are assumed to follow an exponential family distribution.

The official reference to ERA is the following paper:

-   Takane and Hwang (2005),

You also can find plenty of extensions:

-   Hwang (2000), Takane and Hwang (2005), Takane and Hwang (2007), Hwang et al. (2012), Tan, Choi, and Hwang (2015), Hwang et al. (2015), DeSarbo et al. (2015), Lee et al. (2016), and Lee et al. (2018)

Authors
-------

-   **Heungsun Hwang** - *Initial work in Matlab* - <heungsun.hwang@mcgill.ca>
-   **Sunmee Kim** - *Code maintainer and translator to R* - <sunmee.kim@mail.mcgill.ca>

Come visit our lab, here: [The Quantitative Methods Lab](https://sites.google.com/view/hwanglab/home?authuser=0)

References
==========

Wayne S. DeSarbo, Heungsun Hwang, Ashley Stadler Blank, and Eelco Kappe. 2015. “Constrained Stochastic Extended Redundancy Analysis.” *Psychometrika* 80 (2). Springer US: 516–34. doi:[10.1007/s11336-013-9385-6](https://doi.org/10.1007/s11336-013-9385-6).

Heungsun Hwang. 2000. “Structural equation modeling by extended redundancy analysis.” PhD thesis, McGill Univ.

Heungsun Hwang, Hye Won Suk, Jang-Han Lee, D. S. Moskowitz, and Jooseop Lim. 2012. “Functional Extended Redundancy Analysis.” *Psychometrika* 77 (3). Springer-Verlag: 524–42. doi:[10.1007/s11336-012-9268-2](https://doi.org/10.1007/s11336-012-9268-2).

Heungsun Hwang, Hye Won Suk, Yoshio Takane, Jang-han Lee, and Jooseop Lim. 2015. “Generalized functional extended redundancy analysis.” *Psychometrika* 80 (1): 101–25. doi:[10.1007/S11336-013-9373-X](https://doi.org/10.1007/S11336-013-9373-X).

Sungyoung Lee, Sungkyoung Choi, Young Jin Kim, Bong-Jo Kim, T2d-Genes Consortium, Heungsun Hwang, and Taesung Park. 2016. “Pathway-based approach using hierarchical components of collapsed rare variants.” *Bioinformatics* 32 (17). Oxford University Press: i586–i594. doi:[10.1093/bioinformatics/btw425](https://doi.org/10.1093/bioinformatics/btw425).

Sungyoung Lee, Yongkang Kim, Sungkyoung Choi, Heungsun Hwang, and Taesung Park. 2018. “Pathway-based approach using hierarchical components of rare variants to analyze multiple phenotypes.” *BMC Bioinformatics* 19 (S4). BioMed Central: 79. doi:[10.1186/s12859-018-2066-9](https://doi.org/10.1186/s12859-018-2066-9).

Yoshio Takane and Heungsun Hwang. 2005. “An extended redundancy analysis and its applications to two practical examples.” *Computational Statistics & Data Analysis* 49: 785–808. doi:[10.1016/j.csda.2004.06.004](https://doi.org/10.1016/j.csda.2004.06.004).

Yoshio Takane and Heungsun Hwang. 2007. “Regularized linear and kernel redundancy analysis.” *Computational Statistics & Data Analysis* 52 (1). Elsevier: 394–405. doi:[10.1016/j.csda.2007.02.014](https://doi.org/10.1016/j.csda.2007.02.014).

Tianyu Tan, Ji Yeh Choi, and Heungsun Hwang. 2015. “FUZZY CLUSTERWISE FUNCTIONAL EXTENDED REDUNDANCY ANALYSIS.” *Behaviormetrika* 42 (1). The Behaviormetric Society: 37–62. doi:[10.2333/bhmk.42.37](https://doi.org/10.2333/bhmk.42.37).
