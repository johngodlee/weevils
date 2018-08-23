% Stats notes for weevils project
% John Godlee

# O'Hara & Kotze 2010

-   Generalised linear models would be better at dealing with count data

-   How do you deal with zeros

    -   Normally done by adding an arbitrary integer (usually 1) to the
        data set to fudge it.

    -   zero-inflated models are also a good way of handling lots of
        zeroes (Zuue, Ieno & Elphick 2010)

-   in comparisons, transformations performed poorly, except when
    dispersion was small and the mean counts were large

    -   Quassi-Poisson and negative binomial models consistently
        performed well, with little bias.

        -   These distributions are often specified as an alternative to
            the traditional poisson distribution when dealing with count
            data as overdispersion is common, they add an additional
            parameter that allow the variance to vary independently from
            the mean.

        -   Quassi-poisson assumes $v = \sigma\lambda$

        -   negative binomial assumes:
            $v = \lambda + \lambda^{2}/\theta$

        -   $\theta$ and $\sigma$ are both overdispersion parameters

        -   Read Ver Hoef & Boveng (2007)

-   Textbooks on statistical methodology in ecology recommend the use of
    the SQRT transformation


# Zeileis *et al.* 2008 - Regression models for count data in R

-   The classic poisson regression is not suitable for count data due to
    over-dispersion and lots of zeroes.

-   glm() - the Generalized linear model function

-   use the models hurdle() and zeroinfl() in the pscl package

-   Although quassi-poisson, negative binomial, and other Generalised
    Linear models are fine for overdispersion, they don’t deal well with
    excess zeroes

-   *zero-augmented* models capture the zero counts in a second model
    component

-   *Hurdle* models combine a left-truncated count component with a
    right-censored hurdle component

-   *zero inflated* models are mixture models that combine a count
    component and a point mass at zero.


# Chang & Pocock 2000 Analysing data with clumping at zero An example demonstration

-   Uses 2 approaches to analyse the relationship of multiple covariates
    to an outcome with lots of zeroes

    -   Categorise the continuous outcome and fit a proportional odds
        model

    -   Use a logistic regression to model the probability of a zero
        response and ordinary linear regression to model the non-zero
        continuous responses

-   The data used here appears to be log-normal.

# Discussion with Kyle - 2015_11_30

Packages

-   Lme4

-   Lmer

-   Glmer

-   Vegan – pcnm()

Incorporating the spatial influence of other cells

-   Progressively increase the radius of influence and evaluate the fit
    with AIC

-   Distance matrix with x,y coordinates of the trees, use this to make
    a pairwise distance matrix and add this as a covariate to determine
    their effect

-   Use a Moran’s correlogram/Semivariogram to determine when the
    spatial autocorrelation ceases

    -   Then use a principal coordinate analysis

-   Remember to incorporate edge effects

Extras

-   Remember to ID the blocks as factors rather than A,B,C,D

Zero-inflated models/Hurdle Models

-   Allows zeroes to be removed from the model-

-   Model in 2 parts, is there a difference between zeroes and ones, if
    a 1, what causes the difference.



