% Stats notes for weevils project
% John Godlee

# Summary

__Model types that deal with zeroes:__

* Zero inflated models (Zuur et al. 2010)
* Quasi-Poisson GLM
* Negative Binomial GLM
* Hurdle (two-stage) models
    * A model in two parts, a binomial model for differences in proportion of zeroes and greater than zero, then if greater than one, a linear model 

R functions to investigate:

* vegan::pcnm()
* pscl::zeroinfl()
* pscl::hurdle()

__Dealing with spatial influence of other cells:__

* Progressively increase the radius of influence then evaluate fit with AIC - but I'm not sure what the fit is meant to be of. 
* Make a pairwise distance matrix and add as a covariate to models to determine their effect
* Remember to incorporate edge effects


# O'Hara & Kotze 2010

- Generalised linear models would be better at dealing with count data

- How do you deal with zeroes?

    - Normally done by adding an arbitrary integer (usually 1) to the data set to fudge it.

    - zero-inflated models are also a good way of handling lots of zeroes (Zuue, Ieno & Elphick 2010)

- In comparisons, transformations performed poorly, except when dispersion was small and the mean counts were large

    - Quassi-Poisson and negative binomial models consistently performed well, with little bias.

        - These distributions are often specified as an alternative to the traditional poisson distribution when dealing with count data as overdispersion is common, they add an additional parameter that allow the variance to vary independently from the mean.

        - Quassi-poisson assumes $v = \sigma\lambda$

        - negative binomial assumes: $v = \lambda + \lambda^{2}/\theta$

        - $\theta$ and $\sigma$ are both overdispersion parameters

        - Read Ver Hoef & Boveng (2007)

- Textbooks on statistical methodology in ecology recommend the use of the SQRT transformation


# Zeileis *et al.* 2008 - Regression models for count data in R

- The classic poisson regression is not suitable for count data due to over-dispersion and lots of zeroes.

- glm() - the Generalized linear model function

- use the models hurdle() and zeroinfl() in the pscl package

- Although quassi-poisson, negative binomial, and other Generalised Linear models are fine for overdispersion, they don’t deal well with excess zeroes

- *zero-augmented* models capture the zero counts in a second model component

- *Hurdle* models combine a left-truncated count component with a right-censored hurdle component

- *zero inflated* models are mixture models that combine a count component and a point mass at zero.


# Chang & Pocock 2000 Analysing data with clumping at zero An example demonstration

- Uses 2 approaches to analyse the relationship of multiple covariates to an outcome with lots of zeroes

    - Categorise the continuous outcome and fit a proportional odds model

    - Use a logistic regression to model the probability of a zero response and ordinary linear regression to model the non-zero continuous responses

- The data used here appears to be log-normal.

- glmer()

- {vegan} – pcnm()


# Incorporating the spatial influence of other cells

- Progressively increase the radius of influence and evaluate the fit
    with AIC

- Distance matrix with x,y coordinates of the trees, use this to make
    a pairwise distance matrix and add this as a covariate to determine
    their effect

- Use a Moran’s correlogram/Semivariogram to determine when the
    spatial autocorrelation ceases

    - Then use a principal coordinate analysis

- Remember to incorporate edge effects

# Zero-inflated models/Hurdle Models

- Allows zeroes to be removed from the model-

- Model in 2 parts, is there a difference between zeroes and ones, if
    a 1, what causes the difference.


# Links

https://stats.stackexchange.com/questions/59879/logistic-regression-anova-chi-square-test-vs-significance-of-coefficients-ano

https://stats.stackexchange.com/questions/46345/how-to-calculate-goodness-of-fit-in-glm-r

https://stats.stackexchange.com/questions/86351/interpretation-of-rs-output-for-binomial-regression

https://stats.stackexchange.com/questions/121593/using-the-same-variable-as-a-fixed-and-random-effect-in-mixed-effect-models

https://stats.stackexchange.com/questions/228800/crossed-vs-nested-random-effects-how-do-they-differ-and-how-are-they-specified

https://stats.idre.ucla.edu/other/mult-pkg/introduction-to-generalized-linear-mixed-models/

https://www.researchgate.net/post/How_to_interpret_interaction_in_a_glmer_model_in_R   - For graph

http://www.mypolyuweb.hk/~sjpolit//logisticregression.html  - Interpreting logistic regression output

http://bbolker.github.io/mixedmodels-misc/glmmFAQ.html#model-definition - Using `lme4`

https://strengejacke.wordpress.com/2017/08/27/marginal-effects-for-negative-binomial-mixed-effects-models-glmer-nb-and-glmmtmb-rstats/ - ggeffects for mixed models

https://www.researchgate.net/post/How_do_I_report_the_results_of_a_linear_mixed_models_analysis - pretentious guy explains reporting glmer

https://cran.r-project.org/web/packages/emmeans/vignettes/comparisons.html - Using Tukey's HSD test to compare groups in a mixed model and plotting pairwise comparisons as P-value charts
