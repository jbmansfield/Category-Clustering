#' ---
#' title: |
#'   | Category Clustering^[This document was `rmarkdown::render`'ed from an `R` script available at  \url{https://github.com/jbmansfield/Category-Clustering}]
#'   | *Supplementary Material 3: Statistical models*
#' author: ""
#' date: ""
#' output:
#'  pdf_document:
#'   fig_crop: true
#'   fig_caption: true
#'   latex_engine: xelatex
#'   toc: true
#'   toc_depth: 4
#'   number_section: true
#'   pandoc_args:
#'     - "--variable=lof"
#'     - "--variable=lot"
#'     - "--bibliography=/Users/bbickel/Documents/Bibliographies/repo/bbbib.bib"
#'     - "--csl=/Users/bbickel/Documents/Bibliographies/repo/unified-style-linguistics.csl"
#' header-includes: \usepackage[width=\textwidth]{caption}
#' ---

#+ setup, include=F
library(ggplot2)
library(dplyr)
library(tidyr)
library(knitr)
library(kableExtra)
library(brms)
library(HDInterval)
library(glmmTMB)
library(extraDistr)
knit_hooks$set(crop = hook_pdfcrop, pars = function(before, options, envir) {if(before) {par(family=my.font)} else NULL}) # trick for keeping par across chunks; define my.font below!
opts_chunk$set(fig.path='figures/',
           dev = 'cairo_pdf', dev.args=list(bg='transparent'), # or quartz_pdf (for lattice)
                fig.height = 7,
                fig.width = 14,
                message = F,
                warning = F,
                autodep=T,
                cache.comments=F,
                crop=T,
                pars=T,
                out.extra = ''
                )
# graphics setup:
my.font = 'Helvetica'
# ggplot
theme_set(theme_bw(base_size = 24) +
           theme(text = element_text(family = my.font),
                 axis.text = element_text(colour = 'black'),
                 plot.background = element_rect(fill = "transparent", colour = NA),
                 legend.text.align=0
           ))

options(width=180, knitr.kable.NA = '', knitr.table.format = "latex")

#' \clearpage
#'
#' Chintang prefix bigrams (Section 3.5)
#' ====================================
#'
#' Data
#' ----
#'
#' The variable names used here correspond to the terminology used in the main paper as follows:
#'
#'   - `ref.category`: "Co-prefix"
#'   - `subj.marker`: "Prefix identity"
#'   - `prime.subj`: "Persistence"
#'
#+ bigram-data, cache=T
subj.bigrams <- read.table("Chintang-prefixes/subj-bigrams.txt", sep = "\t", quote = "", header = TRUE)

#' We set the reference levels as follows:
#+ refs, cache=T
subj.bigrams <- within(subj.bigrams, prime.subj <- relevel(prime.subj, ref = "unknown"))
subj.bigrams <- within(subj.bigrams, subj.marker <- relevel(subj.marker, ref = "u"))
subj.bigrams <- within(subj.bigrams, ref.category <- relevel(ref.category, ref = "neg"))

#'
#' Main model
#' ----------
#'
#' We fit the model reported in the main text with the `brms` interface to the `stan` programming language [@Burkner2018Advanced]. We choose a skeptical Student-*t*(5,0,3) prior for the fixed effects, i.e. centered on 0 and therefore favoring the null hypothesis of no effect. For the random effects we use a half-*Cauchy*(0,1) prior [@McElreath2016Statistical].
#'
# /*
# data.frame(x = seq(-10, 10, 1)) %>% ggplot(aes(x = x)) +
#   stat_function(
#     fun = dstudent_t,
#     args = c(5, 0, 3),
#     geom = "area",
#     color = 'blue',
#     alpha = .2
#   ) +
#    stat_function(
#      fun = dcauchy,
#      args = c(0, 1),
#      geom = "area",
#      color = 'black',
#      alpha = .2
#      ) +
#        xlim(-5,5)
# */
#'
#'
#' We first fit a model that compares each prefix against an intercept fixed at log odds = 0, i.e. as deviations from  a uniform 0.5 probability for each prefix to occur either to the left or the right in the bigram.
#'
#+ bigrams-brms0, cache=T
bigrams.brm0 <- brm(left ~ 0 + subj.marker + ref.category + prime.subj +
                           subj.marker:ref.category + subj.marker:prime.subj +
                           (ref.category|lex.stem) + (ref.category|speaker.id),
                    prior = c(prior(student_t(5, 0, 3), class = b),
                              prior(cauchy(0,1), class = sd)),
                    sample_prior = T,
                    chains = 4,
                    cores = 4,
                    iter = 3000,
                    warmup = 1000,
                    data = subj.bigrams,
                    family = bernoulli(),
                    control = list(adapt_delta = 0.999)
)

#'
#' Visual convergence check:
#'
# plot(bigrams.brm0)
#'
#' Figure \ref{fig:bigrams-brms0-ppc} shows a posterior predictive check. $y_{rep}$ denotes the replicated predictions $\tilde{y}$ from the posterior distribution of the model parameters applied to the original data.
#'
#+ bigrams-brms0-ppc, fig.cap = 'Posterior predictive check of the Chintang prefix bigram model. Predicted values follow the data very well.', cache=T, fig.pos='h', fig.height=5, fig.scap='Posterior predictive check of the Chintang prefix bigram model'
pp_check(bigrams.brm0, n.samples = 100)
#'
#' \newpage
#' **Model summary (with mean estimates and equal-tailed intervals):**
#'
#+ bigrams-brms0-sum, cache=T, echo=F
summary(bigrams.brm0)
#'
#' All interactions include 0 in their 95% credibility interval, while the predicted effects of the clustering hypothesis are firmly away from 0: both prefixes are biased towards the same side dependent on whether the co-prefix is a negation or object marker. The variation by lexical stem is moderate (estimated mean intercept sd = `r round(summary(bigrams.brm0)$random$lex.stem[1,1],2)`; slope sd = `r round(summary(bigrams.brm0)$random$lex.stem[2,1],2)`). The variation by speaker is substantially larger (intercept sd = `r round(summary(bigrams.brm0)$random$speaker.id[1,1],2)`; slope sd = `r round(summary(bigrams.brm0)$random$speaker.id[2,1],2)`), in line with earlier findings [@Bickeletal2007Free].
#'
#' For the figure in the main text we compute the 95% highest density intervals. We do the same for the group level (random) effects. A synopsis is given in Table \ref{tab:bigrams-brms0-tab}.
#'
#+ bigrams-brms0-hdi, cache=F
bigrams.post0.df <- posterior_samples(bigrams.brm0, "^b|^sd")[,-c(8)] %>%
  gather(Term, Estimate)  %>%
  mutate(Term = factor(Term, levels = unique(Term),
                       labels = c("Prefix: u-",
                                  "Prefix: a-",
                                  "Co-Prefix: OBJ (vs. NEG)",
                                  "Persistence:\nprevious token left",
                                  "Persistence:\nprevious token right",
                                  "Prefix × Co-Prefix",
                                  "Prefix × Persistence",
                                  "Std. Dev. Intercept | Stem",
                                  "Std. Dev. Slope | Stem",
                                  "Std. Dev. Intercept | Speaker",
                                  "Std. Dev. Slope | Speaker")))
bigrams.post.sum0.df <- group_by(bigrams.post0.df, Term) %>%
  summarise(lo = hdi(Estimate, credMass = .95)[1],
            hi = hdi(Estimate, credMass = .95)[2],
            Median = median(Estimate),
            Mean = mean(Estimate),
            CI = paste0("[", round(hdi(Estimate, credMass = .95)[1], 2), ", ",
                   round(hdi(Estimate, credMass = .95)[2], 2), "]" ))

#+ bigrams-brms0-tab, cache = F, echo = F
kable(bigrams.post.sum0.df[,c(1,4:6)], col.names = c('Term', 'Median', 'Median', '95%CI (HDI)'),
      booktabs = T, digits = 2, caption.short = 'Estimated coefficients in the Chintang prefix  model',
      caption = 'Estimated coefficients in the Chintang prefix  model, with highest posterior density intervals', escape = T,
      linesep = c(rep('',4), '\\addlinespace', '', '\\addlinespace')
        ) %>% kable_styling(latex_options = "hold_position")


#+ bigrams-brms0-fig, cache=T, echo=F, fig.show = 'hide'
# quartz()
ggplot(bigrams.post.sum0.df[1:5,]) +
  geom_errorbarh(aes(xmin = lo, xmax = hi, y = reorder(Term, desc(Term))),
                 height = 0.3, size = 1.5) + # height = 0.1
  geom_point(aes(x = Median, y = Term), size = 3) +
  geom_vline(xintercept = 0, linetype = 'dashed', size = .3) +
  ylab("") + xlab("Estimated log odds")

#'
#'
#' Alternative model with intercept
#' --------------------------------
#'
#' In order to more specifically assess the difference between *a-* and *u-* we also fit a model with an intercept. We use sum-to-zero contrast coding, so the intercept is the grand mean, in line with the default treatment of intercepts in `brms`. We choose the same Student-*t*(5, 0, 3) prior as for the other parameters.
#'
# /*
# data.frame(x = seq(-10, 10, 1)) %>% ggplot(aes(x = x)) +
#   stat_function(
#     fun = dstudent_t,
#     args = c(5, 0, 3),
#     geom = "area",
#     color = 'blue',
#     alpha = .2
#   ) +
#    stat_function(
#      fun = dstudent_t,
#      args = c(3,0, 10),
#      geom = "area",
#      color = 'black',
#      alpha = .2
#      ) +
#        xlim(-5,5)
# */
#'
#+ bigrams-brms1, cache=T
contrasts(subj.bigrams$subj.marker) <- c(-.5, 5)
contrasts(subj.bigrams$ref.category) <- c(-.5, 5)
contrasts(subj.bigrams$prime.subj) <- matrix(c(-.5, 0, .5, 0, .5, -.5), ncol = 2)

bigrams.brm1 <- brm(left ~ 1 + subj.marker + ref.category + prime.subj +
                      subj.marker:ref.category + subj.marker:prime.subj +
                      (1 + ref.category|lex.stem) + (1 + ref.category|speaker.id),
                    prior = c(prior(student_t(5, 0, 3), class = Intercept),
                              prior(student_t(5, 0, 3), class = b),
                              prior(cauchy(0,1), class = sd)),
                    sample_prior = T,
                    chains = 4,
                    cores = 4,
                    iter = 3000,
                    warmup = 1000,
                    data = subj.bigrams,
                    family = bernoulli(),
                    control = list(adapt_delta = 0.999)
)

#' Visual convergence check:
#'
# plot(bigrams.brm1)
#'
#' **Model summary (with mean estimates and equal-tailed intervals):**
#'
#+ bigrams-brms1-sum, cache=T, echo=F
summary(bigrams.brm1)

#' Credible intervals based on highest posterior probability densities (instead of equal-tailed probabilities in the summary):
#'
#+ bigrams-brms1-hdi, cache=T
t(apply(posterior_samples(bigrams.brm1, '^b|^sd'), 2, function(x) {
  round(hdi(x, credMass = .95), 4) }))

ROPE <- posterior_samples(bigrams.brm1, 'b_subj.marker1', exact_match = T) %>%
  summarize(mean(b_subj.marker1 < 0.1 & b_subj.marker1 > -0.1))

#' The coefficient `subj.marker1` is the estimated difference in log odds of leftward placement between the *a-* and *u-* prefixes. Its 95\%CI includes zero, and in fact `r 100*round(ROPE,2)`% of the estimates are between -0.1 and .1. There is thus no reason to believe the prefixes behave differently relative to what is expected at the mean of the other parameters.
#'
#' Alternative frequentist model
#' -----------------------------
#'
#' For the sake of frequentist interests, we also fitted a model with `glmmTMB` [@Brooksetal2017glmmTMB]. However, the model converged only if limited to random intercepts, excluding random slopes.
#'
#+ bigram-tmb, cache=T

bigrams.TMB <- glmmTMB(left ~ 0 + subj.marker + ref.category + prime.subj +
                         ref.category:subj.marker + subj.marker:prime.subj +
                         (1|lex.stem) + (1|speaker.id),
                       data = subj.bigrams,
                       contrasts = list(subj.marker = "contr.treatment",
                                        ref.category = "contr.treatment",
                                        prime.subj  = "contr.treatment"),
                       family = "binomial")
summary(bigrams.TMB)

#'
#' Despite the difference in the random-effects structure, results are similar to what the Bayesian model showed: both prefixes show the same leftward probability and dependency on the co-prefix, all  significant at a rejection level of $\alpha = .05$.
#'
#'
#' Paradigmatic alignment (Section 4.3)
#' ===================================
#'
#' Data
#' ----
#'
#' The variable names used here correspond to the terminology used in the main text as follows:
#'
#' - `Atr`: "A"
#' - `Category`: "grammatical role (A or P)"
#' - `CLUST.INDEX` = "paradigmatic alignment index"
#' - `Stock`:  "language family"
#'
#+ align-data, cache=T
alignment.df <- read.table("clustering-typology/paradigm-partitions-probs_filt.txt",
                           header = T, sep = "\t")

#'
#'
#' Model
#' ----
#'
#' The paradigmatic index is 1 minus the cumulative probability of allocations (of markers to slots) that are aligned at least as much as the observed system, given all logically possible allocations. Since there exists necessarily at least one allocation (the observed one), this probability cannot be zero, and so the index cannot be 1. However, the cumulative probability can be 1 and the index 0, e.g. when the observed system allocates markers evenly to slots, so that *all* other allocations show more alignment (see Table 6 and 7 in the main text). This is the case in `r sum(alignment.df$CLUST.INDEX==0)` out of `r length(alignment.df$CLUST.INDEX==0)` cases.
#'
# /*
quartz()
alignment.df %>% mutate(Category = factor(Category, labels=c('A','P'))) %>%
ggplot(aes(x=CLUST.INDEX)) + geom_histogram(fill='blue') + facet_wrap(~Category) +
  xlab("Paradigmatic Alignment")
# */
#' Given these considerations, we fit a *mixture model* that allows both a zero vs. non-zero response and any response value in the interval (0,1).^[For a gentle introduction to mixture models in a Bayesian framework, see @McElreath2016Statistical; a helpful tutorial is provided by Matti Vuorre at \url{https://vuorre.netlify.com/post/2019/02/18/analyze-analog-scale-ratings-with-zero-one-inflated-beta-models}.] Such a model corresponds to what is also known as a *hurdle beta model* (or somewhat less accurately as a *zero-inflated beta model*). The first component is a logistic model with a response $\alpha$ (called `zi` in `brms`), the second a beta model with a response $\mu$. Both responses assume a logit link, so the regression coefficients represent the log odds ratio of probabilities (logistic) or proportions (beta) that they induce on the response with one unit change. Under the parametrization we use here [@Ferrarietal2004Beta;@Cribari-Netoetal2010Beta;@Figueroa-Zunigaetal2013Mixed;@Burkner2018Advanced], the beta model has one more parameter beyond the response $\mu$, often symbolized as $\phi$ (or $\nu$).^[The parameters $\mu$ and $\phi$ correspond to the more familiar parametrization in terms of “shape” parameters $\alpha$ and $\beta$ as follows: $\alpha = \mu\phi$ and $\beta=(1-\mu)\phi$.] The parameter is also known as a precision or overdispersion parameter and captures the inverse of the variance in the model fits.  Here, we model $\phi$ in the same way as the logistic and the beta component, albeit on the log scale (following the `brms` default).
#'
#' Like in the Chintang bigram model, we capture the main effects of interest by comparing each marker category, A and P, against an intercept fixed at log odds = 0, i.e. as parallel deviations from equal probabilities/proportions. We control for phylogenetic autocorrelation by including language family as a random intercept.
#'
#' We again assume a weakly informative, regularizing Student-*t*(5, 0, 3) prior for fixed effect predictors and a half-*Cauchy*(0,1) prior for the groups (random effects).
#'
# /*
# data.frame(x = seq(-10, 10, 1)) %>% ggplot(aes(x = x)) +
#   stat_function(
#     fun = dbeta,
#     args = c(2,2),
#     geom = "area",
#     color = 'blue',
#     alpha = .2
#   ) +
#    stat_function(
#      fun = dbeta,
#      args = c(1,1),
#      geom = "area",
#      color = 'black',
#      alpha = .2
#      ) +
#        xlim(-2,2)
# make_stancode(bf(CLUST.INDEX ~ 0 + Category + (1|Stock),
#        zi ~ 0 + Category + (1|Stock)),
#        data = alignment.df,
#        family = zero_inflated_beta(link_zi = 'logit')
# )
# */
#'
#'
#+ alignment-brm, cache = T
alignment.brm <- brm(bf(CLUST.INDEX ~ 0 + Category + (1|Stock), # mu: (0,1) responses
                        zi ~ 0 + Category + (1|Stock), # alpha: zero vs non-zero responses
                        phi ~ 0 + Category + (1|Stock) # phi: precision
                        ),
                    prior = c(prior(student_t(5, 0, 3), class = b),
                              prior(cauchy(0,1), class = sd)
                             ),
                    sample_prior = T,
                    chains = 4,
                    cores = 4,
                    iter = 5000,
                    warmup = 2500,
                    data = alignment.df,
                    family = zero_inflated_beta(),
                    control = list(adapt_delta = 0.99)
)

#'
#'
#' Visual convergence check:
# plot(alignment.brm)
#'
#' Figure \ref{fig:alignment.brm-ppc} shows a posterior predictive check.
#'
#+ alignment.brm-ppc, fig.cap = 'Posterior predictive check of paradigmatic alignment model. Predictions are close to observations', cache=T, fig.pos='h', fig.height=5, fig.scap='Posterior predictive check of paradigmatic alignment model'
pp_check(alignment.brm, n.samples = 100)

#'
#'
#'
#'
#' **Model summary (with mean estimates and equal-tailed intervals):**
#+ alignment.brm-sum, echo=F, cache=T
summary(alignment.brm)

#'
#' To help interpret the estimates we convert them to the response scale by taking the inverse logit for the main coefficients and the inverse log for the precision parameter $\phi$. Results are displayed in Table \ref{tab:alignment.brm-tab-print}.
#' \captionsetup{margin = 2cm}
#+ alignment.brm-tab, cache = T
alignment.post.probs.df <- posterior_samples(alignment.brm, "^b|^sd") %>%
  mutate_at(vars(contains("phi")), exp) %>% # inverse of the log link for phi
  mutate_at(vars(-contains("phi")), plogis) %>% # inverse of the logit link
  gather(Parameter, Estimate) %>% group_by(Parameter) %>%
  summarise(Median = median(Estimate),
            Mean = mean(Estimate),
            lo = hdi(Estimate, credMass = .95)[1],
            hi = hdi(Estimate, credMass = .95)[2])

#+ alignment.brm-tab-print, echo = F
alignment.post.probs.df %>% mutate(CI = paste0("[", round(lo, 2), ", ", round(hi, 2), "]" )) %>%
  mutate(Term = c(rep(c("A","P"),3), rep("Family", 3)),
         Parameter = c(rep('$\\mu$',2), rep('$\\phi$', 2), rep('$\\alpha$ (\\texttt{zi})',2),
                       'sd($\\mu$)', 'sd($\\phi$)', 'sd($\\alpha$)')) %>%
  select(Term, Parameter, Median, Mean, CI) %>%
kable(., digits = 2, col.names=c('Term', 'Parameter', 'Median', 'Mean', '95\\%CI (HDI)'), caption.short = 'Estimated effects in the paradigmatic alignment model',
      caption = 'Estimated effects on the response scale in the paradigmatic alignment model, with highest posterior density intervals', booktabs = T, linesep = '', escape = F) %>%
  kable_styling(latex_options = 'hold_position')

#'
# /*
quartz()
ggplot(droplevels(subset(alignment.post.probs.df, !grepl('phi', Parameter)))) +
  geom_errorbarh(aes(xmin = lo, xmax = hi, y = rev(Parameter)), height = 0.3, size = 1.2) + # .1
  geom_point(aes(x = Median, y = rev(Parameter)), size = 3) +
  geom_vline(xintercept = 0.5, linetype = 'dashed', size = .3, color = 'blue') +
  scale_y_discrete(labels = rev(c(expression(paste('(0,1): ', mu(A))),
                                  expression(paste('(0,1): ', mu(P))),
                                  expression(paste("0 vs (0, 1): ", alpha(A))),
                                  expression(paste("0 vs (0, 1): ", alpha(P))),
                             expression(paste('SD(',mu,'|Family)')),
                             expression(paste('SD(',alpha,'|Family)'))))) +
  ylab("") + xlab("Estimated Paradigmatic Alignment")
# */
#'
#' The estimated alignment index values are high for both A and P (mean and median $\mu$ = `r round(alignment.post.probs.df[1,2],2)`), with 95%CI far away from .5 proportions (i.e. from 0 log odds). The precision parameters $\phi$ are relatively high (median `r round(alignment.post.probs.df[3,2],2)` for A and `r round(alignment.post.probs.df[4,2],2)` for P), suggesting modest variance (cf. Figure \ref{fig:beta-ex} for illustrations of what effect different $\phi$ values have on the shape of the estimates). At the same time, the estimated probabilities of zero values, i.e. minimal alignment, are exceedingly small (median $\alpha$ = `r round(alignment.post.probs.df[5,2],3)` for A and $\alpha$ = `r round(alignment.post.probs.df[6,2],3)` for P). They are notably lower than the marginal estimates in Figure 7 in the main text (where there are `r mean(alignment.df$CLUST.INDEX[alignment.df$Category %in% 'Atr'] == 0)` zeros in A and `r round(mean(alignment.df$CLUST.INDEX[alignment.df$Category %in% 'P'] == 0),3)` in P). This is due to the fact that the model controls for the historical relationships between languages in the random effects. Indeed, many languages with 0 alignment  come from a single family (Berber), so the values are historically dependent on each other (Table \ref{tab:alignment-zero-tab}). In line with this, we note a high standard deviation estimate of the random intercept for $\alpha$ (median `r round(alignment.post.probs.df[9,2],2)`). The variation between families also has a substantial impact on the random intercepts of the beta parameters ($\mu$ and $\phi$ in Table \ref{tab:alignment.brm-tab-print}). As discussed in Section 4.4 of the main text, this is largely due to the overrepresentation of Algonquian and Kiranti languages.
#' \captionsetup{margin = 1cm}
#+ alignment-zero-tab, echo=F
subset(alignment.df, CLUST.INDEX == 0, c('Stock', 'Language', 'Category', 'PARTITION')) %>%
  mutate(Category=factor(Category, labels=c('A','P'))) %>% arrange(Stock) %>%
  kable(.,booktabs=T, caption = 'Languages with no paradigmatic alignment (i.e. index value 0)', col.names = c('Family', 'Language', 'Category', 'Partition'), linesep='', row.names=F, caption.short = 'Languages with no paradigmatic alignment') %>% kable_styling(latex_options = 'hold_position')
# /*
as.data.frame(xtabs(~Stock, droplevels(subset(alignment.df,CLUST.INDEX < 0.5)))) %>% arrange(-Freq)
# */
#'
#+ beta-ex, echo=F, fig.cap = 'Exemplifications of various beta distribution parameters to help interpretation', cache=F,fig.height=5, fig.pos = 'h'
data.frame(x = seq(0, 1, .001)) %>% ggplot(aes(x = x)) +
  stat_function(
    fun = dprop,
    args = list(size=2, mean=.86),
    geom = "line", size = 1,
    aes(color = "A")
  ) +
  stat_function(
    fun = dprop,
    args = list(size=3, mean=.86),
    geom = "line", size = 1,
    aes(color = "B")
  ) +
  stat_function(
    fun = dprop,
    args = list(size=4, mean=.86),
    geom = "line", size = 1,
    aes(color = "C")
  ) +
  stat_function(
    fun = dprop,
    args = list(size=6, mean=.86),
    geom = "line", size = 1,
    aes(color = "D")
  ) +
  stat_function(
    fun = dprop,
    args = list(size=10, mean=.86),
    geom = "line", size = 1,
    aes(color = "E")
  ) +
  stat_function(
    fun = dprop,
    args = list(size=20, mean=.86),
    geom = "line", size = 1,
    aes(color = "F")
  ) +
  scale_color_manual(values = c('black', 'red', 'blue', 'purple', 'darkgreen', 'orange'),
                      labels = c(expression(italic(Beta)(mu==.86, phi==2)),
                                 expression(italic(Beta)(mu==.86, phi==3)),
                                 expression(italic(Beta)(mu==.86, phi==4)),
                                 expression(italic(Beta)(mu==.86, phi==6)),
                                 expression(italic(Beta)(mu==.86, phi==10)),
                                 expression(italic(Beta)(mu==.86, phi==20))
                                 )) +
  labs(x='Proportion',y="Density",color="") +
  theme(legend.position="top", legend.key.width = unit(1.5, "cm"), legend.text.align=0) +
  guides(color = guide_legend(override.aes = list(size = 1.5), nrow = 2, byrow = T))

#' \clearpage
#'
#' Alternative frequentist model
#' -----------------------------
#'
#' We are not aware of a frequentist implementation of zero-inflated beta regressions (as of April 2019), and so we instead smooth the data, avoiding values at zero. Note that $\phi$ is called here an overdispersion parameter and is estimated for the entire model, not allowing variation by predictors.
#'
#+ alignment-rmb, cache=T
alignment.df$CLUST.INDEX.smoothed <- (alignment.df$CLUST.INDEX *
                                     (length(alignment.df$CLUST.INDEX) - 1) + .5) /
                                     length(alignment.df$CLUST.INDEX)

alignment.TMB <- glmmTMB(CLUST.INDEX.smoothed ~ 0 + Category + (1|Stock),
                         data = alignment.df,
                         family = beta_family())
summary(alignment.TMB)

#'
#' Estimates are much lower here since the smoothing moves 0 values to the probability mass near zero. At any rate, the estimated responses of the alignment index (A: logit$^{-1}$(`r round(fixef(alignment.TMB)$cond[1],2)`) = `r round(plogis(fixef(alignment.TMB)$cond[1]),2)`; P: logit$^{-1}$(`r round(fixef(alignment.TMB)$cond[2],2)`) = `r round(plogis(fixef(alignment.TMB)$cond[2]),2)`) are still significantly higher than logit$^{-1}$(0) = .5 at a rejection level of $\alpha = .001$.
#'
#'
#' Featural coherence (Section 4.5)
#' ================================
#'
#' Data
#' ----
#'
#' In this data, `crmv.crct` refers to the bias-corrected version of Cramér's $V$ measure of association [@Bergsma2013A-bias-correction]. The bias correction leads to an `NA` value in the case of the Zuni paradigm, which has just one marker for A and P in our database, and these go in two different slots. We remove this datapoint:
#'
#+ coherence-data, cache=T
coherence.df <- read.table("clustering-typology/distinctive-alignment-crmv.txt", sep = "\t",
                           header = T, stringsAsFactors = F) %>%
                filter(!is.na(crmv.crct))

#' We recode the language families, so that we can zoom in on Kiranti and Algonquian. Note that all Sino-Tibetan languages in the database are in fact from the Kiranti clade, and all Algic languages from the Algonquian clade; we rename the family accordingly.

#+ coherence-recode, cache=T
coherence.df$Stock.grouped[
  coherence.df$Stock %in% names(table(coherence.df$Stock))[table(coherence.df$Stock) < 10]] <- "Other"
coherence.df$Stock.grouped[coherence.df$Stock.grouped %in% "Sino-Tibetan"] <- "Kiranti"
coherence.df$Stock.grouped[coherence.df$Stock.grouped %in% "Algic"] <- "Algonquian"
coherence.df$Stock.grouped <- factor(coherence.df$Stock.grouped)
coherence.df$Stock.grouped <- factor(coherence.df$Stock.grouped, c("Algonquian", "Kiranti", "Other" ))

#'
# /*
sapply(list(m1=matrix(c(1,0,0,5,4,0), nrow=2),
            m2=matrix(c(1,0,4,3,0,2), nrow=2),
            m4=matrix(c(1,1,2,2,2,2), nrow=2)),
       rcompanion::cramerV, bias.correct = T)
quartz()
ggplot(coherence.df, aes(x=crmv.crct)) + geom_histogram(fill='blue') +
  xlab("Featural Coherence")

# */
#' Model
#' ----
#'
#' The data contain both 0s and 1s (`r sum(coherence.df$crmv.crct==0)` and `r sum(coherence.df$crmv.crct==1)` out of `r length(coherence.df$crmv.crct)`, respectively). Therefore, we fit again a mixture model, but we now also model 1s, not only 0s. In `brms` this is parametrized via `zoi`, the probability $\alpha$ of 0 *or* 1, and `coi`, the conditional probability $\gamma$ that 1 occurs rather than 0. Since in this model, there is no predictor term (but only a control for random effects), we assume a single estimate of the precision parameter $\phi$ for all observations, with no further transformation.
#'
#' We choose the same priors as in the alignment model; for $\phi$ we assume the *Gamma*(.01,.01) prior that `brms` takes as a default and which has been found to perform well [@Figueroa-Zunigaetal2013Mixed]. We also tried a slightly less constrained *Gamma*(.1,.1) prior, but found no appreciable difference in the results.^[Note that $\phi$ is naturally constrained to be positive.]
#'
# /*
# data.frame(x = seq(0, 20, 1)) %>% ggplot(aes(x = x)) +
#   stat_function(
#     fun = dgamma,
#     args = list(shape=1, rate=1),
#     geom = "line", size = 1,
#     aes(color = "A")
#   ) +
#   stat_function(
#     fun = dgamma,
#     args = list(shape=.1, rate=.1),
#     geom = "line", size = 1,
#     aes(color = "B")
#   ) +
#   stat_function(
#     fun = dgamma,
#     args = list(shape=.01, rate=.01),
#     geom = "line", size = 1,
#     aes(color = "C")
#   )
#
# get_prior(bf(crmv.crct ~  0 + intercept + (1|Stock),
# zoi ~  0 + intercept + (1|Stock),
# coi ~  0 + intercept + (1|Stock)),
# family = zero_one_inflated_beta(),
# data = coherence.df)
#
# */
#'
#+ coherence-brm, cache=T

coherence.brm <- brm(bf(crmv.crct ~  0 + intercept + (1|Stock), # mu: (0,1) response
                        zoi ~  0 + intercept + (1|Stock), # alpha: Pr(0 or 1) response
                        coi ~  0 + intercept + (1|Stock)), # gamma: Pr(1 | 0 or 1) response
                        prior = c(prior(student_t(5, 0, 3), class = b),
                                  prior(cauchy(0,1), class = sd),
                                  prior_string("gamma(0.01,0.01)", class = "phi")),
                        sample_prior = T,
                        chains = 4,
                        cores = 4,
                        iter = 7000,
                        warmup = 5000,
                        data = coherence.df,
                        family = zero_one_inflated_beta(),
                        control = list(adapt_delta = 0.9999999999999, max_treedepth = 35)
)

#'
#'
#' \captionsetup{margin = 0cm}
#' Visual convergence check:
# plot(coherence.brm)
#'
#' Figure \ref{fig:coherence.brm-ppc} shows a posterior predictive check.
#'
#+ coherence.brm-ppc, fig.cap = 'Posterior predictive check of the coherence model. Predictions are relatively close to observations.', cache=T, fig.pos='h', fig.height=5, fig.scap='Posterior predictive check of the coherence model.'
pp_check(coherence.brm, n.samples = 100)

#'
#' \newpage
#'
#' **Model summary (with mean estimates and equal-tailed intervals):**
#+ coherence.brm-sum, echo=F
summary(coherence.brm)

#'
#' For interpretation we again convert the estimates to the response scale for easier interpretation of the effects. Results are displayed in Table \ref{tab:alignment.brm-tab-print}.
#'
#' \captionsetup{margin = 1cm}
#+ coherence.brm-tab, cache = F
coherence.post.probs.df <- posterior_samples(coherence.brm, "^b|^sd") %>%
  mutate_all(plogis) %>%  # inverse of the logit link
  gather(Parameter, Estimate) %>% group_by(Parameter) %>%
  summarise(Median = median(Estimate),
            Mean = mean(Estimate),
            lo = hdi(Estimate, credMass = .95)[1],
            hi = hdi(Estimate, credMass = .95)[2],
            CI = paste0("[", round(lo, 2), ", ", round(hi, 2), "]" ))

#+ coherence.brm-tab-print, echo = F, include = T
coherence.post.probs.df$Parameter <- factor(coherence.post.probs.df$Parameter,
                                           levels = c(
                                             "b_intercept",
                                             "b_zoi_intercept",
                                             "b_coi_intercept",
                                             "sd_Stock__Intercept",
                                             "sd_Stock__zoi_Intercept",
                                             "sd_Stock__coi_Intercept"), ordered = T)
coherence.post.probs.df[order(coherence.post.probs.df$Parameter),] %>%
  mutate(Term = c(rep("Intercept",3), rep("Family", 3)),
         Parameter = c("$\\mu$ (i.e. $0 > V < 1$)",
                       "$\\alpha$ (\\texttt{zoi}, i.e. $P(V \\in \\{0,1\\}$))",
                       "$\\gamma$ (\\texttt{coi}, i.e. $P(V = 1 | V \\in \\{0,1\\}$))",
                       "sd($\\mu$)",
                       "sd($\\alpha$)",
                       "sd($\\gamma$)")
         ) %>%
select(Term, Parameter, Median, Mean, CI) %>%
kable(., digits = 2, col.names=c('Term', 'Parameter', 'Median', 'Mean', '95\\%CI (HDI)'),
      caption.short='Estimated effects in the featural coherence model',
      caption = 'Estimated effects on the response scale in the featural coherence model, with  highest posterior density intervals.',
      booktabs = T, linesep = '', escape = F) %>%
  kable_styling(latex_options = 'hold_position')

# /*
quartz()

ggplot(coherence.post.probs.df[c(2,3,1),]) +
  geom_errorbarh(aes(xmin = lo, xmax = hi, y = reorder(Parameter, desc(Parameter))),
                 height = 0.3, size = 1.2) +
  geom_point(aes(x = Median, y = reorder(Parameter, desc(Parameter))), size = 3) +
  # geom_vline(xintercept = 0.5, linetype = 'dashed', size = .3, color = 'blue') +
  scale_y_discrete(labels = rev(c(expression({0>V}*phantom()<1),
                                  # expression(italic(P)(V==1)),
                                  # expression(italic(P)(V==0))
                                  expression(italic(P)(V %in% {group("{",list(0,1),"}")})),
                                  expression(italic(P)(V==1~"|"~V %in% {group("{",list(0,1),"}")}))
                                  # expression(paste("1: ", gamma, alpha)),
                                  # expression(paste("0: ", (1 - gamma), alpha)),
                                  # expression(paste('SD(',mu,' | Family)')),
                                  # expression(paste('SD(',alpha,' | Family)')),
                                  # expression(paste('SD(',gamma,' | Family)'))
                                  ))) +
  ylab("") + xlab(expression(paste("Estimated Featural Coherence", italic(V)))) +
  xlim(0,1)
# */
#'
#' The results in Table \ref{tab:coherence.brm-tab-print} evidence a value of Cramér's $V$ considerably above 0.5, with a 95%CI of $\mu$ = `r coherence.post.probs.df[2,6]` and relatively high precision (mean $\phi$ = `r round(summary(coherence.brm)$spec_pars[1],2)`; see Figure \ref{fig:beta-ex}).  Furthermore, of the estimated probability that values are 0 or 1 as opposed to inside the (0,1) interval (95%CI of $\alpha$ = `r coherence.post.probs.df[3,6]`), almost all are estimated to be 1 ($\gamma$ = `r coherence.post.probs.df[1,6]`). This suggests that the apparent high count of 0s in Figure 8 of the main text is an artefact of historical dependencies between languages of the same family. Once this is controled for, the  model estimates are lower than the observed counts. The variation among families is indeed substantial for all responses (Table \ref{tab:coherence.brm-tab-print}) and the effects on 0s chiefly stem from multiple representatives of only two families, Algonquian and Kiranti (Table \ref{tab:coherence-zero-tab}).
#' \captionsetup{margin = 0cm}
#+ coherence-zero-tab, echo=F
subset(coherence.df, crmv.crct == 0, c('Language', 'Stock', 'Atr', 'P')) %>%
  mutate(Family = ifelse(Stock %in% 'Algic', 'Algonquian', # they are all Algonquian
                         ifelse(Stock %in% 'Sino-Tibetan', 'Kiranti', Stock))) %>%  # all Kiranti
  arrange(Family) %>% select(Family, Language:P, -Stock) %>%
  kable(.,booktabs=T, caption = "Languages with no featural coherence (i.e. Cramér's $V = 0$). The ‘Distribution’ columns indicate the number of markers of the respective agreement category that occur in a given position.", col.names = c('Family', 'Language', 'Distribution of A', 'Distribution of P'), linesep='', row.names=F, caption.short = 'Languages with no featural coherence', longtable =T
        ) %>% kable_styling(latex_options = c("repeat_header", "hold_position"),
                            repeat_header_text = "Languages with no featural coherence (\\emph{cont'd})",
                            repeat_header_method = 'replace')
# /*
as.data.frame(xtabs(~Stock, droplevels(subset(coherence.df,crmv.crct < 0.5)))) %>% arrange(-Freq)
# */
#'
#'
#' Alternative frequentist model
#' -----------------------------
#'
#' In the absence of mixture models for beta regressions, we again smooth the response variable at the extremes before fitting the model:
#'
#+ coherence-tmb, cache=T
coherence.df$corr.crmv <- (coherence.df$crmv.crct * (length(coherence.df$crmv.crct) - 1) + .5) /
  length(coherence.df$crmv.crct)

coherence.TMB <- glmmTMB(corr.crmv ~ (1|Stock),
                         data = coherence.df,
                         family = beta_family())

summary(coherence.TMB)

#'
#' The estimated effect of logit$^{-1}$(`r round(fixef(coherence.TMB)$cond,2)`) = `r round(plogis(fixef(coherence.TMB)$cond),2)` is significantly above logit$^{-1}$(0) = .5 at a rejection level of $\alpha = .001$. It is slightly higher than in the Bayesian mixture model (`r round(plogis(fixef(coherence.TMB)$cond),2)` vs. a median of `r round(coherence.post.probs.df[2,2],2)`). This reflects the fact that in the `glmmTMB` model, the single response also absorbs the high count of 1s, which were modeled separately in the mixture model.
#'
#' \clearpage
#'
#'
#' References
#' ==========
#'
