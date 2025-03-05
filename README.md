# Wine Quality Prediction using Bayesian Methods
Ever wondered what makes a great white wine? In this project, we use Bayesian statistics to predict wine quality based on physicochemical properties like acidity, sugar content, and alcohol levels. Using Bayesian logistic regression and MCMC sampling (JAGS), we classify wines into two categories:
  + High quality (Rating ≥ 7)
  + Not high quality (Rating ≤ 6)

By comparing non-informative vs. informative priors, we also explore how expert knowledge influences our predictions

## Dataset
- The dataset originates from Cortez et al. (2009) "Modeling Wine Preferences by Data Mining from Physicochemical Properties."
- It includes 5,000+ observations and 12 variables:
  + 11 Independent Variables: Acidity, pH, Alcohol, Sulphates, Sugar, etc
  + 1 Dependent Variable: Wine quality rating (converted to a binary classification)

## Methodology
### Data Preprocessing
- Convert wine quality scores into binary categories
- Identified outliers and applied Bayesian guess parameters
- Split dataset: 70% training, 30% testing
### Bayesian Modeling - Logistic Regression
- Likelihood Function: Logistic model for binary classification
- Two Priors Tested:
  + Non-informative (Neutral): No prior knowledge with high large variance, which reflects the absence of strong prior beliefs
  + Informative (Expert-Driven): Based on prior research findings about White wine input importances (in %)
- Guess Parameter: Helps the model handle unexpected outliers
### MCMC Sampling & Posterior Inference
- Used JAGS + RStudio for Bayesian inference
- MCMC Settings: Optimized through extensive testing, adjusting parameters such as burn-in steps, adaptation period, thinning intervals, saved steps, number of chains, and total iterations to ensure efficiency and accuracy
- Diagnostics: Checked for convergence, autocorrelation, and effective sample size (ESS)
### Confusion Matrix 
- We used Confusion Matrices to evaluate classification performance through accuracy, precision, recall, F1-score.

  ## Results

  ## References
Cortez, P. et al (2009). Modeling wine preferences by data mining from physicochemical properties. Semantic Scholar. https://www.semanticscholar.org/paper/Modeling-wine-preferences-by-data-mining-from-Cortez-Cerdeira/bf15a0ccc14ac1deb5cea570c870389c16be019c

Cortez, P. et al (2009). Wine Quality [Dataset]. Machine Learning Repository. https://archive.ics.uci.edu/dataset/186/wine+quality

Demirhan, H., Demirhan, K. A Bayesian approach for the estimation of probability distributions under finite sample space. Stat Papers 57, 589–603 (2016). https://doi.org/10.1007/s00362-015-0669-z

Demirhan, H. (n.d.) Beta Distribution Specified by Mean and Concentration. [Shinyapps.io]. https://rmitsam.shinyapps.io/beta_3/

H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2016.

Kruschke, J. K. (2015). Doing Bayesian Data Analysis, Second Edition: A Tutorial with R, JAGS, and Stan. Academic Press / Elsevier. 
