# Tutorial on generalized linear mixed models in R, 
# with an emphasis on proportional count data

# Evan Eskew
# Original version: 05.15.15
# Updated: 05.09.16

#==============================================================================


# Load packages needed for data analysis
library(lme4)
library(AICcmodavg)
library(rethinking)

# Define a logistic() function that is needed for interpretation of 
# parameter output from binomial family generalized linear models 
# (e.g., if the intercept estimate from a binomial GLM is 2, then the 
# maximum likelihood estimate (MLE) for the probability of success will be
# logistic(2) = 0.88)
logistic <- function(x, print = TRUE) {
  result <- exp(x)/(1 + exp(x))
  return(result)
}

#==============================================================================


# Data import and cleaning


# Import the data into R
d <- read.csv("glmm_tutorial_data.csv")

# Examine the data
d
str(d)

# Convert Study to a factor (as opposed to an integer variable)
d$Study <- as.factor(d$Study)
str(d)

# Create a new Prev column that computes the observed prevalence (by row) 
# from the data and a Negative column that equals Num_Tested - Positive 
# for each row
d$Prev <- d$Positive/d$Num_Tested
d$Negative <- d$Num_Tested - d$Positive
d

# Create dataframes subsetted by Habitat to be used later
d_Ag <- d[d$Habitat == "Ag", ]
d_Forest <- d[d$Habitat == "Forest", ]

#==============================================================================


# Generalized linear models

# Begin modeling procedure with binomial generalized linear models (GLMs). 
# We're interested in the effects of Habitat and Season on prevalence, so we 
# fit four models: a null model without any fixed effects (aka an intercept 
# only model), a model including the fixed effect of Habitat, a model including 
# the fixed effect of Season, and an additive model considering the 
# fixed effects of both Habitat and Season.
m1 <- glm(cbind(Positive, Negative) ~ 1, family = binomial, data = d)
m2 <- glm(cbind(Positive, Negative) ~ Habitat, family = binomial, data = d)
m3 <- glm(cbind(Positive, Negative) ~ Season, family = binomial, data = d)
m4 <- glm(cbind(Positive, Negative) ~ Habitat + Season, 
          family = binomial, data = d)


# Conduct model comparison based on AICc (note that the correct number of 
# observations here is the total number of individuals tested)
aictab(list(m1, m2, m3, m4), 
       modnames = c("Intercept Only", "Habitat",
                    "Season", "Habitat + Season"), 
       nobs = sum(d$Num_Tested))


# What are these models giving us?
summary(m2)

# And how do we extract just the parameter information therein? Note this will
# only give you the maximum likelihood estimate for the parameter. It does 
# not account for uncertainty in the parameter estimate.
coef(m2)

# Or you could use coeftab() for comparison of parameter estimates among models
coeftab(m1, m2, m3, m4, nobs = FALSE)


# How does this translate into expected rate of infection (aka the MLE 
# for probability of success)? Let's use the parameter values from model m2.

# For Agricultural habitats
logistic(coef(m2)[1])

# For Forest habitats (you need the intercept value PLUS the effect of Forest)
logistic(coef(m2)[1] + coef(m2)[2])

# How does that compare to the overall prevalence observations in 
# Agricultural vs. Forest Habitat?
sum(d_Ag$Positive)/sum(d_Ag$Num_Tested)
sum(d_Forest$Positive)/sum(d_Forest$Num_Tested)

# Essentially, the binomial GLM with Habitat as a main effect (m2) has 
# fit the OVERALL observed prevalence in Agricultural vs. Forest habitats
# extremely well

#==============================================================================


# Optional (be careful with interpretation of models like this!)


# Alternatively, you can perform what's called a means parameterization 
# Let's do a means parameterization of model m2, which we'll call m2.alt
m2.alt = glm(cbind(Positive, Negative) ~ Habitat - 1, 
             family = binomial, data = d)


# Compare models m2 and m2.alt 
aictab(list(m2, m2.alt), modnames = c("m2", "m2.alt"), 
       nobs = sum(d$Num_Tested))

coef(m2)
coef(m2.alt)

# m2.alt contains the exact same information as m2, except the 
# "HabitatForest" parameter is fit as the actual mean for the Forest 
# category (rather than being the DIFFERENCE between Forest and 
# Agricultural habitats as in m2)


# This sort of parameterization can make model interpretation easier in some
# cases. For example, getting expected rate of infection now requires only 
# one parameter input for each Habitat type.

# For Agricultural habitats
logistic(coef(m2.alt)[1])

# For Forest habitats
logistic(coef(m2.alt)[2])

#==============================================================================


# Back to the main story. Let's examine the data a bit more in depth to 
# see what sorts of variation might be lurking about.


# Plot number of animals tested vs. study
plot(Num_Tested ~ as.numeric(Study), data = d, pch = 19, las = 1, 
     col = ifelse(d$Habitat == "Ag", "brown", "darkgreen"), 
     xlab = "Study", ylab = "Number of Animals Tested")

legend("topright", pch = 19, col = c("brown", "darkgreen"), 
       legend = unique(d$Habitat), bty = "n", cex = 1.2)

# Studies have drastically different numbers of animals tested. You could 
# have also seen this by just looking at the raw data.


# Plot observed prevalence vs. number of animals tested
plot(Prev ~ Num_Tested, data = d, pch = 19, las = 1,
     col = ifelse(d$Habitat == "Ag", "brown", "darkgreen"),
     xlab = "Number of Animals Tested", ylab = "Observed Prevalence")

legend("topleft", pch = 19, col = c("brown", "darkgreen"), 
       legend = unique(d$Habitat), bty = "n", cex = 1.2)

# Seems to be a relationhip here too...


# Critical observation incoming!!! 

# What do MOST of the studies tell us about the effect of Forest?
plot(Prev ~ as.numeric(Study), data = d, pch = 19, las = 1, xaxt = "n",
     col = ifelse(d$Habitat == "Ag", "brown", "darkgreen"), 
     ylim = c(0, 0.8), xlab = "Study", ylab = "Observed Prevalence")
axis(1, at = 1:11, labels = 1:11)

legend("topright", pch = 19, col = c("brown", "darkgreen"), 
       legend = unique(d$Habitat), bty = "n", cex = 1.2)

# So 8/11 studies found LOWER prevalence in forested habitats

# Let's remind ourselves of the predicted mean prevalence by Habitat from 
# model m2
abline(h = logistic(coef(m2)[1]), col = "brown", lwd = 3)
abline(h = logistic(coef(m2)[1] + coef(m2)[2]), col = "darkgreen", lwd = 3)

# So again, just as a reminder, our model m2 suggested that prevalence 
# should be HIGHER in forested habitats, which is not what we see in the 
# data when we look at the results on a per study basis

#==============================================================================


# Generalized linear mixed models


# So perhaps an approach controlling for Study would be appropriate. 
# Enter generalized linear mixed models (GLMMs). Here our random effects 
# (also called varying effects) will be for Study. Each model will 
# essentially fit a different intercept value for every study in 
# our dataset (11 in all). Below we have the same model set as before except
# with the random effect of Study included in every one. Note we now have 
# to fit the models with a different function (glmer), but the general 
# format is similar.
m5 <- glmer(cbind(Positive, Negative) ~ (1|Study), 
            family = binomial, data = d)
m6 <- glmer(cbind(Positive, Negative) ~ Habitat + (1|Study), 
            family = binomial, data = d)
m7 <- glmer(cbind(Positive, Negative) ~ Season + (1|Study), 
            family = binomial, data = d)
m8 <- glmer(cbind(Positive, Negative) ~ Habitat + Season + (1|Study), 
            family = binomial, data = d)


# Model comparison again
aictab(list(m5, m6, m7, m8), 
       modnames = c("Intercept Only (w Ran Effects)", 
                    "Habitat (w Ran Effects)",
                    "Season (w Ran Effects)",
                    "Habitat + Season (w Ran Effects)"), 
       nobs = sum(d$Num_Tested))

# Our inferences about the best model have completely changed due to 
# the inclusion of the random effect!


# Comparison of the parameter output from model m2 and m6. Note the 
# value of the "HabitatForest" coefficient.
summary(m2)
summary(m6)

# So not only have our relative model rankings changed, but the values of 
# the fixed effects coefficients have changed. We now have the complete
# opposite inference about the effect of Habitat! m2 suggests that Forest
# increases prevalence, whereas m6 suggests that being in Forest decreases
# prevalence.


# And for completeness, the actual random effects estimates for a model 
# (so in this case the actual study-level effect on prevalence) can be 
# viewed like so:
ranef(m6)

#==============================================================================


# Plotting model results


# Plot estimated prevalence by Habitat for m2 to remind ourselves what that
# looks like
plot(Prev ~ as.numeric(Study), data = d, pch = 19, las = 1, xaxt = "n",
     col = ifelse(d$Habitat == "Ag", "brown", "darkgreen"), 
     ylim = c(0, 0.8), xlab = "Study", ylab = "Observed Prevalence")
axis(1, at = 1:11, labels = 1:11)

legend("topright", pch = 19, col = c("brown", "darkgreen"), 
       legend = unique(d$Habitat), bty = "n", cex = 1.2)

abline(h = logistic(coef(m2)[1] + coef(m2)[2]*0), col = "brown", lwd = 3)
abline(h = logistic(coef(m2)[1] + coef(m2)[2]*1), col = "darkgreen", lwd = 3)


# Plot estimated prevalence by Habitat for m6. Note we are NOT explicitly
# plotting the effect of the random effects here. We're only looking at 
# what the model says an "average" study would find in Agricultural and 
# Forest habitat. But crucially our estimates of these effects have changed 
# in m6 (relative to m2) because we included the random effect of Study in 
# this model.
plot(Prev ~ as.numeric(Study), data = d, pch = 19, las = 1, xaxt = "n",
     col = ifelse(d$Habitat == "Ag", "brown", "darkgreen"), 
     ylim = c(0, 0.8), xlab = "Study", ylab = "Observed Prevalence")
axis(1, at = 1:11, labels = 1:11)

legend("topright", pch = 19, col = c("brown", "darkgreen"), 
       legend = unique(d$Habitat), bty = "n", cex = 1.2)

abline(h = logistic(fixef(m6)[1] + fixef(m6)[2]*0), col = "brown", lwd = 3)
abline(h = logistic(fixef(m6)[1] + fixef(m6)[2]*1), col = "darkgreen", lwd = 3)


# Now plot estimated prevalence by Habitat AND account for the random effect 
# of study (using estimates from m6)
plot(Prev ~ as.numeric(Study), data = d, pch = 19, las = 1, xaxt = "n",
     col = ifelse(d$Habitat == "Ag", "brown", "darkgreen"),
     ylim = c(0, 0.8), xlab = "Study", ylab = "Observed Prevalence")
axis(1, at = 1:11, labels = 1:11)

legend("topright", pch = 19, col = c("brown", "darkgreen"), 
       legend = unique(d$Habitat), bty = "n", cex = 1.2)

for (x in 1:11){
  points(x, logistic(fixef(m6)[1] + fixef(m6)[2]*0 + ranef(m6)$Study[x,1]), 
         col = "brown", pch = 2)
}
for (x in 1:11){
  points(x, logistic(fixef(m6)[1] + fixef(m6)[2]*1 + ranef(m6)$Study[x,1]), 
         col = "darkgreen", pch = 2)
}

# So the hollow triangles are showing mean predicted prevalence by Habitat
# for each study using estimates derived from m6. Notice that the green
# Forest triangle is always lower than the Ag triangle. This model assumes
# there is a fixed effect of Habitat that holds across all Studies. So while
# the random effects intercepts are pushing overall prevalence higher or lower
# as you look across studies, the fixed effect of Habitat is always going to
# indicate that Forest has lower prevalence than Ag (for this model!). For
# some Studies, these predicted prevalence values are extremely close to the
# observed values.

#==============================================================================


# More advanced plotting


# Loop to create a histogram of predicted prevalence values by Habitat. 
# This is a little complicated, but stay with me!

# First, create a dummy variable representing Forest habitat (equals 1 if
# Habitat is Forest)
d$is_forest <- ifelse(d$Habitat == "Forest", 1, 0)

# Clear the graphics device and set up a counter
dev.off()
counter <- 0

# Our predictions are going to be generated from parameter values drawn 
# from the rethinking package function "sample.qa.posterior." The default 
# is that this function gives us 10,000 sets of parameter samples from a
# given model. Importantly, these parameter sets account for covariance 
# in the parameter estimates. Using these parameter sets will allow us to
# accurately represent the uncertainty in our parameter estimates!
post <- sample.qa.posterior(m6)

# The following loop is used to plot a histogram of prevalence predictions 
# for one habitat, then it adds on the histogram for the other habitat
for (x in levels(d$Habitat)) {
  
  # Calculate the total number of observations for Habitat x in the 
  # original data (or you could choose something arbitrary like 100 or 
  # 1000 observations for each habitat)
  sample_size <- sum(d$Num_Tested[d$Habitat == x])
  
  # We want to generate a number of random binomial samples equal to the 
  # number of parameter sets we generated using "sample.qa.posterior"
  binom_samples <- nrow(post)
  
  # This is just a dummy variable telling us whether we're considering
  # Agricultural or Forest habitat
  forest <- mean(d$is_forest[d$Habitat == x])
  
  # Print info about group x
  print(c(x, sample_size))
  
  # Assign a new variable for predictions regarding habitat x. As you can 
  # see, the predictions are just random samples from the binomial 
  # distribution using probability of success values that are calculated 
  # using our 10,000 parameter sets from "sample.qa.posterior" and size 
  # values = to the "sample_size" variable we defined before. The predictions
  # themselves are raw counts of successes. Divide this vector of 
  # success counts by the size of the sample to get prevalence values.
  assign(paste0("preds.", gsub(" ","_",x)), 
         (rbinom(binom_samples, 
                 prob = logistic(post[ , 1] + forest*post[ , 2]),
                 size = sample_size)
          /sample_size))
  
  # Assign these predictions to "current_sample" each pass through the for loop
  current_sample <- get(paste0("preds.", gsub(" ","_",x)))
  
  # Below are some histogram settings we want to use
  bin.vec <- seq(0, 1, 0.005)
  xlimits <- c(0.1, 0.7)
  ylimits <- c(0, 1000)
  
  # Add one to the counter each time through the for loop
  counter <- counter + 1
  
  # Plot the histograms. The first time through the for loop 
  # (when counter == 1), plot the histogram and legend. Each time after 
  # that just add a histogram of the current predictions to the existing plot.
  if (counter == 1) {
    
    hist(current_sample, breaks = bin.vec, 
         col = adjustcolor("brown", alpha.f = 0.8), 
         xlim = xlimits, ylim = ylimits, main = "",
         xlab = "Prevalence in Simulated Sample (Predictions)")
  
    legend("topright", pch = 19, col = c("brown", "darkgreen"), 
           legend = unique(d$Habitat), bty = "n", cex = 1.3)
  }
  
  else {
    
    hist(current_sample, breaks = bin.vec, 
         col = adjustcolor("darkgreen", alpha.f = 0.8), add = T)
  }
}
# Rerun the loop multiple times (from the dev.off() line) to see the 
# stochasticity in our prevalence predictions! This stochasticity emerges
# because of uncertainty in our parameter estimates and because of randomness
# in the binomial process (the same probability of success value can result
# in either failure or sucess in the binomial process unless the
# probability of success equals 0 or 1). You can also look at the objects 
# "preds.Ag" or "preds.Forest" which hold all the prevalence predictions
# generated by the loop.


# Instead of histograms of just the predictions, we can also plot the 
# predictions along with the observed data. First, we just make a plot 
# similar to thosewe made before with the mean estimated probability of 
# success for each habitat type shown as solid horizontal lines (again 
# parameter estimates are coming from m6).
plot(Prev ~ as.numeric(Study), data = d, pch = 19, las = 1, xaxt = "n",
     col = ifelse(d$Habitat == "Ag", "brown", "darkgreen"),
     xlim = c(0,13), ylim = c(0,0.8),
     xlab = "Study Observations + Predictions", ylab = "Prevalence")
axis(1, at = 1:13, labels = c(1:11, "Ag", "Forest"))

legend("topleft", pch = 19, col = c("brown", "darkgreen"), 
       legend = unique(d$Habitat), bty = "n", cex = 1.2)

abline(h = logistic(fixef(m6)[1]), col = "brown", lwd = 3)
abline(h = logistic(fixef(m6)[1] + fixef(m6)[2]), col = "darkgreen", lwd = 3)

# Plot all predictions for Agricultural habitats along with mean and 
# 95% highest posterior density intervals (HDPI) of these predictions. 
# The beauty of having these prediction vectors is that they can be easily 
# summarized. What do the predictions suggest about the mean prevalence for 
# a habitat? Just calculate the mean of the predictions! What about 
# confidence intervals or similar metrics? Just calculate those intervals 
# on the predictions!
points(preds.Ag ~ rep(12, length(preds.Ag)), col = "brown", pch = 19)
points(12, mean(preds.Ag), pch = 19)
points(12, HPDI(preds.Ag, prob = 0.95)[1], pch = 20)
points(12, HPDI(preds.Ag, prob = 0.95)[2], pch = 20)

# Plot all predictions for Forest habitats along with mean and 95% HDPI of 
# these predictions
points(preds.Forest ~ rep(13, length(preds.Forest)), 
       col = "darkgreen", pch = 19)
points(13, mean(preds.Forest), pch = 19)
points(13, HPDI(preds.Forest, prob = 0.95)[1], pch = 20)
points(13, HPDI(preds.Forest, prob = 0.95)[2], pch = 20)
