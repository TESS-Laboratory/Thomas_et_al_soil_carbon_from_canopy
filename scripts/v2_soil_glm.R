## A script to use GLMs to model and analyse the mechanistic effects of canopy structure on Soil C

#### set up environment ####
library(performance)

#TODO pick method of narrowing values
df_numeric <- converted_data %>%
  select(where(is.numeric))

## plot without categorical data 
rda.out <- vegan::rda(df_numeric, scale = TRUE)
# add scores()
rda_scores <- scores(rda.out)
# add biplot()
biplot(rda.out, type = "text")
## group by depth ### this doesnt work yet
ordihull(rda.out,
         group = Fmetrics_num$Profundidade_cm_,
         col = 1:11,
         lty = 1:11,
         lwd = c(3,6), 
         label = TRUE)



# option 1 frequency density 
df<- read.csv("Data/variable_presence_count.csv")


fixed_vars <- df$value[4:25]

#create a formula for chosen variables
formula_str <- paste("wmean_percC_5 ~", paste(fixed_vars, collapse = " + "))
glm_formula <- as.formula(formula_str)

##compare performance of different family and links
glm.gaus.log <- glm(glm_formula, data = converted_data, family= gaussian('log'))
glm.gaus.id <- glm(glm_formula, data = converted_data, family = gaussian('identity'))
glm.gam.log <- glm(glm_formula, data = converted_data, family = Gamma(link="log"))
glm.gam.id<- glm(glm_formula, data = converted_data, family = Gamma(link="identity"))

## review
compare_performance(glm.gam.log, glm.gaus.id, glm.gaus.log, glm.gam.id, verbose = FALSE, rank = TRUE)
summary(glm.gam.log)

check_model(glm.gam.log)

plot(resid(glm.gam.log, type='response'))
lines(resid(glm.gam.log, type='response'), col='red')

##selected formula
glm_model <- glm(glm_formula, data = converted_data, family = Gamma(link="log"))

summary(glm_model)


# Perform automated backward selection
best_model <- step(glm_model, direction = "backward")

#View the summary of the best model
summary(best_model)

## plot effect sizes 
coef_summary <- summary(best_model)$coefficients
conf_int <- confint(best_model)

# Create a data frame for plotting
effects_df <- data.frame(
  Predictor = rownames(coef_summary),
  Estimate = coef_summary[, "Estimate"],
  CI_low = conf_int[, 1],
  CI_high = conf_int[, 2]
)

# remove the intercept to focus only on predictors
effects_df <- effects_df[effects_df$Predictor != "(Intercept)", ]

# Plot
ggplot(effects_df, aes(x = reorder(Predictor, Estimate), y = Estimate)) +
  geom_point() +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = 0.2) +
  coord_flip() +
  theme_minimal() +
  labs(
    title = "Effect Sizes (Log Scale) from Gamma GLM",
    x = "Predictor",
    y = "Log Coefficient Estimate ± 95% CI"
  )

