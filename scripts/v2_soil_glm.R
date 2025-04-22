## A script to use GLMs to model and analyse the mechanistic effects of canopy structure on Soil C

## pick method of narrowing values
# frequency density 
df<- read.csv("Data/variable_presence_count.csv")
fixed_vars <- df$value[4:25]

formula_str <- paste("wmean_percC_5 ~", paste(fixed_vars, collapse = " + "))
glm_formula <- as.formula(formula_str)


glm_model <- glm(glm_formula, data = converted_data, family = gaussian())

summary(glm_model)


# Perform backward selection
# direction = "backward" ensures the model starts with all predictors
best_model <- step(glm_model, direction = "backward")

#View the summary of the best model
summary(best_model)

coef_summary <- summary(best_model)$coefficients
conf_int <- confint(best_model)

# Create a data frame for plotting
effects_df <- data.frame(
  Predictor = rownames(coef_summary),
  Estimate = coef_summary[, "Estimate"],
  CI_low = conf_int[, 1],
  CI_high = conf_int[, 2]
)

# Optional: remove the intercept if you want to focus only on predictors
effects_df <- effects_df[effects_df$Predictor != "(Intercept)", ]

# Plot
ggplot(effects_df, aes(x = reorder(Predictor, Estimate), y = Estimate)) +
  geom_point() +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = 0.2) +
  coord_flip() +
  theme_minimal() +
  labs(
    title = "Effect Sizes from Best GLM",
    x = "Predictor",
    y = "Coefficient Estimate (with 95% CI)"
  )
