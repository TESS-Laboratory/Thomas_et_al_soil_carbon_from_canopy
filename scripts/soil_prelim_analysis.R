### script for prelim soil analysis
##### set up environment ####
library(readr)
install.packages(c("dplyr", "ggplot2", "patchwork", "broom"))
library(dplyr)
library(ggplot2)
library(patchwork)
library(broom)

##### read in data #####
fp<- "C:/Users/jpt215/OneDrive - University of Exeter/PhD_Data/Soil_manuscript_data"
soil_samples <- read_csv(file.path(fp, "Jess_Plinio_soil_samples.csv"), 
                                     col_types = cols(massa = col_number(), 
                                                      `15N` = col_number(),
                                                      `%N` = col_number(), 
                                                      `13C` = col_number(),
                                                      `%C` = col_number(), 
                                                      `C/N` = col_number(), 
                                                      ...13 = col_skip(), 
                                                      ...14 = col_skip(),
                                                      ...15 = col_skip(), 
                                                      ...16 = col_skip()))
#View(soil_samples)


LAI_table <- read_csv(file.path(fp, "LAI_w_meta_data.csv"), 
                      col_types = cols(ID = col_character()))


#####combine tables ####


clean_column <- function(column) {
  return(gsub("[A-Za-z]", "", column))
}

# Clean the 'id_A' column in Table A
soil_samples$id_clean <- clean_column(soil_samples$Ponto)

##clean Ids in LAI table 
LAI_table <- LAI_table %>%
  filter(!grepl("[a-zA-Z]", ID))

colnames(LAI_table) <- paste0(colnames(LAI_table), "_4")
colnames(soil_samples) <- paste0(colnames(soil_samples), "_5")
# Bind the tables based on the ID
merged_table <- soil_samples%>%
  left_join(LAI_table, by = c("id_clean_5" = "ID_4"))
cleaned_table <- merged_table %>%
  filter(!is.na(effective.LAI_4))

merged_table<- merged_table|> select(-c("...1_4"))

write.csv(merged_table, file =file.path(fp, "soil_meta_table.csv"))


##### box plots #####
# Filter the data for soil depth
filter_values <- c("20-30")
filtered_data <- cleaned_table %>%
  filter(`Profundidade (cm)` %in% filter_values)

# Create box plots of SOC separated by LAI (Low, Middle, High)
box_plot <- ggplot(filtered_data, aes(x = Degradation, y = `%C`)) +
  geom_boxplot() +
  labs(title = "Box Plot of %carbon by degredation history (20-30 cm)",
       x = "LAI",
       y = "%SOC") +
  theme_minimal()

# Print the box plot
print(box_plot)

#create boxplots of SOC separated by degredation history
# Filter the data for soil depth
filter_values <- c("20-30","0-5","5-10","10-20")
filtered_data <- cleaned_table %>%
  filter(`Profundidade (cm)` %in% filter_values)


# Create categories for low, middle, and high values of LAI
filtered_data <- filtered_data %>%
  mutate(C_Category = cut(effective.LAI, 
                          breaks = quantile(effective.LAI, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE), 
                          labels = c("Low", "Middle", "High"), 
                          include.lowest = TRUE))

# Create box plots of SOC separated by LAI (Low, Middle, High)
box_plot <- ggplot(filtered_data, aes(x = C_Category, y = `%C`)) +
  geom_boxplot() +
  labs(title = "Box Plot of %carbon by Low, Middle, and High Values of LAI (0-30cm)",
       x = "LAI",
       y = "%SOC") +
  theme_minimal()

# Print the box plot
print(box_plot)



#### line plots ####
# Filter the data for soil depth
filter_values <- c("0-5","5-10","10-20","20-30")
filtered_data <- cleaned_table %>%
  filter(`Profundidade (cm)` %in% filter_values)

# Plot the data from LAI and SOC
y_limits <- c(0, 10)  
plot <- ggplot(filtered_data, aes(x = effective.LAI, y = `%C`)) +
  geom_point() +                        # Add points for the data
  geom_smooth(method = "lm", se = FALSE, color = "blue") +  # Add line of best fit
  labs(title = "Scatter Plot  of %SOC and LAI with Linear model 20-30cm",
       x = "LAI",
       y = "%SOC") +
  ylim(y_limits) +
  theme_minimal()

# Print the plot
print(plot)

# Fit a linear model (line of best fit)
fit <- lm(`%C` ~ effective.LAI, data = filtered_data)

# Print the summary of the linear model
summary(fit)

## add weights for means
weight_table <- filtered_data %>%
  mutate(weight = case_when(
    `Profundidade (cm)`== "0-5" ~ 1,
    `Profundidade (cm)` =="5-10" ~ 1,
    `Profundidade (cm)` == "10-20" ~ 2,
    `Profundidade (cm)` == "20-30" ~2
  ))

# Summarize the data by unique IDs with a weighted mean of SOC
summary_data <- weight_table %>%
  group_by(id_clean) %>%
  summarize(`%C` = sum(`%C` * weight) / sum(weight),
            effective.LAI = mean(effective.LAI),
            Degradation = modal(Degradation),
            `Profundidade (cm)` = "0-30")

## plot comined graph
plot <- ggplot(summary_data, aes(x = effective.LAI, y = weighted_mean_C)) +
  geom_point() +                        # Add points for the data
  geom_smooth(method = "lm", se = FALSE, color = "blue") +  # Add line of best fit
  labs(title = "Scatter Plot  of %SOC and LAI with Linear model, weighted mean of 0-30cm",
       x = "LAI",
       y = "%SOC") +
  ylim(y_limits) +
  theme_minimal()
# Print the plot
print(plot)

# Fit a linear model (line of best fit)
fit <- lm(weighted_mean_C ~ effective.LAI, data = summary_data)

# Print the summary of the linear model
summary(fit)

# Create box plots of SOC separated by LAI (Low, Middle, High)
box_plot <- ggplot(summary_data, aes(x = Degradation, y = weighted_mean_C)) +
  geom_boxplot() +
  labs(title = "Box Plot of %carbon by degredation history (0-30cm)",
       x = "LAI",
       y = "%SOC") +
  theme_minimal()

# Print the box plot
print(box_plot)


########### functionise ###########
# Function to create box plots separated by a given category (LAI or Degradation)
create_box_plot <- function(data, profundidade, category, value_label) {
  filtered_data <- data %>%
    filter(`Profundidade (cm)` == profundidade)
  
  box_plot <- ggplot(filtered_data, aes(x = category, y = "`%C`")) +
    geom_boxplot() +
    labs(title = paste("Box Plot of %SOC by", value_label, "(Profundidade:", profundidade, "cm)"),
         x = value_label,
         y = "%SOC") +
    theme_minimal()
  
  return(box_plot)
}

# Function to perform t-tests between groups in a given category
perform_t_tests <- function(data, profundidade, category) {
  filtered_data <- data %>%
    filter(`Profundidade (cm)` == profundidade)
  
  category_levels <- unique(filtered_data[[category]])
  t_test_results <- list()
  
  for (i in 1:(length(category_levels) - 1)) {
    for (j in (i + 1):length(category_levels)) {
      group1 <- filtered_data %>%
        filter(get(category) == category_levels[i]) %>%
        pull(`%C`)
      group2 <- filtered_data %>%
        filter(get(category) == category_levels[j]) %>%
        pull(`%C`)
      
      t_test <- t.test(group1, group2)
      t_test_results[[paste(category_levels[i], "vs", category_levels[j])]] <- tidy(t_test)
    }
  }
  
  return(t_test_results)
}

# Function to create line plots and fit linear models
create_line_plot <- function(data, profundidade) {
  filtered_data <- data %>%
    filter(`Profundidade (cm)` == profundidade)
  
  plot <- ggplot(filtered_data, aes(x = effective.LAI, y = `%C`)) +
    geom_point() + 
    geom_smooth(method = "lm", se = FALSE, color = "blue") +
    labs(title = paste("Scatter Plot of %SOC and LAI with Linear model (Profundidade:", profundidade, "cm)"),
         x = "LAI",
         y = "%SOC") +
    theme_minimal()
  
  fit <- lm(`%C` ~ effective.LAI, data = filtered_data)
  
  return(list(plot = plot, model_summary = summary(fit)))
}

#prep data 
filter_values <- c("0-5","5-10","10-20","20-30")
filtered_data <- cleaned_table %>%
  filter(`Profundidade (cm)` %in% filter_values)
filtered_data <- filtered_data %>%
  mutate(C_Category = cut(effective.LAI, 
                          breaks = quantile(effective.LAI, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE), 
                          labels = c("Low", "Middle", "High"), 
                          include.lowest = TRUE))
## add weights for means
weight_table <- filtered_data %>%
  mutate(weight = case_when(
    `Profundidade (cm)`== "0-5" ~ 1,
    `Profundidade (cm)` =="5-10" ~ 1,
    `Profundidade (cm)` == "10-20" ~ 2,
    `Profundidade (cm)` == "20-30" ~2
  ))

# Summarize the data by unique IDs with a weighted mean of SOC
summary_data <- weight_table %>%
  group_by(id_clean) %>%
  summarize(`%C` = sum(`%C` * weight) / sum(weight),
            effective.LAI = mean(effective.LAI),
            Degradation = modal(Degradation),
            `Profundidade (cm)` = "0-30")


# Apply functions to each value of `Profundidade (cm)` and `LAI`
profundidades <- c("0-5", "5-10", "10-20", "20-30")
box_plots_LAI <- lapply(profundidades, create_box_plot, data = filtered_data, category = "C_Category", value_label = "LAI")
t_tests_LAI <- lapply(profundidades, perform_t_tests, data = filtered_data, category = "C_Category")
line_plots_LAI <- lapply(profundidades, create_line_plot, data = filtered_data)

# Create weighted mean C box plot and perform t-tests by LAI
weighted_box_plot_LAI <- ggplot(summary_data, aes(x = C_Category, y = weighted_mean_C)) +
  geom_boxplot() +
  labs(title = "Box Plot of Weighted %SOC by LAI (0-30cm)",
       x = "LAI",
       y = "%SOC") +
  theme_minimal()


# Apply functions to each value of `Profundidade (cm)` and `Degradation`
box_plots_Degradation <- lapply(profundidades, create_box_plot, data = cleaned_table, category = "Degradation", value_label = "Degradation")
t_tests_Degradation <- lapply(profundidades, perform_t_tests, data = cleaned_table, category = "Degradation")

# Create weighted mean C box plot and perform t-tests by Degradation
weighted_box_plot_Degradation <- ggplot(summary_data, aes(x = Degradation, y = weighted_mean_C)) +
  geom_boxplot() +
  labs(title = "Box Plot of Weighted %SOC by Degradation (0-30cm)",
       x = "Degradation",
       y = "%SOC") +
  theme_minimal()

t_tests_weighted_Degradation <- perform_t_tests(summary_data, profundidade = "0-30", category = "Degradation")

# Display all the line plots on one screen
line_plots_combined <- wrap_plots(lapply(line_plots_LAI, function(x) x$plot))

# Print the box plots, line plots, and t-test results
print(wrap_plots(box_plots_LAI))
print(weighted_box_plot_LAI)
print(t_tests_LAI)
print(t_tests_weighted_LAI)

print(wrap_plots(box_plots_Degradation))
print(weighted_box_plot_Degradation)
print(t_tests_Degradation)
print(t_tests_weighted_Degradation)

print(line_plots_combined)

# Print linear model summaries
lapply(line_plots_LAI, function(x) print(x$model_summary))
