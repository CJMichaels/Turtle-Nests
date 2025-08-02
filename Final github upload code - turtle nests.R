###Clear environment

rm(list = ls())

###Notes: PCA is OK with repeated measure as it is a data reduction tool
### not an inferential tool : https://www.sciencedirect.com/science/article/abs/pii/S0167811602000654

###Packages###
library(readxl)
library(snakecase)
library(dplyr)
library(dplyr)
library(FactoMineR)
library(factoextra)
library(lme4)
library(lmerTest)
library(DHARMa)
library(ggplot2)
library(rstan)
library(bayesplot)
library(loo)
library(robust)
library(dplyr)
library(corrplot)
library(ggcorrplot)
library(patchwork)


###Functions###

# Make colnames snake case, make NRs = NA and remove any blank rows ------------------------------
# x is the dataframe to which you want to apply the function

snake_case_NAs <- function(x){
  dataframe <- x
  names(dataframe)<- to_snake_case(names(dataframe)) #Make colnames snake case
  
  dataframe[dataframe == "NR"] <- NA #make NAs = NR
  
  dataframe <- dataframe %>% filter_all(any_vars(!is.na(.))) #remove any blank rows
  
  return(dataframe)
}

#function to centre cols
center_columns <- function(data, columns_to_center) {
  # Loop through each specified column
  for (col in columns_to_center) {
    # Center the column by subtracting the mean
    data[[col]] <- data[[col]] - mean(data[[col]], na.rm = TRUE)
  }
  return(data)
}

###Import Data###

turtle_data <- read_excel("turtle_data.xlsx")

sand_data <- read_excel("Sand particle_Final.xlsx")

colnames(sand_data) <- sand_data[1,]
colnames(sand_data) <- gsub (" ", "_", colnames(sand_data))

sand_data <- sand_data[-1,]

sand_data <- sand_data %>% select(location, contains("mm"), contains("μ"), contains("Year"), Type) %>%
  rename(location_sand = location, year_sand = Year, Type_sand = Type)

turtle_data_all <- cbind(turtle_data, sand_data) %>% select(-Type_sand, -location_sand, -year_sand)

turtle_data <- turtle_data_all

###Cleaning###

turtle_data <- snake_case_NAs(turtle_data) #fix data headings and NA entries

turtle_data <- turtle_data %>% 
  mutate(location = gsub("Baruli", "Barauli", location)) #correct spelling error

#select columns to centre         
cols_to_centre <- c("slope", "moisture", 
                    "distance_m", "total_nest_depth", "total_water_depth",
                    "width", "ph")
sand_selected <- colnames(turtle_data)[grep("mm|μ", colnames(turtle_data))]

turtle_data <- turtle_data %>%
  mutate(across(.cols = all_of(sand_selected), ~ as.numeric(.) ))

turtle_data <- turtle_data %>%
  mutate(across(.cols = all_of(sand_selected), ~ (./500 *100)))

turtle_data_save <- turtle_data

#turtle_data <- center_columns(turtle_data, c(cols_to_centre, sand_selected)) #centre selected columns - this helps with model fit and interpretation

turtle_data %>% group_by(location, type )%>% #summarise outcome by location
  summarise(number = n())

ggplot(turtle_data, aes(x = location, fill = type))+
  geom_histogram(stat = "count")+
  theme_bw()



###PCA###

# Subset active variables for the PCA 
df_pca <- turtle_data %>% select(-c(year, type, location, nest_no, wet, dry, total_nest_depth)) #negative selection to remove vars not for PCA

#check correlations to justify PCA

# Compute the standard correlation matrix using Pearson's method
correlation_matrix_all_start <- cor(df_pca, use = "pairwise.complete.obs")

# Print the correlation matrix
print(correlation_matrix_all_start)

# Define custom variable labels

variable_labels <- c("Slope", "Moisture Content", "Distance", "Water Depth", "Nest Site Width", "pH", 
                     "% 6.3mm",  "% 2.7mm",  "% 2.36mm", "% 1.18mm", "% 600μ",   "% 300μ",   "% 150μ",   "% 90μ",    "% 75μ",    "% 45μ"   )


colnames(correlation_matrix_all_start) <- variable_labels
rownames(correlation_matrix_all_start) <- variable_labels

# Visualize the correlation matrix with specified customizations
corrplot(correlation_matrix_all_start, method = "color", 
         type = "upper",   # Show only upper triangle
         addCoef.col = "black",  # Add coefficient labels in black
         tl.col = "black", tl.srt = 45, # Text label color and rotation
         diag = FALSE)+ # Use custom labels
  guides(fill = guide_legend(title = "Pearson's R")) 


# PCA (FactoMineR automatically standardizes the data so st.dev=1 and mean=0)
res.pca <- PCA(df_pca, graph = FALSE, ncp = Inf) #run PCA
print(res.pca) #print results
eig.val <- get_eigenvalue(res.pca) #extract eigenvalues
fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 50)) # scree plot to help inform dim selection

#Extract results from PCA output (list of matrices combining results for included vars)
var <- get_pca_var(res.pca)

#Extract relevant matrices See http://factominer.free.fr/question/FAQ.html
coords <- var$coord    #principal component scores
loadings <- sweep(res.pca$var$coord,2,sqrt(res.pca$eig[1:ncol(res.pca$var$coord),1]),FUN="/") #loadings
contributions <-  var$contrib #contributions

#select dimensions that have ev >=1 and cum var >70%
dimensions_selected <- as.data.frame(eig.val) %>%
  filter(cumulative.variance.percent <70) %>% #threshold set a bit higher - check manually and adjust if required
  filter(eigenvalue >=1)%>%
  mutate(dimensions = row.names.data.frame(.))%>%
  pull(dimensions) 

#decide groupings by allocating variable to dimension where it has highest loadings
loadings_selected <- as.data.frame(loadings) %>% 
  select(all_of(dimensions_selected))

variables <- as.vector(row.names(loadings_selected))

groupings_table_loadings <- data.frame(var = variables, dimension = NA)

for(i in 1: length(variables)){
  variable <- variables[i]
variable_of_interest <- loadings_selected[variable, ]
highest_loading_dimension <- as.vector(which.max(abs(variable_of_interest)))
groupings_table_loadings[groupings_table_loadings$var == variable, 2] <- highest_loading_dimension
}

groupings_table_loadings <- groupings_table_loadings %>%
  arrange(dimension, var)


#decide groupings by allocating variable to dimension where it has highest contributions
contributions_selected <- as.data.frame(contributions) %>% 
  select(all_of(dimensions_selected))

variables <- as.vector(row.names(contributions_selected))

groupings_table_contributions <- data.frame(var = variables, dimension = NA)

for(i in 1: length(variables)){
  variable <- variables[i]
  variable_of_interest <- contributions_selected[variable, ]
  highest_contributions_dimension <- as.vector(which.max(abs(variable_of_interest)))
  groupings_table_contributions[groupings_table_contributions$var == variable, 2] <- highest_contributions_dimension
}

groupings_table_contributions <- groupings_table_contributions %>%
  arrange(dimension, var)

# COmpare
groupings_table_contributions == groupings_table_loadings

##plots
# Extract PCA scores
#pca_scores <- cbind(turtle_data$type, as.data.frame(res.pca$ind$coord))

pca_scores <- as.data.frame(res.pca$ind$coord)[, 1:5]
loadings_for_graph <- as.data.frame(loadings)[, 1.:5]

# Load necessary libraries

library(gridExtra)
# Function to select relevant variables based on a threshold
# Function to clean variable names
clean_variable_name <- function(name) {
  # Convert to title case and remove underscores
  name_clean <- gsub("_", " ", name)
  name_clean <- tools::toTitleCase(name_clean)
  
  # Replace spaces between numbers with dashes
  name_clean <- gsub("([0-9]+) ([0-9]+)", "\\1-\\2", name_clean)
  
  # Specific replacements
  name_clean <- gsub("\\bm\\b", "", name_clean, ignore.case = TRUE) # Remove lone 'm'
  name_clean <- gsub("\\bPh\\b", "pH", name_clean, ignore.case = TRUE) # Convert 'ph' to 'pH'
  name_clean <- gsub("\\bMm\\b", "mm", name_clean) # Keep 'mm' as 'mm'
  
  return(trimws(name_clean))
}

# Function to select relevant variables based on a threshold
select_relevant_vars <- function(loadings, threshold = 0.3) {
  # Identify columns (dimensions) in the loadings table
  dimensions <- colnames(loadings)
  relevant_vars_list <- list()
  
  # Iterate through each dimension
  for (dim in dimensions) {
    relevant_vars <- rownames(loadings)[abs(loadings[[dim]]) >= threshold]
    relevant_vars_list[[dim]] <- relevant_vars
  }
  
  return(relevant_vars_list)
}

# Customizable threshold value for loadings
threshold_value <- 0.3  # Adjust based on your analysis needs
# Identify relevant variables with a threshold
relevant_vars_list <- select_relevant_vars(loadings_for_graph, threshold = threshold_value)

# Function to plot correlation matrices using ggcorrplot with cleaned variable names
plot_correlation_matrix <- function(data, variables, dim_name) {
  if (length(variables) > 1) {
    # Clean variable names
    cleaned_names <- sapply(variables, clean_variable_name)
    
    # Update column names in the data
    cleaned_data <- data
    colnames(cleaned_data)[match(variables, colnames(data))] <- cleaned_names    
    # Calculate the correlation matrix
    correlation_matrix <- cor(cleaned_data[, cleaned_names, drop = FALSE])
    
    # Plot the correlation matrix
    ggcorrplot(correlation_matrix, hc.order = TRUE, type = "upper",
               lab = TRUE, lab_col = "white", lab_size = 3, 
               colors = c("#1E90FF", "white", "#FF4500"),
               title = paste("Correlation Matrix for", dim_name)) +
      theme(plot.title = element_text(hjust = 0.5, size = 10))
  } else {
    ggplot() + 
      annotate("text", x = 0.5, y = 0.5, label = paste("Not enough variables for", dim_name), size = 5) + 
      theme_void()
  }
}

# Create a list to hold ggplot objects
plot_list <- list()

# Assuming you have the original dataset used for PCA in data frame `df_pca`
# Generate correlation matrix plots for each dimension with cleaned variable names
for (dim in names(relevant_vars_list)) {
  p <- plot_correlation_matrix(df_pca, relevant_vars_list[[dim]], dim)
  plot_list <- c(plot_list, list(p))
}

# Filter out NULL elements that may have been added for dimensions with insufficient variables
plot_list <- Filter(Negate(is.null), plot_list)

# Arrange plots in a 3x2 grid and save as a grob (graphical object)
combined_plot <- do.call(gridExtra::grid.arrange, c(plot_list[1:min(6, length(plot_list))], ncol = 3, nrow = 2))


####Logistic Regression####

#extract reduced data

ind <- get_pca_ind(res.pca) #pca data
reduced_data <- ind$coord
reduced_data <- cbind(turtle_data, reduced_data)

keep_cols <- colnames(turtle_data)[!colnames(turtle_data) %in% colnames(df_pca)] #identify some columns to retain based on PCA df
keep_cols <- c(keep_cols, dimensions_selected) #update to include columns for retained dimensions, too
keep_data <- reduced_data %>% select(all_of(keep_cols), -wet, -dry) #create retained dataset with only retained columns
colnames(keep_data) <- to_snake_case(colnames(keep_data)) #fix headings


###Logistic regression
turt_centred <- center_columns(keep_data, c("dim_1", "dim_2", "dim_3", "dim_4", "dim_5")) #re-centre dim cols

analysis_data_agg <- turt_centred %>% #create aggregate dataframe; given structure of data collection
  group_by(location) %>%              # this helps to address combination of random effect of location
  summarise(                          # and complete separation of outcome by location
    trials = n(),                     # boils data down to trials and successes with means of covariates
    success = sum(type == "Preferred"),
    dim_1_mean = mean(dim_1, na.rm = TRUE),
    #dim_2_mean = mean(dim_2, na.rm = TRUE), #note dim_1 only retained for final model (others n.s.)
    #dim_3_mean = mean(dim_3, na.rm = TRUE),
    dim_1_sd = sd(dim_1, na.rm = TRUE)#,
    #dim_2_sd = sd(dim_2, na.rm = TRUE),
    #dim_3_sd = sd(dim_3, na.rm = TRUE),
    #dim_4_mean = mean(dim_4, na.rm = TRUE),
    #dim_5_mean = mean(dim_5, na.rm = TRUE),
    #dim_6_mean = mean(dim_6, na.rm = TRUE),
    #dim_4_sd = sd(dim_4, na.rm = TRUE),
    #dim_5_sd = sd(dim_5, na.rm = TRUE)#,
    #dim_6_sd = sd(dim_6, na.rm = TRUE)
   ) %>%
  mutate(success = ifelse(success == 0, 0, 1))


# Prepare data for Stan
# Prepare the data list for Stan

# Prepare the data list for Stan
stan_data <- list(
  N = nrow(analysis_data_agg), # Number of data points (sites)
  y = analysis_data_agg$success, # Binary outcome indicating overall success (1) or failure (0)
  K = length(grep("dim_.*_mean", colnames(analysis_data_agg))), # Number of covariates
  x_mean = as.matrix(analysis_data_agg[, grep("dim_.*_mean", colnames(analysis_data_agg))]), # Matrix of covariate means
  x_sd = as.matrix(analysis_data_agg[, grep("dim_.*_sd", colnames(analysis_data_agg))]), # Matrix of covariate standard deviations
  n_obs = matrix(
    analysis_data_agg$trials, 
    nrow = nrow(analysis_data_agg), 
    ncol = length(grep("dim_.*_mean", colnames(analysis_data_agg)))
  ) # Matrix of sample sizes
)


stan_model_code <- "
data {
  int<lower=0> N;                     // Number of observations
  int<lower=0, upper=1> y[N];         // Binary outcome (1 for nesting, 0 for no nesting)
  int<lower=1> K;                     // Number of covariates
  matrix[N, K] x_mean;                // Mean values of covariates for each observation
  matrix<lower=0>[N, K] x_sd;         // Standard deviations of covariates
  matrix<lower=1>[N, K] n_obs;        // Sample sizes for each covariate
}
parameters {
  real alpha;                         // Intercept
  vector[K] beta;                     // Coefficients for the covariates
}
model {
  vector[N] eta;                      // Linear predictor

  // Compute linear predictor
  eta = alpha + x_mean * beta;

  // Adjust for measurement uncertainty with inverse variance weighting
  for (n in 1:N) {
    for (k in 1:K) {
      real variance = square(x_sd[n, k]) / n_obs[n, k]; // Variance of the mean
      eta[n] += 0.5 * beta[k] * beta[k] * variance;     // Weighted adjustment
    }
  }

  // Likelihood
  y ~ bernoulli_logit(eta);

  // Priors
  alpha ~ normal(0, 10);
  beta ~ normal(0, 1.5);              // Prior for all beta coefficients
}
generated quantities {
  vector[N] log_lik;                  // Store log likelihood of each observation for model evaluation

  for (n in 1:N) {
    real adjustment = 0;

    // Recalculate adjustment term for generated quantities with weighting
    for (k in 1:K) {
      real variance = square(x_sd[n, k]) / n_obs[n, k]; // Variance of the mean
      adjustment += 0.5 * beta[k] * beta[k] * variance; // Weighted adjustment
    }

    // Recompute linear predictor
    real eta_n = alpha + dot_product(beta, x_mean[n]) + adjustment;

    // Log likelihood
    log_lik[n] = bernoulli_logit_lpmf(y[n] | eta_n);
  }
}



"


# 
#Fit the Stan model
### Note earlier STan model showed no sig effect od dims 2 to 6, so these removed and refitted. LOOIC lower for reduced model

 fit_fe <- stan(
   model_code = stan_model_code,  # Stan model code
   data = stan_data,
  iter = 80000, warmup = 20000,
   chains = 4,
   control = list(adapt_delta = 0.99, max_treedepth = 15)
 )

# Save the fit object to a file
# saveRDS(fit_fe, file = "final_turtle_nest_model.rds")

pars = c("beta", "alpha", "lp__")

fit <- fit_fe # update name to match code 

# Print the summary of the model
fit_df <- as.data.frame(summary(fit, pars = pars)$summary)

# Extract log-likelihood values from the fitted model
log_lik <- extract_log_lik(fit, parameter_name = "log_lik", merge_chains = FALSE)

# Calculate LOO
loo_result <- loo(log_lik)

# Calculate WAIC
waic_result <- waic(log_lik)

###DIagnositics###

posterior_samples <- rstan::extract(fit)
alpha_draws <- posterior_samples$alpha
beta_draws <- posterior_samples$beta

alpha_mean <- mean(alpha_draws)
beta_mean <- mean(beta_draws)
linear_predictor <- alpha_mean + beta_mean * analysis_data_agg$dim_1_mean
predicted_prob <- 1 / (1 + exp(-linear_predictor))

plot(analysis_data_agg$dim_1_mean, predicted_prob, type = "b",
     xlab = "PC1",
     ylab = "Predicted Probability of Nesting",
     main = "Predicted Probability vs PC1")
points(analysis_data_agg$dim_1_mean, analysis_data_agg$success, col = "red", pch = 19)
legend("topright", legend = c("Predicted Probability", "Observed Outcome"), 
       col = c("black", "red"), pch = c(1, 19))


# Load ggplot2
library(ggplot2)
library(dplyr)

# Sort your data by PC1 for a smooth line
plot_data <- analysis_data_agg %>%
  mutate(predicted_prob = 1 / (1 + exp(-(mean(alpha_draws) + mean(beta_draws) * dim_1_mean)))) %>%
  arrange(dim_1_mean)

# Create the plot
ggplot(plot_data, aes(x = dim_1_mean)) +
  geom_line(aes(y = predicted_prob), color = "black", size = 1.2) +
  geom_point(aes(y = success), color = "red", size = 3) +
  labs(
    title = "Predicted Probability of Nesting vs PC1",
    x = "PC1 (dim_1_mean)",
    y = "Predicted Probability / Observed Outcome"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 12)
  ) +
  scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(0, 1))



posterior <- as.array(fit) #extract posterior
lp <- log_posterior(fit) #extract log posterior
np <- nuts_params(fit) #extract NUTS parameters (mainly for divergences)

#Trace plots
# Extract parameter names
param_names <- dimnames(posterior)$parameters
# Filter parameters to include only those with "alpha", "beta", or "lp_" in their names
selected_params <- grep("alpha|beta|lp_", param_names, value = TRUE)

# Subset the posterior samples to include only the selected parameters
posterior_filtered <- posterior[, , selected_params, drop = FALSE]

# Create trace plots for the selected parameters
trace_plot <- mcmc_trace(posterior_filtered) + 
  xlab("Post-warmup iteration")

#Rhat plot
summary_fit <- summary(fit)
rhats <- summary_fit[[1]][, "Rhat"]
rhat_plot <- mcmc_rhat(rhats)+ yaxis_text(hjust = 1)

#Effective sample size
ratios <- neff_ratio(fit)
ess_plot <- mcmc_neff(ratios, size = 2)+ yaxis_text(hjust = 1)

#Autocorrelation
autocorr_plot <- mcmc_acf(posterior_filtered,lags = 10)

#Density plots
dens_plot <- mcmc_dens(fit)



###Plots

library(ggplot2)
library(dplyr)
library(tidyr)


# Use `turtle_data_save` for variable data and `turt_centred` for dimension data
dim1_data <- turt_centred %>% select(dim1 = dim_1)

significant_vars <- loadings_selected %>%
  filter(abs(Dim.1) >0.3)

significant_vars <- rownames(significant_vars)

# Combine dimension data with variable data
real_data_combined <- cbind(dim1_data, turtle_data_save %>%
                              select(all_of(significant_vars), location, type))

# Pivot longer for plotting
real_data_long <- real_data_combined %>%
  #mutate(type = ifelse(type == "Preferred", 1, 0))%>%
  pivot_longer(cols = -c(dim1, type, location), names_to = "variable", values_to = "value") %>%
  mutate(variable = gsub("slope", "Slope (deg.)", variable),
         variable = gsub("moisture", "Moisture Content (%)", variable),
         variable = gsub("distance_m", "Distance (m)", variable),
         variable = gsub("total_nest_depth", "Total Nest Depth (cm)", variable),
         variable = gsub("total_water_depth", "Total Water Depth (cm)", variable),
         variable = gsub("width", "Nest Width (cm)", variable),
         variable = gsub("ph", "pH", variable),
         variable = gsub("1_18_", "% 1.18", variable),
         variable = gsub("150_μ", "% 150μm", variable),
         variable = gsub("2_36_", "% 2.36", variable))

# Plot variable impacts with actual data
variable_impact_plot <- ggplot(real_data_long, aes(x = dim1, y = value)) +
  geom_point(aes(shape = as.factor(type), colour = location),  alpha = 0.6) +  # Different shapes for types
  geom_smooth(aes(group = variable), colour = "grey", linewidth = 0.5, method = "loess", se = FALSE) +  # Add a smoothing line per variable
  #stat_ellipse(aes(fill = as.factor(location), group = as.factor(location)),
              # geom = "polygon", size = 0.5, alpha = 0.2, level = 0.95) +  # Add ellipses for location categories
  facet_wrap(
    ~variable, 
    ncol = 1, 
    scales = "free_y", 
    labeller = labeller(variable = label_wrap_gen(width = 15)),  # Wrap labels at 15 characters
    strip.position = "left"
  ) +
  labs(x = "Dimension 1 (Centered Scale)", y = NULL) +  # Remove joint y-axis label
  theme_minimal() +
  theme(
    panel.spacing = unit(1, "lines"),  # Increase spacing between panels
    axis.text.y = element_text(size = 8),  # Adjust y-axis text size for clarity
    strip.placement = "outside",  # Place facet labels outside the axis
    strip.text.y.left = element_text(size = 10, face = "bold", hjust = 0.5),  # Align facet labels to the left
    panel.grid.major.y = element_line(color = "grey90"),  # Optional: Enhance grid visibility
    legend.position = "bottom"
    ) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) +  # Limit number of breaks for better spacing
  scale_color_discrete(name = "") +  # Remove color legend title
  scale_shape_discrete(name = "") +  # Remove shape legend title
  scale_fill_discrete(name = "")  # Remove fill legend title for stat_ellipse


# Extract posterior samples
posterior_samples <- rstan::extract(fit)
beta_samples <- posterior_samples$beta[, 1] # Select the coefficient for dimension 1
alpha_samples <- posterior_samples$alpha

# Calculate predicted probabilities for a range of dimension 1 values
dim1_range <- seq(min(real_data_combined$dim1), max(real_data_combined$dim1), length.out = 100)
pred_probs_matrix <- sapply(dim1_range, function(dim1_val) {
  linear_predictor <- alpha_samples + beta_samples * dim1_val
  plogis(linear_predictor)
})

pred_probs_mean <- colMeans(pred_probs_matrix)
pred_probs_ci <- apply(pred_probs_matrix, 2, quantile, probs = c(0.025, 0.975))

# Create a data frame for probability plotting
prob_data <- data.frame(dim1 = dim1_range, prob_mean = pred_probs_mean, prob_lower = pred_probs_ci[1,], prob_upper = pred_probs_ci[2,])

# Plot probability vs. dimension 1
probability_plot <- ggplot(prob_data, aes(x = dim1, y = prob_mean)) +
  geom_ribbon(aes(ymin = prob_lower, ymax = prob_upper), fill = "lightblue", alpha = 0.5) +
  geom_line(color = "blue") +
  labs(x = "", y = "Predicted Probability \n of Nesting") +
  theme_minimal() +
  scale_y_continuous(labels = scales::percent_format())


# Combine both plots using patchwork
combined_plot <- probability_plot / variable_impact_plot + plot_layout(heights = c(1, 4))

# Print the combined plot
print(combined_plot)


# TaBle of variables

summary_table <- turtle_data_save %>% 
  select(type, location,  slope, moisture, distance_m, total_nest_depth, 
         total_water_depth, width, ph, contains(c("mm", "μ"))) %>%
  group_by(type, location) %>%
  summarise(
    across(.cols = c("slope", "moisture", "distance_m", "total_nest_depth", "total_water_depth", "width", "ph", contains(c("mm", "μ"))), 
           .fns = list(mean = ~mean(.x, na.rm = TRUE), sd = ~sd(.x, na.rm = TRUE))
    ) %>%
      mutate(across(.cols = contains("mean"), .fns = ~round(., 2)),
             across(.cols = contains("sd"), .fns = ~round(., 2)))
  ) %>% t() %>%
  as.data.frame()

rownames(summary_table) <-snakecase::to_title_case(gsub("_", " ", rownames(summary_table)))
rownames(summary_table) <- gsub("Sd", "SD", rownames(summary_table))
rownames(summary_table) <- gsub("1 1", "1.1", rownames(summary_table))
rownames(summary_table) <- gsub("2 3", "2.3", rownames(summary_table))
rownames(summary_table) <- gsub("2 7", "2.7", rownames(summary_table))
rownames(summary_table) <- gsub("6 3", "6.3", rownames(summary_table))



turtle_data_long <- turtle_data_save %>%
  select(type, slope, moisture, distance_m, total_nest_depth, total_water_depth, width, ph, contains(c("mm", "μ"))) %>%
  pivot_longer(cols = -type, names_to = "variable", values_to = "value")

# Create the facetted dot plot
ggplot(turtle_data_long, aes(x = type, y = value)) +
  geom_point() +
  facet_wrap(~ variable, scales = "free_y", ncol = 4) +  # Use free_y to allow different scales in y axis and adjust ncol as needed
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Type", y = "Value", title = "Facetted Dot Plot of Variables by Type")



# Export combined plot to a file (optional)
# ggsave("combined_plot.png", combined_plot, width = 20, height = 12, dpi = 300)
# 
# # Export combined plot to a file (optional)
# ggsave("regression_results.pdf", p, width = 10, height = 6, dpi = 300)
# 
# rownames(correlation_matrix_all_start) <- gsub("_", " ", rownames(correlation_matrix_all_start))
# rownames(correlation_matrix_all_start) <- clean_variable_name(rownames(correlation_matrix_all_start))
# 
# colnames(correlation_matrix_all_start) <- gsub("_", " ", colnames(correlation_matrix_all_start))
# colnames(correlation_matrix_all_start) <- clean_variable_name(colnames(correlation_matrix_all_start))
# 
# 
# # Export combined plot to a file (optional)
# png("correlation_plot_all.png", width = 3200, height = 3200, res = 300)
# corrplot(correlation_matrix_all_start, method = "color", 
#          type = "upper",    # Show only upper triangle
#          addCoef.col = "black",  # Add coefficient labels in black
#          tl.col = "black", tl.srt = 45, # Text label color and rotation
#          diag = FALSE)  # Hide the diagonal
# dev.off()

