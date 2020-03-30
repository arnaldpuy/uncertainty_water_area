## ----setup, include=FALSE-----------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ----load_packages, results="hide", message=FALSE, warning=FALSE--------

# LOAD PACKAGES ---------------------------------------------------------------

# Function to read in all required packages in one go:
loadPackages <- function(x) {
  for(i in x) {
    if(!require(i, character.only = TRUE)) {
      install.packages(i, dependencies = TRUE)
      library(i, character.only = TRUE)
    }
  }
}

loadPackages(c("data.table", "ggplot2", "sensobol", "scales",
               "ncdf4", "rworldmap", "sp", "countrycode", 
               "dplyr", "IDPmisc", "boot", "parallel", 
               "MASS", "doParallel", "complmrob", 
               "mvoutlier", "sandwich", "lmtest", "mice", 
               "ggridges", "broom", "naniar", "cowplot", 
               "tidyr", "benchmarkme"))

# SET CHECKPOINT --------------------------------------------------------------

dir.create(".checkpoint")

library("checkpoint")

checkpoint("2019-09-10", 
           R.version ="3.6.1", 
           checkpointLocation = getwd())

# CUSTOM FUNCTION TO DEFINE THE PLOT THEMES -----------------------------------

theme_AP <- function() {
  theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.background = element_rect(fill = "transparent",
                                           color = NA),
          legend.key = element_rect(fill = "transparent",
                                    color = NA))
}


## ----functions_data, cache=TRUE-----------------------------------------

# CREATE FUNCTIONS ------------------------------------------------------------

# Function to obtain UN code, Continent and Country names
country_code <- function(dt) {
  dt[, `:=` (Code = countrycode(dt[, Country], 
                                origin = "country.name", 
                                destination = "un"), 
             Continent = countrycode(dt[, Country], 
                                     origin = "country.name", 
                                     destination = "continent"))]
  dt[, Country:= countrycode(dt[, Code], 
                             origin = "un", 
                             destination = "country.name")]
  setcolorder(dt, c("Country", "Continent", "Code", "Water.Withdrawn"))
  return(dt)
}

## Function to transform longitude and latitude to country.
# It is borrowed from Andy:
# https://stackoverflow.com/questions/14334970/convert-latitude-and-longitude-coordinates-to-country-name-in-r)
coords2country = function(points) {  
  countriesSP <- rworldmap::getMap(resolution = 'low')
  pointsSP = sp::SpatialPoints(points, proj4string=CRS(proj4string(countriesSP)))  
  indices = sp::over(pointsSP, countriesSP)
  indices$ADMIN  
  #indices$ISO3 # returns the ISO3 code 
  #indices$continent   # returns the continent (6 continent model)
  #indices$REGION   # returns the continent (7 continent model)
}

# Function to load and extract data from .nc files
get_nc_data <- function(nc_file) {
  nc <- nc_open(nc_file)
  ww <- ncvar_get(nc, "withd_irr")
  lon <- ncvar_get(nc, "lon")
  lat <- ncvar_get(nc, "lat")
  water <- rowSums(ww[, 469:ncol(ww)]) # Obtain year values for 2010 only
  ww.df <- data.frame(cbind(lon, lat, water)) 
  countries <- coords2country(ww.df[1:nrow(ww.df), 1:2])
  df <- cbind(countries, ww.df)
  setDT(df)
  final <- df[, .(Water.Withdrawn = sum(water)), countries]
  setnames(final, "countries", "Country")
  country_code(final)
  out <- na.omit(final[order(Continent)])
  out[, Water.Withdrawn:= Water.Withdrawn / 1000] # From mm to m
  return(out)
}


## ----water_with_dataset, cache=TRUE, warning=FALSE----------------------

# READ IN DATASETS ON IRRIGATION WATER WITHDRAWAL -----------------------------

# FAO data (Table 4) ----------------------------
# http://www.fao.org/nr/water/aquastat/water_use_agr/IrrigationWaterUse.pdf

# UNIT IS KM3/YEAR
table4.tmp <- fread("table_4.csv", skip = 3, nrows = 167) %>%
  .[, .(Country, Year, Water.withdrawal)] %>%
  setnames(., "Water.withdrawal", "Water.Withdrawn")

# Extract the selected years
table4.dt <- country_code(table4.tmp[Year %in% 1999:2012])[
  , Water.Dataset:= "Table 4"][
    , Year:= NULL]

# Liu et al. dataset ----------------------------
#UNIT IS 10^9 m3/year = km3/year
liu.dt <- fread("liu.csv")[, .(country, irr)] %>%
  setnames(., c("country", "irr"), c("Country", "Water.Withdrawn")) %>%
  country_code(.) %>%
  .[, Water.Dataset:= "Liu et al. 2016"]

# Huang et al datasets --------------------------
names_nc_files <- c("withd_irr_lpjml.nc", "withd_irr_pcrglobwb.nc", 
                    "withd_irr_h08.nc", "withd_irr_watergap.nc")
out.nc <- lapply(names_nc_files, function(x) get_nc_data(x))
names(out.nc) <- c("LPJmL", "PCR-GLOBWB", "H08", "WaterGap")

GHM.dt <- rbindlist(out.nc, idcol = "Water.Dataset") %>%
  .[order(Country)]


## ----arrange_total_countries, cache=TRUE, dependson="water_with_dataset"----

# ARRANGE THE TOTAL NUMBER OF COUNTRIES ---------------------------------------

# Read in list of countries in UN
DT <- fread("UN_countries2.csv", select = "Country")

# Give standard country names, UN codes and link with Continent
DT <- DT[, `:=` (Code = countrycode(DT[, Country], 
                                    origin = "country.name", 
                                    destination = "un"), 
                 Continent = countrycode(DT[, Country], 
                                         origin = "country.name", 
                                         destination ="continent"))]

# Manually add Micronesia
DT <- DT[, Continent:= ifelse(Country %in% "Micronesia", "Oceania", Continent)]


## ----final_water_dataset, cache=TRUE, dependson=c("water_with_dataset" ,"arrange_total_countries")----

# CREATE THE FINAL IRRIGATION WATER WITHDRAWAL DATASET ------------------------

# Merge with the Country vector
tmp <- GHM.dt[, merge(.SD, DT, by = c("Country", "Code", "Continent"), 
                      all.y = TRUE), by = Water.Dataset]

# Check whether there are duplicated Countries
tmp[tmp[, duplicated(Country), Water.Dataset][, V1]]

# Get mean values for the duplicated Countries
GHM.dt.full <- tmp[, .(Water.Withdrawn = mean(Water.Withdrawn)), 
                   .(Water.Dataset, Country, Code, Continent)]

# Arrange Liu data set
liu.dt.full <- merge(DT, liu.dt, 
                     by = c("Country", "Code", "Continent"), 
                     all.x = TRUE) %>%
  .[, Water.Dataset:= ifelse(is.na(Water.Dataset), 
                             "Liu et al. 2016", 
                             Water.Dataset)]

# Arrange Table 4 dataset
table4.dt.full <- merge(DT, table4.dt, 
                        by = c("Country", "Code", "Continent"), 
                        all.x = TRUE) %>%
  .[, Water.Dataset:= ifelse(is.na(Water.Dataset), 
                             "Table 4", 
                             Water.Dataset)]

# Obtain final irrigation water withdrawal dataset
water.dt <- rbind(liu.dt.full, table4.dt.full, GHM.dt.full)


## ----area_dataset, cache=TRUE-------------------------------------------

# READ IN IRRIGATED AREA DATASETS ---------------------------------------------

meier.dt <- fread("meier.csv") %>%
  setnames(., "Codes", "Code")


## ----merge_with_area, cache=TRUE, dependson=c("area_dataset", "water_with_dataset", "final_water_dataset")----

# MERGE DATASETS --------------------------------------------------------------

irrigated.area.datasets <- colnames(meier.dt)[-c(1:3)]

irrigated.dt <- melt(meier.dt, measure.vars = irrigated.area.datasets) %>%
  .[, merge(.SD, DT, by = c("Country", "Code", "Continent"), 
            all.y = TRUE), by = variable] %>%
  setnames(., c("variable", "value"), c("Area.Dataset", "Irrigated.Area"))

tmp.dt <- merge(irrigated.dt, water.dt, on = c("Continent", "Country", "Code"), 
                allow.cartesian = TRUE) %>%
  .[!Continent == "Oceania"] # Drop Oceania

# Vector with the countries to drop
countries.drop <- tmp.dt[Water.Withdrawn == 0 & is.na(Irrigated.Area) == TRUE] %>%
  .[, unique(Country)]

# Drop the countries
full.dt <- tmp.dt[!Country %in% countries.drop]


## ----export.irrigated.dt, cache=TRUE, dependson="merge_with_area"-------

# EXPORT IRRIGATED AREA DT ----------------------------------------------------

fwrite(irrigated.dt, "irrigated.dt.csv")


## ----plot_merged, cache=TRUE, dependson="merge_with_area", dev="tikz", fig.height=8, fig.width=7, fig.cap="Scatter plots of irrigated areas (x-axis) against irrigation water withdrawal (y-axis). All the possible combinations of datasets are shown."----

# PLOT ------------------------------------------------------------------------

full.dt %>%
  ggplot(., aes(Irrigated.Area, Water.Withdrawn, 
                color = Continent)) +
  geom_point(size = 0.8) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10 ^ (2 * x)),
                labels = trans_format("log10", math_format(10 ^ .x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10 ^ x),
                labels = trans_format("log10", math_format(10 ^ .x))) +
  labs(x = "Irrigated area (Mha)", 
       y = expression(paste("Water withdrawal", "", 10^9, m^3/year))) +
  facet_grid(Water.Dataset ~ Area.Dataset) +
  theme_AP() +
  theme(legend.position = "top", 
        strip.text = element_text(size = 7))


## ----log10, cache=TRUE, dependson="merge_with_area"---------------------

# TRANSFORM DATASET -----------------------------------------------------------

cols <- c("Water.Withdrawn", "Irrigated.Area")
col_names <- c("Continent", "Water.Dataset", "Area.Dataset", "Regression", 
               "Imputation.Method", "Iteration")
cols_group <- c("Continent", "Area.Dataset", "Water.Dataset")
full.dt <- full.dt[, (cols):= lapply(.SD, log10), .SDcols = (cols)]


## ----export_dataset_log10, cache=TRUE, dependson="log10"----------------

# EXPORT FULL DATASET WITH MISSING VALUES --------------------------------------------

fwrite(full.dt, "full.dt.csv")


## ----plot_missing, cache=TRUE, dependson="log10", dev="tikz", fig.height=8, fig.width=7----

# SCATTERPLOT SHOWING MISSING VALUES ------------------------------------------

scatter.na <- copy(full.dt)
cols_transform <- c("Irrigated.Area", "Water.Withdrawn")
scatter.na[, (cols_transform):= lapply(.SD, function(x) 10 ^ x), .SDcols = cols_transform]
tmp <- split(scatter.na, scatter.na$Continent)

gg <- list()
for(i in names(tmp)) {
  gg[[i]] <- ggplot(tmp[[i]], aes(Irrigated.Area, Water.Withdrawn)) +
    geom_miss_point() +
    facet_grid(Water.Dataset ~ Area.Dataset) + 
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10 ^ (2 * x)),
                labels = trans_format("log10", math_format(10 ^ .x))) +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10 ^ (2 * x)),
                  labels = trans_format("log10", math_format(10 ^ .x))) +
    labs(x = "Irrigated area (Mha)", 
         y = expression(paste("Water withdrawal", "", 10^9, m^3/year))) +
    theme_AP() +
    theme(legend.position = "top",
          strip.text = element_text(size = 7)) +
    ggtitle(names(tmp[i]))
}

gg


## ----plot_missing2, cache=TRUE, dependson="log10", fig.height=8, fig.width=4----

# PLOT PERCENTAGE OF MISSING --------------------------------------------------

na.dt <- full.dt[, .(Continent, Country, Irrigated.Area, Water.Withdrawn)]
tmp <- split(na.dt[, !"Continent", with = FALSE], na.dt$Continent)

gg <- list()
for(i in names(tmp)) {
  gg[[i]] <- gg_miss_fct(x = tmp[[i]], fct = Country) +
    coord_flip() +
    labs(x = "Parameter", 
         y = "") +
    viridis::scale_fill_viridis(name = "Missing") +
    theme(legend.position = "top") +
    ggtitle(names(tmp[i]))
}

gg


## ----outliers, cache=TRUE, fig.keep="none", echo=FALSE,results='hide', dependson="log10"----

# CHECK OUTLIERS -------------------------------------------------------------

tmp <- NaRV.omit(full.dt)[, dd.plot(.SD), .SDcols = cols, cols_group]


## ----plot_outliers, cache=TRUE, dependson="outliers", dev="tikz", fig.height=8.5, fig.width=6.5, fig.cap="Bivariate outliers."----

# PLOT OUTLIERS ---------------------------------------------------------------

ggplot(tmp, aes(md.cla, md.rob, 
                shape = outliers, 
                color = Continent)) +
  geom_point() +
  scale_shape_manual(name = "Outlier", 
                     labels = c("No", "Yes"), 
                     values = c(16, 4)) +
  facet_grid(Water.Dataset ~ Area.Dataset) +
  labs(x = "Mahalanobis distance", 
       y = "Robust distance") +
  theme_AP() +
  theme(legend.position = "top", 
        strip.text = element_text(size = 7))


## ----hetero, cache=TRUE, dependson="merge_with_area"--------------------

# CHECK HETEROSKEDASTICITY ----------------------------------------------------

hetero <- NaRV.omit(full.dt)[, .(model = .(lm(Water.Withdrawn ~ Irrigated.Area))), cols_group]


## ----hetero_plot, cache=TRUE, dependson="hetero", dev="tikz", fig.height=8.5, fig.width=6.5, fig.cap="Residual versus fitted values."----

# PLOT HETEROSKEDASTICITY ----------------------------------------------------

tmp <- hetero[, c("resid" = lapply(model, residuals), 
                  "pred" = lapply(model, fitted)), cols_group]

hetero.plot <- split(tmp, tmp$Continent)

gg <- list()
for(i in names(hetero.plot)) {
  gg[[i]] <- ggplot(hetero.plot[[i]], aes(pred, resid)) +
    geom_point() +
    facet_grid(Area.Dataset ~ Water.Dataset) +
    theme_AP() +
    labs(x = "Fitted", 
         y = "Residual") +
    geom_hline(yintercept = 0, 
               lty = 2) +
    ggtitle(names(hetero.plot[i]))
}

gg


## ----missing_values, cache=TRUE, dependson="log10", echo=FALSE, message=FALSE, results='hide', warning=FALSE----

# IMPUTATION OF MISSING VALUES -------------------------------------------------

# Substitute Inf values for NA
for (j in 1:ncol(full.dt)) set(full.dt, which(is.infinite(full.dt[[j]])), j, NA)

full.dt[, lapply(.SD, function(x) sum(is.infinite(x)))] # Check

# Imputation settings
m.iterations <- 20
imputation.methods <- c("rf", "norm.boot", "norm.nob")
set.seed(10) # Set seed to allow for replication

# Run
imput <- full.dt[, .(Group = lapply(imputation.methods, function(x) 
  mice(.SD, m = m.iterations, maxit = m.iterations, method = x, seed = 500, 
       print = FALSE))), 
  cols_group]

imput <- imput[, Imputation.Method:= rep(imputation.methods, .N / length(imputation.methods))]

# Extract iterations
imput <- imput[, Datasets:= lapply(Group, function(x) 
  lapply(1:m.iterations, function(y) data.table(mice::complete(x, y))))] %>%
  .[, Data:= lapply(Datasets, function(x) rbindlist(x, idcol = "Iteration"))]

# Vector to loop onto
columns_add <- c("Country", "Iteration", "Irrigated.Area", "Water.Withdrawn")
tmp <- as.list(columns_add)
names(tmp) <- columns_add

# Extract columns
for(i in names(tmp)) {
  imput <- imput[, tmp[[i]]:= lapply(.SD, function(x) 
    lapply(x, function(y) y[, ..i])), .SDcols = "Data"]
}

# Unlist
full.imput <- imput[, lapply(.SD, unlist), 
                    .SDcols = columns_add, 
                    .(Continent, Area.Dataset, Water.Dataset, Imputation.Method)]


## ----conduct_lm, cache=TRUE, dependson="missing_values"-----------------

# LINEAR REGRESSIONS AND PULL RESIDUALS ---------------------------------------

# Set grouping vectors
cols_group <- c("Continent", "Area.Dataset", "Water.Dataset")
all.cols <- c(cols_group, "Imputation.Method", "Iteration")

dt <- full.imput[, .(fit = list(list(lm(Water.Withdrawn ~ Irrigated.Area)))), all.cols] %>%
  .[, lapply(fit, function(x) lapply(x, function(y) augment(y))), all.cols]

resid.dt <- dt[, lapply(V1, function(x) unlist(x$.resid)), all.cols] %>%
  .[, Country:= full.imput$Country] %>%
  .[, Method:= "Normal"]

dtR <- full.imput[, .(fit = list(list(rlm(Water.Withdrawn ~ Irrigated.Area, 
                                          maxit = 60)))), all.cols] %>%
  .[, lapply(fit, function(x) lapply(x, function(y) augment(y))), all.cols]


resid.dtR <- dtR[, lapply(V1, function(x) unlist(x$.resid)), all.cols] %>%
  .[, Country:= full.imput$Country] %>%
  .[, Method:= "Robust"]

final.resid <- rbind(resid.dt, resid.dtR)


## ----plot_residuals, cache=TRUE, dependson="conduct_lm", fig.height=8, fig.width=6, warning=FALSE----

 # PLOT RESIDUALS -------------------------------------------------------------

# function to plot
gg_bar <- gg_ridge <- list()
for(i in c("Africa", "Americas", "Asia", "Europe")) {
  country.vector <- final.resid[Continent == i] %>%
    .[order(V1)] %>%
    .$Country %>%
    unique(.)
  gg_ridge[[i]] <- final.resid[Continent == i] %>%
    .[, Country:= factor(Country, levels = country.vector)] %>%
    ggplot(., aes(V1, Country, fill = ..x..)) +
    geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01) +
    geom_vline(xintercept = 0, 
               lty = 2) +
    scale_fill_gradient2(name = "Residual", 
                         low = "blue", mid = "white",
                         high = "red", midpoint = 0) +
    labs(x = "Residual", 
         y = "") +
    theme_AP() + 
    theme(legend.position = "top")
  gg_bar[[i]] <- final.resid[Continent == i] %>%
    .[, sum(V1 > 0) / .N, Country] %>%
    ggplot(., aes(reorder(Country, V1), V1, 
                  fill = ..y..)) + 
    geom_bar(stat = "identity") +
    coord_flip() + 
    labs(y = "Residuals > 0 (%)", 
         x = "") +
    theme_AP() + 
    theme(legend.position = "top")
}

gg_ridge


## ----export_dataset_imputation, cache=TRUE, dependson="missing_values"----

# EXPORT FULL DATASET WITHOUT MISSING VALUES --------------------------------------------

fwrite(full.imput, "full.imput.csv")



## ----compute_beta, cache=TRUE, dependson="missing_values"---------------

# COMPUTE BETA AND R^2 --------------------------------------------------------

regressions <- full.imput[, .(Normal = coef(lm(Water.Withdrawn ~ Irrigated.Area)), 
                              Robust = coef(rlm(Water.Withdrawn ~ Irrigated.Area, maxit = 60)), 
                              r.squared = summary(lm(Water.Withdrawn ~ Irrigated.Area))$r.squared), 
                      all.cols] 

results <- regressions[, Type:= rep(c("Intercept", "Beta"), times = nrow(.SD) / 2)] %>%
  melt(., measure.vars = c("Normal", "Robust"), variable.name = "Regression") %>%
  tidyr::spread(., Type, value) %>%
  .[, index:= paste(Continent, Area.Dataset, Water.Dataset, Regression, 
                   Imputation.Method, Iteration, sep = "_")]



## ----export_beta_results, cache=TRUE, dependson="compute_beta"----------

# EXPORT BETA AND R^2 RESULTS -------------------------------------------------

results <- regressions[, Type:= rep(c("Intercept", "Beta"), times = nrow(.SD) / 2)] %>%
  melt(., measure.vars = c("Normal", "Robust"), variable.name = "Regression") %>%
  spread(., Type, value) %>%
  .[, index:= paste(Continent, Area.Dataset, Water.Dataset, Regression, 
                   Imputation.Method, Iteration, sep = "_")]

fwrite(results, "results.csv")


## ----lookup, cache=TRUE, dependson="compute_beta"-----------------------

# CREATE LOOKUP TABLE ---------------------------------------------------------

lookup <- setkey(results, index)

# Export lookup
fwrite(lookup, "lookup.csv")


## ----set_sample_matrix, cache=TRUE--------------------------------------

# DEFINE THE SETTINGS OF THE SAMPLE MATRIX ------------------------------------

Continents <- c("Africa", "Americas", "Asia", "Europe")

# Create a vector with the name of the columns
parameters <- paste("X", 1:5, sep = "")

# Select sample size
n <- 2 ^ 13

# Define order
order <- "third"


## ----sample_matrix, cache=TRUE, dependson="set_sample_matrix"-----------

# CREATE THE SAMPLE MATRIX ----------------------------------------------------

# Create an A, B and AB matrices for each continent
sample.matrix <- lapply(Continents, function(Continents) 
  sobol_matrices(N = n,
                 params = parameters,
                 order = order) %>%
    data.table())

# Name the slots, each is a continent
names(sample.matrix) <- Continents

# Name the columns
sample.matrix <- lapply(sample.matrix, setnames, parameters)


## ----transform_sample_matrix, cache=TRUE, dependson=c("sample_matrix", "set_sample_matrix")----

# TRANSFORM THE SAMPLE MATRIX -------------------------------------------------

# Transform the sample matrix
transform_sample_matrix <- function(dt) {
  dt[, X1:= floor(X1 * (6 - 1 + 1)) + 1] %>%
    .[, X1:= ifelse(X1 == 1, "Aquastat", 
                    ifelse(X1 == 2, "FAOSTAT", 
                           ifelse(X1 == 3, "Siebert.et.al.2013", 
                                  ifelse(X1 == 4, "Meier.et.al.2018", 
                                         ifelse(X1 == 5, "Salmon.et.al.2015", 
                                                "Thenkabail.et.al.2009")))))] %>%
    .[, X2:= floor(X2 * (6 - 1 + 1)) + 1] %>%
    .[, X2:= ifelse(X2 == 1, "LPJmL", 
                    ifelse(X2 == 2, "H08", 
                           ifelse(X2 == 3, "PCR-GLOBWB", 
                                  ifelse(X2 == 4, "WaterGap", 
                                         ifelse(X2 == 5, "Table 4", 
                                                "Liu et al. 2016")))))] %>%
    .[, X3:= floor(X3 * (2 - 1 + 1)) + 1] %>%
    .[, X3:= ifelse(X3==1, "Normal", "Robust")] %>%
    .[, X4:= floor(X4 * (length(imputation.methods) - 1 + 1)) + 1] %>%
    .[, X4:= ifelse(X4 == 1, imputation.methods[1], 
                    ifelse(X4 == 2, imputation.methods[2], imputation.methods[3]))] %>%
    .[, X5:= floor(X5 * (m.iterations - 1 + 1)) + 1] 
}

sample.matrix <- lapply(sample.matrix, transform_sample_matrix)
sample.matrix.dt <- rbindlist(sample.matrix, idcol = "Continent")

fwrite(sample.matrix.dt, "sample.matrix.dt.csv")


## ----print_matrix)------------------------------------------------------

# PRINT SAMPLE MATRIX ---------------------------------------------------------

print(sample.matrix.dt)


## ----define_model, cache=TRUE-------------------------------------------

# THE MODEL -------------------------------------------------------------------

model <- function(X) lookup[.(paste0(X[, 1:6], collapse = "_"))][, c(Intercept, Beta, r.squared)]


## ----run_model, cache=TRUE, dependson=c("define_model", "set_boot", "lookup", "transform_sample_matrix")----

# RUN THE MODEL---------------------------------------------------------------

# Set number of cores at 75%
n_cores <- floor(detectCores() * 0.75)

# Create cluster
cl <- makeCluster(n_cores)
registerDoParallel(cl)

# Run model in parallel
Y <- foreach(i=1:nrow(sample.matrix.dt),
             .packages = "data.table") %dopar%
  {
    model(sample.matrix.dt[i])
  }

# Stop parallel cluster
stopCluster(cl)


## ----arrange_output, cache=TRUE, dependson="run_model"------------------

# ARRANGE MODEL OUTPUT -------------------------------------------------------

sample.matrix.dt <- cbind(sample.matrix.dt, data.table(do.call(rbind, Y))) %>%
  setnames(., c("V1", "V2", "V3"), c("Intercept", "Beta", "r.squared"))

# Select the A and B matrix only (for uncertainty analysis)
AB.dt <- sample.matrix.dt[, .SD[1:(n * 2)], Continent]

# Export results
fwrite(sample.matrix.dt, "sample.matrix.dt.csv")
fwrite(AB.dt, "AB.dt.csv")


## ----plot_uncertainty, cache=TRUE, dependson="arrange_output", dev="tikz", fig.height=4, fig.width=5, fig.cap="Uncertainty in the model output. a) Uncertainty in the model goodness of fit. All sources of uncertainty have been taken into account except the selection between OLS and OLS robust (trigger X2). Robust OLS does not allow to compute $r^2$. b) Uncertainty in the slope."----

# PLOT UNCERTAINTY ------------------------------------------------------------

# Plot r2
a <- ggplot(AB.dt, aes(r.squared)) + 
  geom_histogram(color = "black", fill = "white") + 
  theme_AP() +
  labs(x = expression(r ^ 2), 
       y = "Density") +
  scale_y_continuous(breaks = pretty_breaks(n = 3)) +
  scale_x_continuous(breaks = pretty_breaks(n = 3)) +
  facet_wrap(~Continent, ncol = 4) +
  theme(panel.spacing.x = unit(4, "mm"))

# Plot beta
b <- ggplot(AB.dt, aes(Beta)) + 
  geom_histogram(color = "black", fill = "white") + 
  theme_AP() +
  geom_vline(xintercept = 1, 
             lty = 2, 
             color = "red") +
  labs(x = "$\\beta$", 
       y = "") +
  scale_y_continuous(breaks = pretty_breaks(n = 3)) +
  scale_x_continuous(breaks = pretty_breaks(n = 3)) +
  facet_wrap(~Continent, ncol = 4) +
  theme(panel.spacing.x = unit(4, "mm"))


# Merge
plot_grid(a, b, ncol = 1, align = "hv", labels = "auto")


## ----check_order_r2, cache=TRUE, dependson="arrange_output", fig.height=7, fig.width=4.5----

# CHECK THE COMBINATIONS LEADING TO HIGHEST R2 VALUES (first 200) -------------

# Check min and max of first 200 r.squared values
AB.dt[order(-r.squared), head(.SD, 200), Continent][
  , .(min = min(r.squared), max = max(r.squared)), Continent]

# Create plot
dd <- AB.dt[order(-r.squared), head(.SD, 200), Continent] %>%
  .[, ID:= paste(X1, X2, X4, X5, sep = "_")] %>%
  .[, .N, .(Continent, ID)]

tmp <- split(dd, dd$Continent)
gg <- list()
for(i in names(tmp)) {
  gg[[i]] <- ggplot(tmp[[i]], aes(reorder(ID, N), N)) +
    geom_bar(stat = "identity", 
             color = "black", 
             fill = "white") + 
    coord_flip() +
    labs(y = expression(italic(N)), 
         x = "") +
    theme_AP() +
    ggtitle(names(tmp[i]))
}

gg


## ----super_sub, cache=TRUE, dependson="arrange_output"------------------

# SUPERLINEAR OR SUBLINEAR REGIME? --------------------------------------------

AB.dt[, .(sublinear = sum(Beta < 1, na.rm = TRUE) / .N, 
          superlinear = sum(Beta > 1, na.rm = TRUE) / .N), 
      Continent]

AB.dt[, .(">90" = sum(r.squared > 0.9) / .N,
            "75-90" = sum(r.squared > 0.75 & r.squared < 0.9) / .N, 
            "50-75" = sum(r.squared > 0.50 & r.squared < 0.75) / .N, 
            "25-50" = sum(r.squared > 0.25 & r.squared < 0.5) / .N), Continent]


## ----scatterplots, cache=TRUE, dependson="arrange_output", fig.height=5, fig.width=7----

# PLOT SCATTERPLOTS OF PARAMETERS VS MODEL OUTPUT -----------------------------

AB.dt <- AB.dt[, X5:= factor(X5, levels = as.factor(1:m.iterations))]

scatter.dt <- melt(AB.dt, measure.vars = paste("X", 1:5, sep = "")) %>%
  .[, Regime:= ifelse(Beta > 1, "Super-linear", "Sub-linear")]

# Beta
ggplot(scatter.dt, aes(value, Beta, color = Regime)) +
  geom_point(alpha = 0.3, size = 0.5) +
  facet_grid(Continent ~ variable,
             scales = "free_x") +
  geom_hline(yintercept = 1,
             lty = 2) +
  scale_color_manual(values = c("#00BFC4", "#F8766D")) +
  theme_AP() +
  labs(x = "", y = expression(beta)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top")

# R squared
ggplot(scatter.dt, aes(value, r.squared, color = Regime)) +
  geom_point(alpha = 0.3, size = 0.5) +
  facet_grid(Continent ~ variable,
             scales = "free_x") +
  labs(x = "", 
       y = expression(italic(r)^2)) +
  scale_color_manual(values = c("#00BFC4", "#F8766D")) +
  theme_AP() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top")


## ----sensitivity_settings, cache=TRUE-----------------------------------

# SENSITIVITY SETTINGS --------------------------------------------------------

# Number of bootstrap replicas
R <- 1000


## ----sensitivity, cache=TRUE, dependson=c("arrange_output", "set_sample_matrix")----

# SENSITIVITY ANALYSIS --------------------------------------------------------

# Beta
indices <- sample.matrix.dt[, sobol_indices(Y = Beta, 
                                            N = n,
                                            params = parameters, 
                                            first = "jansen",
                                            R = R, 
                                            boot = TRUE,
                                            parallel = "multicore", 
                                            ncpus = n_cores, 
                                            order = "third"), 
                            Continent]

# r squared
indicesR <- sample.matrix.dt[, sobol_indices(Y = r.squared, 
                                             N = n,
                                             params = parameters, 
                                             first = "jansen",
                                             R = R, 
                                             boot = TRUE,
                                             parallel = "multicore", 
                                             ncpus = n_cores, 
                                             order = "third"), 
                             Continent]

# Compute Sobol' indices for the dummy parameter (Beta)
indices.dummy <- sample.matrix.dt[, sobol_dummy(Y = Beta, 
                                                N = n,
                                                params = parameters, 
                                                R = R, 
                                                boot = TRUE,
                                                parallel = "multicore", 
                                                ncpus = n_cores), 
                                  Continent]

# Compute Sobol' indices for the dummy parameter (r squared)
indices.dummyR <- sample.matrix.dt[, sobol_dummy(Y = r.squared, 
                                                 N = n,
                                                 params = parameters, 
                                                 R = R, 
                                                 boot = TRUE,
                                                 parallel = "multicore", 
                                                 ncpus = n_cores), 
                                   Continent]


## ----plot_sobol, cache=TRUE, dependson="ci", dev="tikz", fig.cap="Sobol' indices.", fig.width = 5----

# PLOT SOBOL' INDICES ---------------------------------------------------------

a <- ggplot(indices[sensitivity %in% c("Si", "Ti")], aes(parameters, original, fill = sensitivity)) +
  geom_bar(stat = "identity",
           position = position_dodge(0.6),
           color = "black") +
  geom_errorbar(aes(ymin = low.ci,
                    ymax = high.ci),
                position = position_dodge(0.6)) +
  scale_y_continuous(breaks = pretty_breaks(n = 3)) +
  facet_wrap(~Continent,
             ncol = 4) +
  labs(x = "",
       y = "Sobol' index") +
  scale_fill_discrete(name = "Sobol' indices",
                      labels = c(expression(S[italic(i)]),
                                 expression(T[italic(i)]))) +
  theme_AP() +
  theme(legend.position = "none")


b <- ggplot(indicesR[sensitivity %in% c("Si", "Ti")], aes(parameters, original, fill = sensitivity)) +
  geom_bar(stat = "identity",
           position = position_dodge(0.6),
           color = "black") +
  geom_errorbar(aes(ymin = low.ci,
                    ymax = high.ci),
                position = position_dodge(0.6)) +
  scale_y_continuous(breaks = pretty_breaks(n = 3)) +
  facet_wrap(~Continent,
             ncol = 4) +
  labs(x = "",
       y = "Sobol' index") +
  scale_fill_discrete(name = "Sobol' indices",
                      labels = c(expression(S[italic(i)]),
                                 expression(T[italic(i)]))) +
  theme_AP() +
  theme(legend.position = "none")

legend <- get_legend(a + theme(legend.position = "top"))

all <- plot_grid(a, b, ncol = 1, align = "hv", labels = "auto")

plot_grid(legend, all, ncol = 1,  rel_heights = c(0.1, 1))


## ----sum_si, cache=TRUE, dependson="ci"---------------------------------

# CHECK SUM OF FIRST-ORDER INDICES --------------------------------------------

lapply(list(indices, indicesR), function(x) 
  x[sensitivity == "Si"] %>%
      .[, sum(original), Continent])


## ----plot_sobol_second_third, cache=TRUE, dependson="ci", dev = "tikz", fig.cap="High-order interactions between the triggers."----

# PLOT SOBOL' INDICES (SECOND AND THIRD ORDER) --------------------------------

lapply(c("Sij", "Sijk"), function(x)
  ggplot(indices[sensitivity == x], aes(reorder(parameters, original), original)) + 
  geom_point() + 
  geom_errorbar(aes(ymin = low.ci, 
                    ymax = high.ci)) + 
  facet_wrap(~Continent) +
  theme_bw() + 
  labs(x = "", y = "Sobol' index") + 
  geom_hline(yintercept = 0, lty = 2, color = "red") + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.background = element_rect(fill = "transparent", 
                                         color = NA), 
        legend.key = element_rect(fill = "transparent", 
                                  color = NA), 
        axis.text.x = element_text(angle = 45, 
                                   hjust = 1)))


## ----session_information------------------------------------------------

# SESSION INFORMATION ---------------------------------------------------------

sessionInfo()

## Return the machine CPU
cat("Machine: "); print(get_cpu()$model_name)

## Return number of true cores
cat("Num cores: "); print(detectCores(logical = FALSE))

## Return number of threads
cat("Num threads: "); print(detectCores(logical = TRUE))

## Return the machine RAM
cat("RAM: "); print (get_ram()); cat("\n")


