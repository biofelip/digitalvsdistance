library(dsims)
# Find the file path to the example shapefile in dssd
shapefile.name <- system.file("extdata", "StAndrew.shp", package = "dssd")
# Create the survey region object
region <- make.region(
    region.name = "St Andrews bay",
    shape = shapefile.name,
    units = "m"
)
plot(region)

# defning population
# We first create a flat density grid
density <- make.density(
    region = region,
    x.space = 500,
    constant = 1
)

# Now we can add some high and low points to give some spatial variability
density <- add.hotspot(
    object = density,
    centre = c(-170000, 6255000),
    sigma = 8000,
    amplitude = 4
)

density <- add.hotspot(
    object = density,
    centre = c(-160000, 6275000),
    sigma = 6000,
    amplitude = 4
)

density <- add.hotspot(
    object = density,
    centre = c(-155000, 6260000),
    sigma = 3000,
    amplitude = 2
)

density <- add.hotspot(
    object = density,
    centre = c(-150000, 6240000),
    sigma = 10000,
    amplitude = -0.9
)

density <- add.hotspot(
    object = density,
    centre = c(-155000, 6285000),
    sigma = 10000,
    amplitude = -1
)

# I will choose to plot in km rather than m (scale = 0.001)
plot(density, region, scale = 0.001)


# add covariates
# Create a covariate list describing the distribution of cluster sizes
covariates <- list(
    cluster_size = list(distribution = "ztruncpois", mean = 3),
    animal_size = list(distribution = "normal", mean = 150, sd = 50),
    training_score = list(distribution = "lognormal", meanlog = 0, sdlog = 1),
    sea_state_score = list(distribution = "normal", mean = 50, sd = 30)
)


# Define the population description
pop.desc <- make.population.description(
    region = region,
    density = density,
    covariates = covariates,
    N = 1000,
    fixed.N = TRUE
)

# defining the survey design
parallel.design <- make.design(
    region = region,
    design = "systematic",
    spacing = 2500,
    edge.protocol = "minus",
    design.angle = 90,
    truncation = 1200
)

zigzag.design <- make.design(
    region = region,
    design = "eszigzag",
    spacing = 2233,
    edge.protocol = "minus",
    design.angle = 0,
    bounding.shape = "convex.hull",
    truncation = 1200
)
p.survey <- generate.transects(parallel.design)
plot(region, p.survey)
z.survey <- generate.transects(zigzag.design)
plot(region, z.survey)


# make detectability
# Define the covariate parameters on the log scale
cov.param <- list(
    cluster_size = log(1.08),
    animal_size = log(2.5),
    training_score = log(5),
    sea_state_score = log(10)
)

# Create the detectability description
detect <- make.detectability(
    key.function = "uf",
    scale.param = 1,
    shape.param = 1,
    cov.param = cov.param,
    truncation = 1200
)

# Plot the simulation detection functions
plot(detect, pop.desc)


analyses <- make.ds.analysis(
    dfmodel = list(~1, ~ cluster_size + animal_size + training_score + sea_state_score),
    key = c("hr", "hr"),
    truncation = 1200,
    er.var = "R2",
    criteria = "AIC"
)


# simulate that shit
sim.parallel <- make.simulation(
    reps = 100,
    design = parallel.design,
    population.description = pop.desc,
    detectability = detect,
    ds.analysis = analyses
)

sim.zigzag <- make.simulation(
    reps = 100,
    design = zigzag.design,
    population.description = pop.desc,
    detectability = detect,
    ds.analysis = analyses
)

# Generate a single instance of a survey: a population, set of transects
# and the resulting distance data
eg.parallel.survey <- run.survey(sim.parallel)

# # Plot it to view a summary
plot(eg.parallel.survey, region)


# # run that shi?
# # Running the simulations
# sim.parallel <- run.simulation(sim.parallel)
# sim.zigzag <- run.simulation(sim.zigzag)

# summary(sim.parallel, usa.max.reps = TRUE)
