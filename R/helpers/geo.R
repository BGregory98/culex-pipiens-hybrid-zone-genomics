# =============================================================================
# R/helpers/geo.R
# =============================================================================
#
# Great-circle distance and the two city centres used for the distance-to-city
# covariate.
#
#   source("R/helpers/geo.R")
#   d_km <- hav(lat1, lon1, lat2, lon2)
# =============================================================================

# Haversine distance in km (Earth radius 6,371 km).
hav <- function(lat1, lon1, lat2, lon2) {
  R <- 6371; p <- pi / 180
  a <- sin((lat2 - lat1) * p / 2)^2 +
       cos(lat1 * p) * cos(lat2 * p) * sin((lon2 - lon1) * p / 2)^2
  2 * R * asin(pmin(1, sqrt(a)))
}

# City centres (latitude, longitude).
BAL <- c(39.2904, -76.6122)   # Baltimore
DC  <- c(38.9072, -77.0369)   # Washington
