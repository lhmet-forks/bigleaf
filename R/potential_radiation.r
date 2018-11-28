
#' extraterrestrial solar radiation
#'
#' Compute the extraterrestrial solar radiation with the
## eccentricity correction.
#'
#' @param doy integer vector with day of year (DoY)
#'
#' @section Details:
#' Computation follows Lanini, 2010 (Master thesis, Bern University)
#'
#' @return numeric vector of extraterrestrial radiation (W_m-2)
#' @export
#'
#' @examples
#' plot(1:365, computeExtRadiation(1:365), type = "l"
#'   , ylab = "radiation (W m-2)", xlab = "day of year")
computeExtRadiation <- function(doy) {
  # Fractional year in radians
  FracYearRad <- 2 * pi * (doy - 1) / 365.24
  # Total solar irradiance
  SolarIrr_Wm2.c <- 1366.1 #W / m-2
  #Eccentricity correction
  ExtRadiation.V.n <- SolarIrr_Wm2.c * (
    1.00011 + 0.034221 * cos(FracYearRad) + 0.00128 * sin(FracYearRad)
     + 0.000719 * cos(2 * FracYearRad) + 0.000077 * sin(2 * FracYearRad)
     )
  ExtRadiation.V.n
}

#' potential radiation
#'
#' Compute potential radiation for given geolocation and day of yeat
#'
#' @param doy integer vector with day of year (start at 1),
#'   same length as hour or length 1
#' @param hour numeric vector with daytime as decimal hour of local time zone
#' @param latDeg Latitude in (decimal) degrees
#' @param longDeg Longitude in (decimal) degrees
#' @param timezone Time zone (in hours)
#' @param useSolartime by default corrects hour (given in local winter time)
#'   for latitude to solar time (where noon is exactly at 12:00).
#'   Set this to FALSE to directly use local winter time
#' @return vector of potential radiation (PotRad, W_m-2)
#' @export
#' @importFrom solartime computeSunPositionDoyHour
#' @examples
#' hour <- seq(5, 18, by = 0.1)
#' potRadApparentLocal <- computePotentialRadiation(
#'   160, hour, 39.94, -5.77, timezone = +1)
#' potRadTimezone <- computePotentialRadiation(
#'   160, hour, 39.94, -5.77, timezone = +1, useSolartime = FALSE)
#' plot(potRadApparentLocal ~ hour, type = 'l'
#'   , ylab = 'potential radiation (W m-2)')
#' lines(potRadTimezone ~  hour, col = "blue")
#' abline(v = 12, col = "blue", lty = "dotted")
#' legend("bottomright", legend = c("solar time", "local winter time")
#' , col = c("black", "blue"), inset = 0.05, lty = 1)
computePotentialRadiation <- function(
  doy, hour, latDeg, longDeg, timezone, useSolartime = TRUE
) {
  # Calculate potential radiation from solar elevation and extraterrestrial
  # solar radiation
  solElevRad <- computeSunPositionDoyHour(
    doy, hour, latDeg, longDeg, timezone
    , isCorrectSolartime = useSolartime)[,"elevation"]
  extRadiation <- computeExtRadiation(doy)
  potRad <- ifelse(
    solElevRad <= 0, 0, extRadiation * sin(solElevRad) )
  potRad
}
