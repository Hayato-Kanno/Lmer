
# Function to derive LMER
calc_LMER <- function(Pars) {
  
  # Extract parameters from Pars
  MK <- Pars$MK
  Lm <- Pars$Lm
  Steepness <- Pars$Steepness
  SRR <- Pars$SRR
  FM <- Pars$FM
  
  # Setting initial value for LMER
  init_LMER <- Lm + (1 - Lm) / 2 
  # Internal function to calculate unexploited biomass per recruit larger than the given length, assuming that weight is proportional to the cube of length
  calc_B0 <- function(L, MK) {
    return((1 - L)^MK * (1 - 3 * (1 - L)^1 / (1 + 1 / MK)
                           + 3 * (1 - L)^2 / (1 + 2 / MK)
                           - 1 * (1 - L)^3 / (1 + 3 / MK)))
  }
  # Spawning stock biomass per recruit at maximum excess recruitment
  SBMER <- switch(SRR,
                  "BH" = (sqrt(4 * Steepness / (1 - Steepness)) - 1) * 
                    (1 - Steepness) / (5 * Steepness - 1),
                  "RI" = (1 - sqrt((1 - Steepness) / (4 * Steepness))) / 
                    log(4 * Steepness / (1 - Steepness)),
                  "HS" = 1 - Steepness,
                  stop("Invalid SRR. Must be \"BH\" or \"RI\" or \"HS\"")
  )
  # Deriving LMER
  result <- optim(
    par = init_LMER,
    fn = function(L) {
      SB_Finf <- 1 - calc_B0(L, MK) / calc_B0(Lm, MK)  # SB/SB0 for F=inf and Lc = L
      
      return(abs(SB_Finf - SBMER))  # Objective function
    },
    method = "Brent",  # One-dimensional non-linear optimization algorithm with bounds
    lower = Lm,
    upper = 1
  )
  # Check for convergence
  if (result$convergence != 0) stop(result$message)
  
  # Return the LMER
  LMER <- result$par
}
