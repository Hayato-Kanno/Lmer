
# Function to derive length at first capture that achieve MER (Lc_MER) for a given fishing pressure (F/M)
calc_LcMER <- function(Pars) { 
  
  # Extract parameters from Pars
  MK <- Pars$MK
  Lm <- Pars$Lm
  Steepness <- Pars$Steepness
  SRR <- Pars$SRR
  Lc <- Pars$Lc
  FM <- Pars$FM
  
  # Setting initial value for Lc_MER
  init_LcMER <- Lm + (1 - Lm) / 2 # Arbitrarily any value
  
  # Function to calculate yield per recruit with knife-edge Lc and FM
  calc_Y <- function(Lc, MK, FM) {
    factor <- MK * (1 + FM)
    FM / (1 + FM) * (1 - Lc)^MK * (1 - 3 * (1 - Lc)^1 / (1 + 1 / factor)
                                     + 3 * (1 - Lc)^2 / (1 + 2 / factor)
                                     - 1 * (1 - Lc)^3 / (1 + 3 / factor))
  }
  
  # Function to calculate equilibrium recruitment relative to virgin state that arbitrarily set to 1
  calc_Rec <- function(SB, SRR, Steepness) {
    if (SRR == "BH") {
      return(max((4 * Steepness * SB) / ((1 - Steepness) + (5 * Steepness - 1) * SB), 0))
    } else if (SRR == "RI") {
      return(max(SB * (4 * Steepness / (1 - Steepness))^(1 - SB), 0))
    } else if (SRR == "HS") {
      return(min(SB * (1 / (1 - Steepness)), 1))
    } else {
      stop("Invalid SRR: Must be 'BH' or 'RI' or 'HS'.")
    }
  }
  
  # Function to calculate unexploited biomass per recruit larger than the given length
  calc_B0 <- function(L, MK) {
    (1 - L)^MK * (1 - 3 * (1 - L)^1 / (1 + 1 / MK)
                    + 3 * (1 - L)^2 / (1 + 2 / MK)
                    - 1 * (1 - L)^3 / (1 + 3 / MK))
  }
  
  # Function to calculate exploited biomass per recruit larger than the given length with knife-edge Lc and FM
  calc_B <- function(Lc, L, MK, FM) {
    factor <- MK * (1 + FM)
    (1 - Lc)^(-MK * FM) * (1 - L)^factor * 
      (1 - 3 * (1 - L)^1 / (1 + 1 / factor) 
         + 3 * (1 - L)^2 / (1 + 2 / factor)
         - 1 * (1 - L)^3 / (1 + 3 / factor))
  }
  
  # Calculate exploited spawning stock biomass relative to unexploited states (SB/SB0)
  calc_SB <- function(Lc, MK, FM, Lm) {
    B0_Lm <- calc_B0(Lm, MK)
    B0_Lc <- calc_B0(Lc, MK)
    B_Lc <- calc_B(Lc, Lc, MK, FM)
    B_Lm <- calc_B(Lc, Lm, MK, FM)
    
    if (Lc >= Lm) {
      return((B0_Lm - B0_Lc + B_Lc) / B0_Lm)
    } else {
      return(B_Lm / B0_Lm)
    }
  }
  
  # Calculate spawning stock biomass per recruit at maximum excess recruitment (MER)
  SBMER <- switch(SRR,
                  "BH" = (sqrt(4 * Steepness / (1 - Steepness)) - 1) * 
                    (1 - Steepness) / (5 * Steepness - 1),
                  "RI" = (1 - sqrt((1 - Steepness) / (4 * Steepness))) / 
                    log(4 * Steepness / (1 - Steepness)),
                  "HS" = 1 - Steepness,
                  stop("Invalid SRR. Must be \"BH\" or \"RI\" or \"HS\"")
  )
  
  # Objective function for optimization
  objective_function <- function(Lc) {
    SB <- calc_SB(Lc, MK, FM, Lm)
    # Check if SB is NA or infinite
    if (is.na(SB) || is.infinite(SB)) {
      return(0)  # Return 0 if SB is invalid
    }
    return(abs(SB - SBMER))
  }
  
  # Optimize to find the FM that achieves MER for a given Lc
  result <- optim(
    par = init_LcMER,
    fn = objective_function,
    method = "Brent",
    lower = 0,
    upper = 1
  )
  
  # Check for convergence
  if (result$convergence != 0) stop(result$message)
  
  # Return the FM_MER
  LcMER <- result$par
  SB_LcMER <- calc_SB(LcMER, MK, FM, Lm)
  Rec_LcMER <- calc_Rec(SB_LcMER, SRR, Steepness)
  SPR_LcMER <- SB_LcMER / Rec_LcMER
  YPR_LcMER <- calc_Y(LcMER, MK, FM) 
  Yield_LcMER <- YPR_LcMER * Rec_LcMER
  
  list(LcMER = LcMER, Yield_LcMER = Yield_LcMER, SPR_LcMER = SPR_LcMER)
}
