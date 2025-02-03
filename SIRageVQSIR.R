# Vital Base Model functions
# An R script to solve ODE's of an SIR model with age structure
# updated January 2023 - February 2025 for
#               - waning immunity
#               - seasonal R
#               - vacination system
#               - Asymptomatic
#               - Symptomatic Home, Outpatient and Inpatient
#               - Nursing homes dynamics
#               - periodical population growth
#
# by:          Arnold Hagens
#              Groningen University
#              a.j.j.hagens@umcg.nl
#
#
library("deSolve")


# Global variables for cache Nursing Homes
Sirbeta <- 0
Sirgamma <- 0
SirS0 <- 0
SirI0 <- 0
SirR0 <- 0
Sirtimes <- seq(0:100)
Sirout <- 0

# Global variables (some variables are global to not pass them to functions)
nage    <- 10      # number of age groups
tc      <- 0      # counter for screen disply progress
C       <- matrix(0, nrow = nage, ncol = nage)
beta    <- 0
gamma   <- 0
wii     <- 0
Wiv     <- 0
V_0     <- rep(0, 10)
QALY_0  <- rep(0, 10)
HC_0    <- rep(0, 10)
ProdLoss_0 <- rep(0, 10)
Icum_0  <- rep(0, 10)
NHdeaths_0 <- 0
NHImmune_0 <- 0
AgeGroups <- c()
Prev_t <-0


calculate_derivativVac <-function(t, x, vparameters){

  if(Nursing_enabled==1) {DoNH <- TRUE} else {
    DoNH<-FALSE
    CFR <- CFR_withNH}          # adjust CFR rate to include elderly in NH

  # Vaccinated part
  Sv    = as.matrix(x[1:nage])                # Susceptible
  Iv    = as.matrix(x[(nage+1):(2*nage)])     # infectious
  Av    = as.matrix(x[(2*nage+1):(3*nage)])   # Asymptomatic
  Vhm    = as.matrix(x[(3*nage+1):(4*nage)])  # Vaccinated sick home
  Vin    = as.matrix(x[(4*nage+1):(5*nage)])  # Vaccinated sick outpatient
  Vicu    = as.matrix(x[(5*nage+1):(6*nage)]) # Vaccinated sick inpatient
  Rv    = as.matrix(x[(6*nage+1):(7*nage)])   # Recovered vaccinated

  # non vaccinated part
  S    = as.matrix(x[(7*nage+1):(8*nage)])
  I   = as.matrix(x[(8*nage+1):(9*nage)])
  An    = as.matrix(x[(9*nage+1):(10*nage)])
  Nhm  =  as.matrix(x[(10*nage+1):(11*nage)])
  Nin =  as.matrix(x[(11*nage+1):(12*nage)])
  Nicu =  as.matrix(x[(12*nage+1):(13*nage)])
  R    = as.matrix(x[(13*nage+1):(14*nage)])

  # common death stock
  D = as.matrix(x[(14*nage+1):(15*nage)]) # Deaths
  V = as.matrix(x[(15*nage+1):(16*nage)]) # vaccinated
  Icum = as.matrix(x[(16*nage+1):(17*nage)]) # Cum Infected

  # stocks for CEA
  QALY     = as.matrix(x[(17*nage+1):(18*nage)])
  HC       = as.matrix(x[(18*nage+1):(19*nage)])
  ProdLoss = as.matrix(x[(19*nage+1):(20*nage)])

  ## Nursing homes [total numbers]
  NHdeaths       = x[(20*nage+1)]  # deaths in nursing hom,e
  NHImmune       = x[(20*nage+2)]  # immune nursing homes, attack will lead to nothing....

  # negative value protection, especially required with high step sizes or agressive diseases.
  I[I<0] = 0
  Iv[Iv<0] =0

  S[S<0] = 0
  Sv[Sv<0] =0

  with(as.list(vparameters),{
    # note that because S, I and R are all vectors of length nage, so will N,
    # and dS, dI, and dR

    # helper variables

    Nalive  <- Sv+Iv+Av+Vhm+Vin+Vicu+Rv + S+I+An+Nhm+Nin+Nicu+R       # total individuals alive
    N       <- Nalive + D   # total N

    # periodical booster vaccination moment
    real_t <- t*StepSize                           # the real time due to step size
    nn     <- rep(0,10)
    dVavS  <- rep(0,10)
    dVavR  <- rep(0,10)
    dVavRv <- rep(0,10)
    dVavAv <- rep(0,10)
    dVavAn <- rep(0,10)

    dt   <- as.Date(real_t, origin = StartDate)

    Yr   <- floor(as.numeric(format(dt,'%Y')))     # real year count

    if (Yr>= YearVacStart) {
      dt1 <- as.Date(paste(YearVacStart, "-1-1",sep=""))
      Vac_t <-  length(seq(from = dt1, to = dt, by = 'day'))-1

      for (n in seq(1, nage)) {
        #define if vaccination is occurring at this time
        if (VacStart[n] <= Vac_t &
            VacStart[n] + Vac_speed + 1 > Vac_t) {
          nn[n] <- 1
        }

        # define next vaccination start
        tn <- (VacStart[n] + Vac_speed + 1)
        if (tn < Vac_t) {
          VacStart[n] <<- VacStart[n] + Vac_Interval * YearDays
        }
      }


    # Define shifts from S to Sv (number of vaccinated persons)
      VC <- (1-Uptake)
      VC <-  1-(VC^(1/(Vac_speed)))   #  day vaccination period

      dVavS <- S*VC*nn        # part that moves per day from S to Sv (vaccinated persons)

      dVavR <- R*VC*nn        # part that is vaccianted per day from R to Sv (vaccinated persons) no move
      dVavRv <- Rv*VC*nn        # part that vaccinated per day from Rv to Sv (vaccinated persons) no move
      dVavAv <- Av*VC*nn        # part that vaccinated per day from Av to Sv (vaccinated persons) no move
      dVavAn <- An*VC*nn        # part that vaccinated per day from An to Sv (vaccinated persons) no move

    }

    dV <- dVavS + dVavR + dVavRv + dVavAv + dVavAn    # log the vaccinations sunk stock (cumulative vaccinations)

    # no moves due to vaccination (only costs) since natueral acquired immunity is preferred
    dVavR <- 0
    dVavRv <- 0
    dVavAn <- 0
    dVavAv <-0

    #cospi(180/365*2)+1

    # seasonal effect on R and beta value after callibration
    s      <- cospi(real_t/YearDays*2)+1   # AH changed from -sinpi  to adjust for yearly variations
    betaS  <- beta * (s*SeasonA + SeasonB)

    #Demographic Change per year for a day
    DemoChange <- c(0+(NetMigrationRate[1]-MortalityRate[1])*Nalive[1]-Nalive[1]/Pop_groupwidth[1],
                    Nalive[1]/Pop_groupwidth[1]-Nalive[2]/Pop_groupwidth[2]+(NetMigrationRate[2]-MortalityRate[2])*Nalive[2],
                    Nalive[2]/Pop_groupwidth[2]-Nalive[3]/Pop_groupwidth[3]+(NetMigrationRate[3]-MortalityRate[3])*Nalive[3],
                    Nalive[3]/Pop_groupwidth[3]-Nalive[4]/Pop_groupwidth[4]+(NetMigrationRate[4]-MortalityRate[4])*Nalive[4],
                    Nalive[4]/Pop_groupwidth[4]-Nalive[5]/Pop_groupwidth[5]+(NetMigrationRate[5]-MortalityRate[5])*Nalive[5],
                    Nalive[5]/Pop_groupwidth[5]-Nalive[6]/Pop_groupwidth[6]+(NetMigrationRate[6]-MortalityRate[6])*Nalive[6],
                    Nalive[6]/Pop_groupwidth[6]-Nalive[7]/Pop_groupwidth[7]+(NetMigrationRate[7]-MortalityRate[7])*Nalive[7],
                    Nalive[7]/Pop_groupwidth[7]-Nalive[8]/Pop_groupwidth[8]+(NetMigrationRate[8]-MortalityRate[8])*Nalive[8],
                    Nalive[8]/Pop_groupwidth[8]-Nalive[9]/Pop_groupwidth[9]+(NetMigrationRate[9]-MortalityRate[9])*Nalive[9],
                    Nalive[9]/Pop_groupwidth[9]+(NetMigrationRate[10]-MortalityRate[10])*Nalive[10])


    # by day
    DemoRatio <- rep(0, 10)
    BirthsSusceptible <- rep(0, 10)

    if (Pop_growth_step == 0) {
        DemoChange <- DemoChange / YearDays
        DemoRatio  <- ((DemoChange + Nalive) / Nalive) - 1
        BirthsSusceptible <-
          c(BirthRate * sum(Nalive) / YearDays, 0, 0, 0, 0, 0, 0, 0, 0, 0)   # new individuals arrive only to Susceptible compartment
    } else
    {
      tf <- floor(t / Pop_growth_step)
      if (tf != Prev_t) {
        DemoChange <- DemoChange / YearDays
        DemoRatio  <- ((DemoChange + Nalive) / Nalive) - 1
        DemoRatio <- (1 + DemoRatio) ^ Pop_growth_step - 1

        BD <- (1 + BirthRate / YearDays) ^ Pop_growth_step - 1
        BirthsSusceptible <-
          c(BD * sum(Nalive), 0, 0, 0, 0, 0, 0, 0, 0, 0)   # new individuals arrive only to Susceptible compartment
        Prev_t <<- tf
      }
    }


    # Differentials

    # Incidence from Susceptible and Susceptible Vaccinated
    Inc_Sv <- as.matrix((1-Vac_Inf_Eff)*Sv*betaS)*(as.matrix(C)%*%as.matrix(((1-Vac_Inf_Eff)*Iv+I)/Nalive))
    Inc_S  <- as.matrix(S*betaS)*(as.matrix(C)%*%as.matrix((I+(1-Vac_Inf_Eff)*Iv)/Nalive))

    # Forensics 1 ####
    #new_row <- data.frame(time=t, betaS = betaS)
    #write.table(new_row, "forensics.csv", sep = ",", col.names = FALSE, row.names = FALSE, append = TRUE)

    dSv   = StepSize * ((-Inc_Sv
                                + dVavS + dVavR +dVavRv +dVavAn+dVavAv- Wiv*Sv) +
                          Sv * DemoRatio)

    dIcum = StepSize * (Inc_Sv + Inc_S)

    dIv   = StepSize *((Inc_Sv - gamma*as.matrix(Iv))+ Iv *DemoRatio)

    dAv   = StepSize * (
      (gamma*Iv * (1-Frac_Sympt) - 1/Duration_Asympt*Av - dVavAv)
      + Av * DemoRatio)

    dVhm  = gamma*Iv * Frac_Sympt *
      (Frac_Home + (1-Frac_Home)*Vac_Hosp_Eff) -
      Vhm/(Duration_Home*(1-Vac_Los_Eff))
       + Vhm * DemoRatio

    dVin = StepSize * ((gamma*Iv * Frac_Sympt * Frac_In *(1-Vac_Hosp_Eff)
            - Vin/(Duration_In*(1-Vac_Los_Eff)))+ Vin * DemoRatio)
    dVicu  = StepSize * ((gamma*Iv * Frac_Sympt * Frac_Icu *(1-Vac_Hosp_Eff)
              - Vicu/(Duration_Icu*(1-Vac_Los_Eff))) + Vicu * DemoRatio)


    dRv = StepSize * ((
      (1/(Duration_In*(1-Vac_Los_Eff))*Vin)*(1-CFR*(1-Vac_CFR_Eff))+
      (1/(Duration_Icu*(1-Vac_Los_Eff))*Vicu) *(1-CFR*(1-Vac_CFR_Eff)) +
      (1/(Duration_Home*(1-Vac_Los_Eff))*Vhm) *(1-CFR*(1-Vac_CFR_Eff)) +
      (1/Duration_Asympt*Av) -
      Wii*Rv -
      dVavRv
      ) +
      Rv * DemoRatio)    # Demographic change


    # differentials non-vaccinated
    dS    = StepSize *
                ((-Inc_S -
                dVavS + Wii*(R+Rv) + Wiv*Sv) +
                S * DemoRatio +BirthsSusceptible)
    dI    = StepSize *
                ((Inc_S -
               gamma*as.matrix(I)) + I * DemoRatio)

    dAn   = StepSize *((gamma*I * (1-Frac_Sympt) - An/Duration_Asympt - dVavAn)+An * DemoRatio)

    dNhm  = StepSize *((gamma*I * Frac_Sympt * Frac_Home - Nhm/Duration_Home) + Nhm * DemoRatio)
    dNin = StepSize *((gamma*I * Frac_Sympt * Frac_In- Nin/Duration_In) + Nin * DemoRatio)
    dNicu  = StepSize * ((gamma*I * Frac_Sympt * Frac_Icu- Nicu/Duration_Icu) + Nicu * DemoRatio)

    dR    = StepSize * ((
      (1/Duration_In*Nin)*(1-CFR)+
        (1/Duration_Icu*Nicu)*(1-CFR)+
        (1/Duration_Home*Nhm)*(1-CFR) +
        (1/Duration_Asympt*An)
      -Wii*R-dVavR) +
        R * DemoRatio)      # Demographic change

    # Nursing home dynamic ####
    if (DoNH==TRUE) {
      Iday  <- sum(I)      # infectious non vac
      Ivday <- sum(Iv)     # infectious vac

      # number of infectious persons that might visit to nursing homes
      NinfecPeople <-  Ivday +  Iday

      # R0 is the reproduction rate for naive and unvaccinated population
      # vaccination assumes that vaccinated cannot get infected according to effectiveness average of agegroup 8,9,10
      R0adjusted  <- Nursing_R0 * (1-Nursing_uptake* mean(Vac_Inf_Eff[8:10]))   # adjust R0 for uptake NursingR0*(1-uptake*VacEff)
      #cat("R0_adjusted: ",R0adjusted,"\n")

      # Forensics 1 ####
      #new_row <- data.frame(time=t, R0adjusted = R0adjusted)
      #write.table(new_row, "forensicsNH.csv", sep = ",", col.names = FALSE, row.names = FALSE, append = TRUE)

      # get Nursing attack results deaths and immune centers
      Nat <- NH_Attack2(R0adjusted,Nursing_CFR,Nursing_visitfraction,Ninfec =  NinfecPeople,Nalive=Nalive,NHImmune)

      dNHdeaths <-  Nat$deaths  # creating periodical deaths
      dNHImmune <-  (Nat$NewImmune - NHImmune / Nursing_LoI)

      # compensate background deaths and deaths by disease with new habitants from S,Sv, R, Rv
      # all Nursing homes are assumed to be always full
      BackgroundDeaths <-  Nursing_number*Nursing_persons*Nursing_BGRate/YearDays   #including background deaths nursing homes, these are not in the normal demographics
      AllNHdeaths <-  dNHdeaths + BackgroundDeaths
      Age9tot <- S[9] + Sv[9] + R[9] + Rv[9]

      dS[9]  <- dS[9] - S[9]/Age9tot * AllNHdeaths
      dSv[9] <- dSv[9] - Sv[9]/Age9tot * AllNHdeaths
      dR[9]  <- dR[9] - R[9]/Age9tot * AllNHdeaths
      dRv[9] <- dRv[9] - Rv[9]/Age9tot * AllNHdeaths

    }else
    {dNHdeaths <- 0
    dNHImmune <-0
      }

    ### Reporting indicators (deaths, NHdeaths, QALYs, costs) only calculated when in reporting period

    if (real_t > Sim_duration-Report_duration) {
      # deaths
      # Single death stock for all age groups and Nursing homes
      dD <-   StepSize * (
        (1/(Duration_In*(1-Vac_Los_Eff))*Vin)*(CFR*(1-Vac_CFR_Eff)) +
          (1/(Duration_Icu*(1-Vac_Los_Eff))*Vicu) *(CFR*(1-Vac_CFR_Eff)) +
          (1/(Duration_Home*(1-Vac_Los_Eff))*Vhm) *(CFR*(1-Vac_CFR_Eff)) +
          (Nin/Duration_In)*(CFR)+
          (Nicu/Duration_Icu)*(CFR) +
          (Nhm/Duration_Home)*(CFR)
      )
            # add NH deaths over agegroup 8,9,10 (as they also move out of the susceptible groups)
      dD[9] <- dD[9] + Nalive[9]/(Nalive[9]+Nalive[10])* dNHdeaths
      dD[10] <- dD[10] + Nalive[10]/(Nalive[9]+Nalive[10])* dNHdeaths

      DayInflation <- ((1+Inflation)^(1/YearDays)) -1
      CostInflator <- (1+DayInflation)^real_t

      # add QALYs
      dQALY <- StepSize * ((
        (S+Sv+I+Iv+An+Av+R+Rv)*Utility_Healthy/YearDays  +
          (Nhm+Vhm)*Utility_Rel_Home*Utility_Healthy/YearDays +
          (Nin+Vin) *Utility_Rel_Inpatient*Utility_Healthy/YearDays +
          (Nicu+Vicu) *Utility_Rel_ICU*Utility_Healthy/YearDays
      )/(1+Re)^(real_t/YearDays))

      # assign NH QALYS over agegroup 8,9,10
      n <- ((Nursing_persons * Nursing_number) * Nursing_utility)/YearDays
      dQALY[9] <-
        dQALY[9] + ((Nalive[9] / ( Nalive[9]+ Nalive[10])) * n) / (1 + Re) ^ (real_t / YearDays)
      dQALY[10] <-
        dQALY[10] + ((Nalive[10] / (Nalive[9]+ Nalive[10])) * n) / (1 + Re) ^ (real_t / YearDays)

      # Health care costs (+ vaccination cost)
      dHC  <- StepSize * (
        (CostInflator*(
          (Nhm+Vhm) * Cost_home +
          (Nin+Vin) * Cost_inpatient +
          (Nicu+Vicu) * Cost_ICU +
          (dV) * Cost_Vac)
        )
        /(1+Rf)^(real_t/YearDays))

      # Productivity loss
      dProdLoss  <- StepSize *(
        CostInflator*(
          (Nhm+Vhm) * ProdLoss_home +
          (Nin+Vin) * ProdLoss_inpatient +
          (Nicu+Vicu) * ProdLoss_ICU +
            dD * ProdLossDeath
      )/
        (1+Rf)^(real_t/YearDays))
    } else {
      #cat ("outside")
      dQALY <- rep(0,nage)
      dHC <- rep(0,nage)
      dProdLoss <- rep(0,nage)
      dD <- rep(0,nage)
      dNHdeaths <- 0
    }

    out=c(dSv,dIv,dAv,dVhm,dVin,dVicu,dRv,dS,dI,dAn,dNhm,dNin,dNicu,dR,dD,dV,dIcum,dQALY,dHC,dProdLoss,dNHdeaths,dNHImmune)
    list(out)
  })
}

BaselineParam <- function(parameters,OSA=c(0,0,0), Scenario=0){
  cat("Parse parameters","\n")

  lbl <- parameters[,1]
  description <-parameters[,2]
  ParaLabels <<- data.frame(lbl,description)

  for(i in 3:123) {                    # assign function within loop
    NewValue <- parameters[i,4]
    for(j in 5:13) {
       if (!(is.na(parameters[i,j]))) {
        NewValue <- append(NewValue,parameters[i,j])

      }
    }

    if(i==OSA[1]) {
      if (OSA[3]>0) {
        n1 <- rep(1,10)
        n1[OSA[3]] <- OSA[2]    # only change specific age group
        n <- n1
      }else
      {n <- OSA[2]  }        # change all parameters of age group

      NewValue=as.numeric(NewValue)*n #OSA[2]
      cat("  New value",NewValue,"\n")
      }   # for one way sensitivity analysis

    if (Scenario>0) {  #if scenario
      scen <- 0
      try(scen <- as.numeric(parameters[i,3]))
      if(is.na(scen)){scen<-0}

      if (scen==Scenario) {  # only if parameter matches for specific scenario
        assign(paste0(parameters[i,1]), as.numeric(NewValue), envir = .GlobalEnv)   # assign values
        cat(" New param", parameters[i,1],NewValue,"\n")
      }
      } else{
      assign(paste0(parameters[i,1]), as.numeric(NewValue), envir = .GlobalEnv)   # assign values
      }

  }

  # Error adjustments if parameters are beyond acceptable limits
  Uptake[Uptake>=1] <<- 0.99
  Uptake[Uptake<=0] <<- 0.01

  AgeGroups <<- c("0-9","10-19","20-29","30-39","40-49","50-59","60-69","70-79","80-89","90+")

  nage   <<- length(Pop_frac)
  C      <<- matrix(c(C1,C2,C3,C4,C5,C6,C7,C8,C9,C10), ncol = 10)   # built contact matrix
  N      <<- Pop_size * Pop_frac
  gamma  <<- 1/inf_period
  Wii    <<- 1/LoIm_Inf
  Wiv    <<- 1/LoIm_Vac
  V_0    <<- rep(0,10)
  NHdeaths_0 <<- 0
  NHImmune_0 <<- 0
  StartDate <<- as.Date(StartDate,origin = "1900-01-01")

  betaCalc()
}

betaCalc <- function(){
  # calculate beta and reverse engineer beta from the R0 and gamma
  nage   <<- length(Pop_frac)
  M   <- C
  for (i in 1:nage){
    for (j in 1:nage) {
      M[i,j] = C[i,j]*Pop_frac[i]/Pop_frac[j]
    }
  }

  # reverse engineer beta from the R0 and gamma
  eig = eigen(M)
  beta <<- R0*gamma/max(Re(eig$values))
  cat(" beta  ",beta,"\n")
  cat(" gamma ",gamma,"\n")
  cat(" eigen ",max(Re(eig$values)),"\n")
  cat(" R0    ",R0,"\n")

  Wii    <<- 1/LoIm_Inf
  Wiv    <<- 1/LoIm_Vac

}

#NH_Attack2(1.15,0.33,0.002,100000,0,1000)
NH_Attack2 <- function(R0, cfr, VisitFraction, Ninfec, Nalive,Immune) {
  # persons visiting homes
  Ninfec <- as.integer(Ninfec * VisitFraction)

  # number homes visited
  NHv <- NH_visited(Ninfec, Nursing_number)
  #cat("NHv : " , NHv, "\n")


  # new susceptible number nursing homes attacked
  nAttacked <- NHv * (1 - Immune / Nursing_number)
  ARn <- AR(R0, inf_period, I0 = 1, Nursing_persons)

  totalInfected <- nAttacked *
    ARn *    # attack rate estimate
    Nursing_persons                           # total infected of all nursing homes

  deaths <- totalInfected * cfr                # total deaths


  if (deaths < 0) {
    deaths = 0
  }
  NH <- list(deaths = deaths, NewImmune = nAttacked)

  return(NH)
}

NH_visited <- function(visits,NHn,Reg=TRUE){
  # estimates the number of uniquely visited nursing homes depending on the number of visitors to a number
  # random nursing homes
  if (Reg==TRUE){
    # 2nd degree polynomal estimate
    b1=-0.309064998654273
    b2=0.931230691993834
    b3=0.00555750457116663
    nv <- visits/NHn
    return(NHn*(b1*nv^2+b2*nv+b3))
  }else{

    if(visits/NHn<8) {
      visits <- as.integer(visits)
      cnt <-
        function(x) {
          length(unique(sample(1:NHn, visits, replace = T)))
        }
      b <- data.frame(lapply(seq(1:1000), cnt))

      return(mean(as.numeric(b[1, 1:1000])))
    } else{
      return(1)
    }
  }

}

ForRegression <- function(){
  NH_visited(620, 500, FALSE)
  NH_visited(620, 500, TRUE)
  df <- data.frame(x = as.numeric(), y = as.numeric())

  for (x in seq(0:500)) {
    print(x)
    df[x, ] <- c(x, NH_visited(x, 500))
  }

}



AR <- function(R0,duration,I0,N) {
  # estimates the attack rate (number of infected persons)
  Results <- sir(R0/duration,1/duration,N-I0,I0,0,seq(1:120))
  n <- nrow(Results)
  AR <- (Results$R[n]+Results$I[n])/N
  return(AR)
}

sir <- function(bet, gamm, S0, I0, R0, times) {
  #cachecheck
  if (Sirbeta==bet & Sirgamma==gamm & SirS0==S0 & SirI0==I0 & SirR0==R0){
    return(Sirout)
  } else {

    # standard SIR model
    require(deSolve) # for the "ode" function

    sir_equations <- function(time, variables, parameters) {
      with(as.list(c(variables, parameters)), {
        N = S+I+R
        dS <- -bet * I * S/N
        dI <-  bet * I * S/N - gamm * I
        dR <-  gamm * I

        #dS <- -bet * I * S^2/N^2
        #dI <-  bet * I * S^2/N^2 - gamm * I
        #dR <-  gamm * I


        return(list(c(dS, dI, dR)))
      })
    }

    # the parameters values:
    parameters_values <- c(bet  = bet, gamm = gamm)

    # the initial values of variables:
    initial_values <- c(S = S0, I = I0, R = R0)

    # solving
    out <- ode(initial_values, times, sir_equations, parameters_values)

    # for cache
    Sirbeta <<- bet
    Sirgamma <<- gamm
    SirS0 <<- S0
    SirI0 <<- I0
    SirR0 <<- R0
    Sirtimes <<- times
    Sirout <<- as.data.frame(out)

    # returning the output:
    as.data.frame(out)
  }
}
