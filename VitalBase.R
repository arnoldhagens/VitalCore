# Vital Model Base Model
#
# An R script to solve ODE's of an SIR based model with age structure and nursing homes
# Influenza
# updated January 2024 for
#               - waning immunity
#               - seasonal R
#               - vacination system
#               - Asymptomatic
#               - Symptomatic Home, Outpatient and Inpatient
#               - Nursing homes dynamics
#
# by:          Arnold Hagens
#              Groningen University
#              arnoldjjhagens@gmail.com, a.j.j.hagens@umcg.nl
#
#
#
# Base package contains three files
# 1. VitalBase.R                    R code to run
# 2. SIRageVQSIR.R                  R code functions
# 3. VitalSimulationBase.xlsx       Excel input parameters

library("readxl")
source("SIRageVQSIR.R")

# load parameters
SimFile <- "VitalSimulationBase.xlsx"    #simulation file
parameters <- as.data.frame(read_excel(SimFile, sheet = paste("Parameters_",DiseaseID=1,sep="")))  # read parameter
BaselineParam(parameters,OSA=c(0,0))  # load parameters
betaCalc()                            # calculate beta

# Start simulation
vparameters = c(gamma=gamma,beta=beta,C=C,CFR=CFR,Uptake=Uptake,VacStart=VacStart,YearVacStart=YearVacStart,StepSize=StepSize)
vt = seq(0,Sim_duration/StepSize,1)
inits = c(Sv=Sv_0,Iv=Iv_0,Av=Av_0,Vhm=Vhm_0,Vin=Vin_0,Vicu=Vicu_0,Rv=Rv_0,S=S_0,I=I_0,An=An_0,Nhm=Nhm_0,Nin=Nin_0,Nicu=Nicu_0,R=R_0,D=D_0,V=V_0,Icum=Icum_0,QALY=QALY_0, HC=HC_0, ProdLoss=ProdLoss_0,NHdeaths=NHdeaths_0,NHImmune=NHImmune_0)
Simresults = as.data.frame(rk4(inits, vt,  calculate_derivativVac, vparameters))

# store simulation results
write.csv(Simresults,"Modelresults.csv")


