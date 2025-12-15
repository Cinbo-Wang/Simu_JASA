This repository contains the  R code to reproduce the main simulation studies from the paper “Covariate Adjustment in Randomized Clinical Trials with GLMs: On the Role of Jackknife and Practical Recommendations.”

**File Overview**

• `Main.R`  

  Defines the simulation settings for both the main text and supplementary material. It compares the following estimation methods:

  • Traditional estimators: $\hat{\tau}_{\texttt{OLS}}$, $\hat{\tau}_{\texttt{gOB}}$, and its calibrated version $\hat{\tau}_{\texttt{gOB{-}cal}}$.

  • Standard leave-one-out (LOO) estimators: $\hat{\tau}_{\texttt{LOO}}$, along with two calibrated versions $\hat{\tau}_{\texttt{LOO-cal}}$ and $\hat{\tau}_{\texttt{LOO-cal-out}}$.

  • U-statistic-based estimators: $\hat{\tau}_{\texttt{adj,2}}$, $\hat{\tau}_{\texttt{adj,2}}^{\dag}$, $\hat{\tau}_{\texttt{adj,2-cal}}$, $\hat{\tau}_{\texttt{adj,2-cal}}^{\dag}$.

  • Newly proposed estimators: $\hat{\tau}_{\texttt{JASA}}$ and $\hat{\tau}_{\texttt{JASA-cal}}$.

• `fun_base.R`  

  Contains the function to compute the traditional estimators: $\hat{\tau}_{OLS}$, $\hat{\tau}_{gob}$, and its calibrated counterpart $\hat{\tau}_{gbcal}$.

• `fun_LeaveOneOut.R` 

  Implements functions for calculating the standard LOO estimators: $\hat{\tau}_{loo}$, $\hat{\tau}_{loo-cal}$, and $\hat{\tau}_{loo-cal-out}$.

• `fun_utils-simudata.R`

  Provides the utility function to generate the simulation data used across all experiments.