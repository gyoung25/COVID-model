# COVID-model
Some code used in Young, Xiao, Newcomb, Michael 2021

covidParamVaried.m is a MATLAB file to simulate our covid model with baseline parameters, then simulates it 1000 more times, each time varying model parameters by up to 10% baseline value (uniformly). The code then plots the baseline infected-class solution curve along with two more curves, between which 80% of all simulates fall, thereby giving a representative range of outcomes

model_Vax_ac.ode is our baseline model written for XPPAUT
Rt_BVP_ac.ode and Ic_BVP_ac.ode are both boundary value problems used to study how certain quantities change as model parameters are varied. Rt_BVP_ac.ode tracks changes in the effective reproduction ratio and Ic_BVP_ac.ode tracks the peak (maximum) ICU caseload. Both are written for XPPAUT 
