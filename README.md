# Timeseries_Analysis_GMO_applicaitons_au
Timeseries analysis and modelling of gmo applications in Australia over the last 15 years


# Summary
As the global biotechnology sector rapidly expands, Australia’s contribution to the industry has so far been considered amongst the top in the world. Therefore, the purpose of this project is to investigate the past trends of Australia’s biotechnology sector in order to find a suitable time series model that could be used to forecast its future trend. This investigation is undertaken by examining the number of monthly Notifiable Low Risk Dealing (NLRD) applications that have been made to the Office of the Gene Technology Regulator between the years 2003 to 2019. 
Both deterministic and stochastic time series modelling methods were utilized to analyze the number of NLRD applications to determine the most suitable model for forecasting. The results indicated that the stochastic seasonal autoregressive integrated moving average (SARIMA) models were more suitable for this particular dataset. In particular, SARIMA(2,1,3)X(1,0,0)12 model was chosen as the most suitable, and was consequently used to forecast ten months ahead of the dataset. The forecasted result indicated that there is a flat trend stabilization in the number of NLRD applications within the period. 

# Coding
Coding is in R.   Run finalcode.R to load the data and produce a list of all the graphs used in the report

# Report
Report is saved as GMO application gmo timeseries Report- Final.docx
