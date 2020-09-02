library(dplyr)  
library(TSA)  
library(stringr)
library(ggplot2)
library(forecast)
library(ggpubr)
library(readxl)
library(lubridate)
library(tseries)
library(urca)
library(x12)
library(tidyverse)
library(lmtest)
library(broom)
library(gridExtra)

#suppress warnings
options(warn=-1)

################################
# Self defined functions #
################################

# Addapted from Ditiatkovski M, 2020
# Custom function to extract summary statistics for a lm
# type model and put them into a single line
summaryTable <- function(model){
  
  residuals <- rstudent(model)
  
  tempSummary <- summary(model)
  
  outputDump <- capture.output(tempSummary)
  adjR <- tempSummary$adj.r.squared
  resErr <- tempSummary$sigma
  
  if(is.null(resErr)){
    resErr <- str_extract(outputDump, "standard error: \\d+.\\d")
    resErr <- str_extract(resErr[!is.na(resErr)],"\\d+.\\d") %>% as.numeric()
  }
  
  # extract p-value from summary dump
  pVal <- str_extract(outputDump, "p-value: \\d+.\\d")
  
  if(sum(!is.na(pVal)) == 0){
    pVal <- str_extract(outputDump, "p-value: < \\d+.\\d+e-\\d+")
    pVal <- str_extract(pVal[!is.na(pVal)], "\\d+.\\d+e-\\d+") %>% as.numeric()
  }else{
    pVal <- str_extract(pVal[!is.na(pVal)], "\\d+.\\d") %>% as.numeric()
  }
  
  # calculate p value for normality
  shw <- shapiro.test(residuals)
  #this does not make sence for for power 1 no intercept this fails if combined with previou line
  shw <- shw$p.value
  resErr <- tempSummary$sigma
  # make a nice summary table
  BIC <- BIC(model) 
  MASE <- dLagM::GoF(model)$MASE
  summaryTable <- tibble("Adjusted R" = adjR, "Residual Error" = resErr,
                         "p value" = pVal, "Shapiro-Wilks p" = shw, "BIC" = BIC, "MASE"= MASE)
  return(summaryTable)
}


# Addapted from Ditiatkovski M, 2020
# Supporting function that ranks model on MASE and BIC. Expects table in format of sumaryTable() output
# Returns the rank table and the model in a list
modelRank <- function(sumTable){
  
  #Add rank for R^2 and residual error
  sumTable <-dplyr::arrange(sumTable, `BIC`) %>% 
    mutate(BICRank = cumsum(index/index))
  
  sumTable <-dplyr::arrange(sumTable, `MASE`) %>%
    mutate(MASERank = cumsum(index/index),
           RankSum = MASERank + BICRank) %>% 
    dplyr::arrange(`RankSum`, `MASE`)
  
  # Check if residuals are normal and p value for model is significant
  norm <- which(sumTable$`Shapiro-Wilks p`>0.05)
  sig <- which(sumTable$`p value`<0.05)
  
  # find the highest ranking model that has normal residuals p<0.05
  # In case of a tie the one with higher adjusted R is selected
  if(length(norm)>0 | length(sig)>0){
    suppressWarnings(modelNo <- min(intersect(norm, sig)))
    if(is.infinite(modelNo)){
      modelNo <- 1
    }
  }
  else{
    modelNo <- 1
  }
  
  # Return updated table and the index of the model in a list
  modelAndTable <- list("table" = sumTable, "selectedModel" = sumTable$index[modelNo])
  
  modelAndTable$table$index <- NULL
  
  return (modelAndTable )
}


# Addapted from Ditiatkovski M, 2020
# Supporting function for Model picking. Writes permutations of formula for lm with permutations
# up to maximum power of 3 (ie cubic). Expects a vector with list of variable names in formula
# sorted in ascending order, max power (int) and optional seasonal/harmonic vector  
formulaList <- function(powers = NULL, maxPower = 1, season = NULL){  
  
  formulaStr <- NULL
  # make the begining of the formula
  if(!is.null(season)){
    
    formulaStr <- paste0("data~", str_c(season, collapse = " + ")) 
  }else{
    formulaStr <- "data ~"
  }
  
  
  formulaList <- NULL
  
  # loop for formula createion (can be combined with later loop but makes it too complex)
  # make this a separate function that returns a list of formula
  for(p1 in 1:2){
    for(p2 in 1:2){
      for(p3 in 1:2){
        for(inter in 0:1){
          # skip if all powers all are 2 (ie season only component)
          if(p3+p2+p1 == 6){
            next() 
          }
          
          #skip if not max power
          if(maxPower == 1 & (p1 == 2 | (p1 == 1 & p3+p2 <= 3))){
            next() 
          }
          if(maxPower == 2 & p3 == 1){
            next() 
          }
          
          
          
          # reset the formula to starting version
          tempFormula <- formulaStr
          # add items as appropriate where 1=true and 2 =false
          if(p1 == 1){
            tempFormula = paste0(tempFormula,"+", powers[1])
          }
          if(p2 == 1){
            tempFormula = paste0(tempFormula,"+", powers[2])
            
          }      
          if(p3 == 1){
            tempFormula = paste0(tempFormula,"+", powers[3])
            
          }
          
          tempFormula = paste0(tempFormula, "+", inter)
          
          # add created formula to a list of formulas  
          formulaList <- c(formulaList, as.formula(tempFormula))
        }
      }
    }
  }
  return(formulaList)
  
}



# Addapted from Ditiatkovski M, 2020
# Function will cycle through combination of models up to power of 3
# and selects one with the highest R^2 and lowest residual error. Expects orignal data, maximum power
# integer that will be constant in all models and seasonal variable.
modelSelect <- function(data, maxPower, season = NULL){
  
  # setupfor formula prep and data managment
  powers <- seq(1:maxPower)
  powerNamesForm <- paste0("Power_", powers)
  
  powerNames <- powerNamesForm
  
  dates <- matrix(time(data), ncol = maxPower, nrow = (length(data)))
  dates <-as_tibble(t(t(dates)^powers), .name_repair = "minimal")
  
  seasonCol <- NULL
  
  
  # deal with seasonal variable if needed
  if(!is.null(season)){
    dimensions <- dim(season)
    
    if(is.null(dimensions)){
      seasonCol <- "season"
    }else{
      seasonCol <- colnames(season)
      seasonCol <- str_replace_all(seasonCol,"\\*|\\(|\\)","")
    }
    dates <- cbind(season, dates)
    powerNames <- c(seasonCol, powerNames)
  }
  # Make the formula list for all models
  formulaList <- formulaList(powerNamesForm, maxPower, season = seasonCol)
  data <- cbind(as.vector(data), dates)
  
  powerNames <- c("data", powerNames)
  
  # create final frame for model fitting
  
  colnames(data) <- powerNames
  data <- as_tibble(data, .name_repair = "minimal")
  
  modelList <- list()
  sumTable <- NULL
  
  for (i in 1:length(formulaList)) {
    modelList[[i]] <- lm(formulaList[[i]], data)
    tempRow <- cbind("index" = i, "Model" = deparse(formulaList[[i]]), summaryTable(modelList[[i]]))
    
    sumTable <- rbind(sumTable, tempRow)
    
  }
  
  # get the ranking table to select model, then return both
  modelRank <- modelRank(sumTable)
  modelAndTable <- list("Table" = as_tibble(modelRank$table, .name_repair = "minimal"), "model" = modelList[[modelRank$selectedModel]])
  return(modelAndTable)
  
}


# Addapted from Ditiatkovski M, 2020
# Function makes and arranges Diagnostic plots for Deterministic type variables.
# Expects input of data and the fitted model. If BoxCox transformation has been applied
# supplying lambda will reverse the transformation
plotRes <- function(data, model, lambda = 1, vertical = TRUE) {
  
  # If needed backtransform data and fitted values
  if(lambda!=1){
    fittedVals <- data.frame("data" = InvBoxCox(fitted(model), lambda = lambda))
    data <- InvBoxCox(data, lambda = lambda)
  }else{
    fittedVals <- data.frame("data" = fitted(model))
  }
  
  residuals <- data.frame("residuals" = rstudent(model))
  
  
  time <- as.vector(time(data))
  
  # Original vs fitted values
  dataVsModel <- ggplot(data = data, aes(x = time, y = data)) + 
    geom_line()  +
    geom_line(data = fittedVals, aes(y = data), colour = "mediumpurple") +
    geom_point(size = 0.4)+
    xlab("Years") +
    ylab("Number of NLRD applications") +
    theme_bw()
  
  # residuals ploted
  resid <- ggplot(data = residuals, aes(x = time, y = residuals))+
    geom_line() +
    geom_point(size = 0.2) +
    theme_bw() +
    xlab("Years") +
    ylab("Standardised residuals")
  
  # Fitted vs residuals
  fitRes <- data.frame("residuals" = unlist(residuals), "fitted" = unlist(fittedVals))
  
  fitVsRes <- ggplot(data = fitRes, aes(x = fitted, y = residuals)) + 
    geom_point() + 
    xlab("Fitted model values") +
    ylab("Standardised residuals") +
    theme_bw()
  
  # create ACF plot
  ACF <- ggAcf(residuals, lag.max = 48)+
    theme_bw() +
    theme(plot.title = element_blank())
  
  # create histogram with normal overlay (this is more for aesthetic reasons in the figure)
  hist <- forecast::gghistogram(residuals[[1]], add.normal = T) +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          plot.title = element_blank(),)
  
  # create qq plot for the residuals
  qq <- ggplot(data = residuals, aes(sample = residuals)) +
    stat_qq() +
    stat_qq_line(linetype = "dashed", color = "mediumpurple") +
    theme_bw()
  
  
  
  ## arrange everything nicely and add labels
  if (vertical){
    plotsComb = ggarrange(dataVsModel, fitVsRes, resid, ACF, 
                          hist, qq, labels = c("A","B","C","D","E","F"), ncol = 2, nrow = 3)
  }else{
    
    plotsComb = ggarrange(dataVsModel, fitVsRes, resid, ACF, 
                          hist, qq, labels = c("A","B","C","D","E","F"), ncol = 3, nrow = 2)
  }
  
  
  return(plotsComb)
  
}

################################
###Determinisitc model codes####
################################

# define significant figure numbers
options(pillar.sigfig = 5)

gmo <- read_xlsx("NLRD list December 2019.xlsx")

colnames(gmo)[5] <- "dates"

gmo1 <- gmo %>% group_by(year(dates), month(dates))
gmo2 <- gmo1 %>% summarise(n())
gmo2 <- gmo2 %>% arrange(.,`year(dates)`,`month(dates)`)
gmo.ts <- ts(gmo2$`n()`, frequency = 12, start = c(2001,7))

#Plot Figure 1
plot(gmo.ts, type = 'l', main = "GMO monthly research application", ylab = "Number of applications")

#Plot Figure 2 
plot(stl(gmo.ts, s.window = "periodic"),cex.axis=20,main="STL decomposition of the GMO Time Series  ",cex.main=10)

#reduce seasons
gmo.ts.short <- ts(gmo.ts[25:length(gmo.ts)], start=c(2003,07), frequency = 12)

#Plot Figure 3
plot(gmo.ts.short, type = 'l', main = "GMO monthly research application", ylab = "Number of applications")
paste("p-value for shapirowilk test of normality for the short series is ",
      shapiro.test(gmo.ts.short)$p.value)

#Plot figure 4
#Box Cox tranform the data
bc <- BoxCox.ar(gmo.ts.short, method = 'mle')

#set up for lambda
lambda = bc$mle  #lambda = 0.1
gmo.bc = (gmo.ts.short^lambda-1)/lambda

paste("p-value for shapirowilk test of normality for the short series is ",
      shapiro.test(gmo.bc)$p.value)

#Plot Figure 5 
plot(gmo.bc, type = 'l', main = "GMO monthly research application", ylab = "Transformed values")

#Output Table 1
shortGmoSeasons <- season(gmo.ts.short)
shortSeasonalDeterministic <-  modelSelect(gmo.bc, 3, shortGmoSeasons)
shortSeasonalDeterministic$Table

# Output Table 2
detSum <- summary(shortSeasonalDeterministic$model)

detSum$coefficients

write.table(detSum$coefficients, "Deterministic Coefs.csv", sep = ",", row.names = FALSE)



#Plot Figure 6

#png("DeterministicHoriz.png", width = 1280, height = 720)
plotRes(gmo.bc, shortSeasonalDeterministic$model, 0.1, F)
#dev.off()


################################################
####Self defined function for Stochastic part###
################################################

acfpacf <- function(Data){
  par(mfrow=c(1,2))
  acf(Data,  ci.type= "ma", main = 'ACF test', lag.max = 60)
  pacf(Data, main = 'PACF test', lag.max = 60)
  par(mfrow=c(1,1))
}


# addapted from Ditiatkovski M. time series assignment 2
# assumes time series as an input, returns ACF PACF and data plots
plotDiag<- function(data, ylab = ""){
  
  x <- as.vector(time(data))
  
  diffData <- ggplot(data = data, aes(x = x, y = as.vector(data))) + 
    geom_line() +
    geom_point() +
    xlab("Years") +
    ylab(ylab) +
    theme_bw()
  
  # create ACF plot
  ACF <- ggAcf(data, lag.max = 48) +
    theme_bw() +
    theme(plot.title = element_blank())
  
  # create PACF plot
  PACF <- ggPacf(data, lag.max = 48) +
    theme_bw() +
    theme(plot.title = element_blank())
  
  # return diagnostic plot
  row2 <- ggarrange(ACF, PACF, labels = c("B","C"))
  return(ggarrange(diffData, row2, labels = "A", nrow = 2))
}


# This function, given the time series data, ARIMA orders and seasonal orders will fit ARIMA models
# The list of models will then be analysed by another function to find the best one 
modelFit <- function(ts, order, seasonal=c(0,0,0), lambda = 1, methods = c("ML","CSS","CSS-ML"), freq = 12){
  
  modelList <- list()
  
  # Try to fit the model with every single method by iterating though method list
  for (i in methods){
    # Error catch in case arima function fails
    errorCheck <- tryCatch({
      # Attempt to fit Model
      model <- Arima(ts, order = order, lambda = lambda, method = i, seasonal=list(order = seasonal, period = 12))
      # Get k of a model and error summary for fitting
      k = sum(order, seasonal)
      modelAccuracy <- as_tibble(accuracy(model))
      # make an entry to ModelList type of an object
      modelList = rbind(modelList, list("order" = order,
                                        "seasonal" = seasonal,
                                        "model" = model,
                                        "accuracy" = modelAccuracy,
                                        "k" = k,
                                        "method" = i))
    }, error=function(er){  
      # if arima function fails show a warning, continue on to the next model
      warning(paste0(c("ARIMA(", paste0(seasonal, collapse=","),")X(", paste0(order, collapse=","), ") with ", i, " method could not be fitted")))})
  }
  return(modelSelectError(modelList))
}

# expects a list models and a accuracy (from model_fit()) and finds the model with the 
# smallest error, default being MAE. Returns a list with that model specifications 
modelSelectError <- function(modelList, error = "MAE"){
  combAccuracy <- bind_rows(modelList[,"accuracy"])
  lowestError <- which(combAccuracy[error] == min(combAccuracy[error]))
  return(modelList[lowestError,])
}

# Calculate BIC for Arima type model
BIC <- function(model, k){
  return(-2*model$loglik + k*log(length(model$residuals)) %>% round(3))
}

# Calculate AIC for Arima type model
AIC <- function(model, k){
  return(-2*model$loglik + 2*k %>% round(3)) 
}


# Model to fit models given a list of lists, the latter of which contain a single order named "seasonal" for seasonal
# data and a list named "order" for ARIMA orders to try for that seasonal model.
# Other inputs are needed to send to model_fit function. Returns a list of models with their diagnostic and a table

modelListFit <- function(ts, orders, lambda = 1, methods = c("ML","CSS","CSS-ML"), freq = 12){
  modelList = list()
  modelCounter = 1
  
  # outer loop to itterate through seasonal part of models
  for(i in seq(1, length(orders))){
    
    
    seasonal = orders[[i]]$seasonal
    orderList = orders[[i]]$order
    
    # Iterate through the normal Arima orders, keeping seasonal orders constant to fit all the models 
    for(j in orderList){
      # for diagnostics
      print(paste0(c("seas: ", seasonal, " Order: ", j), collapse = ""))
      # some models fail spactacularly 
      errorCheck <- tryCatch({
        # Find the best model for a particular order
        temp <- modelFit(ts = ts,
                         order = j, 
                         seasonal = seasonal, 
                         lambda = lambda,
                         methods = methods,
                         freq = freq)
        if(is.null(temp)){
          next()
        }
        # add model, and model index to the modelList object
        modelList <- rbind(modelList,
                           c(temp, "ModelIndex" = modelCounter))
        modelCounter = modelCounter + 1  
      }, error=function(er){  
        # usually shows if all models fail.
        warning(paste0(c("ARIMA(", paste0(seasonal, collapse=","),")X(", paste0(i, collapse=","), ") could not be fitted")))})
      
    }
  } 
  # the return fucntion finalises the modelList type output to contain all models and a summary table
  return(list("ModelList" = modelList, "Table" = modelSummaryTable(modelList)))  
  
}

# Create a summary table of the supplied models. Expects input generated in model fit, 
# or combination of those in a list
modelSummaryTable <- function(modelList){
  
  # setup data frame / columns for the table
  BIC <- mapply(BIC, modelList[,"model"], modelList[,"k"])
  AIC <- mapply(AIC, modelList[,"model"], modelList[,"k"])
  accuracy <- bind_rows(modelList[,"accuracy"])
  
  # Get all orders
  order <- Reduce(rbind, modelList[,"order"]) %>% as_tibble() %>% 
    unite(., "ARIMA order", sep = ",")
  
  seasonalOrder <- Reduce(rbind, modelList[,"seasonal"]) %>% as_tibble() %>% 
    unite(., "Seasonal order", sep = ",")
  
  # Add information about method and number of parameters 
  method <- Reduce(rbind, modelList[,"method"])[,1]
  k <- Reduce(rbind, modelList[,"k"])[,1]
  
  # Find model normality and proportion of coefficents that are significant
  normality <- sapply(modelList[,"model"], residualNormPval)
  propSigCoefs <- sapply(modelList[,"model"], significantCoefProp)
  
  # Populate model index in the table, so it can easily taken from model list
  index <- Reduce(rbind, modelList[,"ModelIndex"])[,1]
  
  #make the table
  results <- tibble(order, 
                    seasonalOrder,
                    "k" = k,
                    "Method" = method,
                    "MASE" = accuracy$MASE,
                    "RMSE" = accuracy$RMSE,
                    "MAE" = accuracy$MAE,
                    "AIC" = AIC,
                    "BIC" = BIC,
                    "Res. Normality" = normality,
                    "Prop. Sig. coefs" = propSigCoefs,
                    "ModelIndex" = index)
  
  return(results %>% arrange(`MASE`, `BIC`, desc(`Prop. Sig. coefs`)))
  
}

# Calculate p value for shapiro test of an Arima model
residualNormPval <- function(model){
  return(shapiro.test(model$residuals)$p.value)
}

# Calculate proportion of significant coeficents in Arima model
significantCoefProp <- function(model){
  coefs <- coeftest(model) %>% broom::tidy() %>% .$p.value
  return(sum(coefs<0.05, na.rm = T)/length(coefs))
}


# This function is adapted from Ditiatkovski M, 2020, Timeseries Assignment 2
# Custom functions that plots diagnostic plots for model residuals
# expects original ts data and model as an input, vertical controls if the returned
# image is vertical (TRUE) or horizontal (FALSE). Kcalc controls if
# number of parameters are calculated for LBtest.

plotResArima <- function(data, model, title, vertical = T, kcalc = T) {
  
  combTable <- tibble("time" = as.vector(time(data)),
                      "data" = data, 
                      "residuals" = model$residuals,
                      "fitted" = fitted(model)) 
  # get number of coefficients assumes ar fits all parameters
  
  dump <- capture.output(model)
  
  # if needed calculate k
  if(kcalc){
    
    k<- dump[2] %>% str_extract_all(., "[[:digit:]]") %>%
      .[[1]] %>%
      as.numeric() %>%
      sum()
  }else{
    k = 0
  }
  
  lag.max <-if_else(length(combTable$data)>40, 40, length(combTable$data)-1) 
  
  # Original vs fitted values
  fittedPlot <- ggplot(data = combTable, aes(x = time, y = data)) + 
    geom_line() +
    geom_line(aes(y = fitted), colour = "mediumpurple", 
              linetype = 2) +
    geom_point(size = 0.4) +
    ylab(title) +
    theme_bw()
  
  
  # residuals plotted
  residualsPlot <- ggplot(data = combTable, aes(x = time, y = residuals))+
    geom_line() +
    geom_point(size = 0.2) +
    theme_bw() +
    ylab("Standardised residuals")
  
  
  # create ACF plot
  ACF <- ggAcf(combTable$residuals, lag.max = lag.max) +
    theme_bw() +
    theme(plot.title = element_blank())
  
  # create ACF plot
  PACF <- ggPacf(combTable$residuals, lag.max = lag.max) +
    theme_bw() +
    theme(plot.title = element_blank())
  
  # create qq plot for the residuals
  qq <- ggplot(data = combTable, aes(sample = residuals)) +
    stat_qq() +
    stat_qq_line(linetype = "dashed", colour = "mediumpurple") +
    theme_bw()
  
  # create Ljung Box plot in ggplot 
  LBTPval <- FitAR::LjungBoxTest(combTable$residuals, k = 0,
                                 lag.max = lag.max) %>% 
    as.tibble()
  
  LBplot <- ggplot(LBTPval, aes(x = m, y = pvalue)) +
    geom_point() +
    geom_abline(aes(intercept = 0.05, slope = 0), 
                colour = "dodgerblue4", linetype = 2) +
    xlab("Lag") +
    ylab("p-value") +
    theme_bw()
  
  ## arrange everything nicely and add labels
  if (vertical){
    plotsComb = ggarrange(fittedPlot, residualsPlot, ACF, PACF, 
                          qq, LBplot, labels = c("A","B","C","D","E","F"), ncol = 2, nrow = 3)
  }else{
    
    plotsComb = ggarrange(fittedPlot, residualsPlot, qq, ACF, PACF, 
                          LBplot, labels = c("A","B","C","D","E","F"), ncol = 3, nrow = 2)
  }
  
  
  return(plotsComb)
}



# Returns coefficent diagnostics in a table for a supplied model
coefTable <- function(model){
  return(coeftest(model) %>% round(., 4) %>% tidy())
}


# for use in presentation, modified from previous assignment
plotAcfPacf<- function(data, ylab){
  time <- as.vector(time(data))
  
  # Plot differenced data
  diffData <- ggplot(data = data, aes(x = time, y = data)) + 
    geom_line() +
    geom_point(size = 0.2) +
    xlab("Years") +
    ylab(ylab) +
    theme_bw()
  
  # create ACF plot
  ACF <- ggAcf(data, lag.max = 40) +
    theme_bw() +
    theme(plot.title = element_blank())
  
  # create PACF plot
  PACF <- ggPacf(data, lag.max = 40) +
    theme_bw() +
    theme(plot.title = element_blank())
  
  # show the diagnostic plot
  row2 <- ggarrange(ACF, PACF)
  return(ggarrange(diffData, row2, nrow = 2))
}

#define a function to forecast based on the model and output the plot and data
forecastPlot <- function(Data, order, seasonal, lambda, number_of_period){
  model = Arima(Data, order=order, seasonal = seasonal, lambda = lambda)
  pred <- forecast(model, h = number_of_period)
  print(pred)
  return(autoplot(pred))
}

######################################################
# Code for tables and figures for Stochastic Section #
######################################################

#Plot Figure 7
#show the acf pacf against the transformed data
plotDiag(gmo.bc)
#test for normality of the transformed data
shapiro.test(gmo.bc)

#Plot figure 8
plotDiag(Arima(gmo.ts.short, order=c(0,0,0), seasonal= list(order = c(1,0,0), period  = 12), lambda = 0.1)$residuals)

#Plot figure 9
plotDiag(Arima(gmo.ts.short, order=c(0,0,0), seasonal= list(order = c(0,1,0), period  = 12, lambda = 0.1))$residuals)
# applying a single difference appears to introduce an MA component to the seasonality, which may indicate overdifferencing 

#Unit root test
ArimaS100Res <- Arima(gmo.ts.short, order=c(0,0,0), seasonal= list(order = c(1,0,0), period  = 12), lambda = 0.1)$residuals
adf.test(ArimaS100Res, k = ar(ArimaS100Res)$order)

#Plot figure 10
Arima010S100 <- Arima(gmo.ts.short, order=c(0,1,0), seasonal= list(order = c(1,0,0), period  = 12), lambda = 0.1)$residuals
adf.test(Arima010S100, k = ar(Arima010S100)$order)
plotDiag(Arima010S100)

# Figure 11
eacf(Arima010S100)

# Figure 12
plot(armasubsets(Arima010S100, nar = 8, nma = 8)) # (6,1,1), (4,1,5)

#Seaonal sarima S(0,1,1)
ArimaS011Res <- Arima(gmo.ts.short, order=c(0,0,0), seasonal= list(order = c(0,1,1), period  = 12), lambda = 0.1)$residuals
#ADF test
adf.test(ArimaS011Res, k = ar(ArimaS011Res)$order)
#  Series is non stationary despite lack of trend. A difference can be applied

#ADF test on the differenced series
Arima010S011Res <- Arima(gmo.ts.short, order=c(0,1,0), seasonal= list(order = c(0,1,1), period  = 12), lambda = 0.1)$residuals
adf.test(Arima010S011Res, k = ar(Arima010S011Res)$order)
#Figure 13
plotDiag(Arima010S011Res)

#Figure 14
eacf(Arima010S011Res)
# EACF is easier to interpret, possible models : (0,1,1), (0,1,2), (1,1,1), (1,1,2)

#Figure 15
plot(armasubsets(Arima010S011Res, nar = 7, nma = 7))

#Seaonal sarima S(3,1,1)
ArimaS311Res <- Arima(gmo.ts.short, order=c(0,0,0), seasonal= list(order = c(3,1,1), period  = 12), lambda = 0.1)$residuals
#ADF test
adf.test(ArimaS311Res, k = ar(ArimaS311Res)$order)

#ADF test on the differenced series
Arima010S311Res <- Arima(gmo.ts.short, order=c(0,1,0), seasonal= list(order = c(3,1,1), period  = 12), lambda = 0.1)$residuals
adf.test(Arima010S311Res, k = ar(Arima010S311Res)$order)

#Figure 16
plotDiag(Arima010S311Res)

#Figure 17
eacf(Arima010S311Res)

#Figure 18
plot(armasubsets(Arima010S011Res, nar = 5, nma = 5))

########### Model Evaluation #####################

#Model setup

S311Models = list("seasonal"= c(3,1,1),
                  "order" = list(c(0,1,1), c(3,1,0), c(5,1,0),
                                 c(3,1,1), c(5,1,1), c(0,1,2), 
                                 c(1,1,1), c(1,1,2)))


S011Models = list("seasonal"= c(0,1,1),
                  "order" = list (c(0,1,1), c(0,1,2), c(1,1,1), c(1,1,2), 
                                   c(3,1,0), c(3,1,1), c(0,1,3), c(5,1,0), 
                                   c(5,1,1) , c(5,1,3) ) ) 

S100MOdels = list("seasonal"= c(1,0,0),
                  "order" = list(c(0,1,1),	c(0,1,2),	c(1,1,2),
                                 c(3,1,0),	c(3,1,1),	c(2,1,3),	
                                 c(6,1,1),	c(4,1,5)))
# combine the models
allModelsOrders = list(S311Models, S011Models, S100MOdels)                 

# Generate all models
allModelFitted =  modelListFit(gmo.ts.short, allModelsOrders, lambda = 0.1)
#Produce all tables
AllTable <- allModelFitted$Table


# original data
ggsave("originalDataPlotAcfPacf.png", plotAcfPacf(gmo.ts, "Number of Apps."),width = 128, height = 72,	units = "mm",
       device = "png")

# seasonal diff
ggsave("SeasonalDiffPlotAcfPacf.png", plotAcfPacf(diff(gmo.ts, lag = 12), "Number of Apps."), width = 128, height = 72,	units = "mm",
       device = "png")

#Write to csv all models and their residual values into a CV
write.table(allModelFitted$Table, "AllModelsDiagnostics.csv", sep = ",", row.names = FALSE)

#Get top 10 models
modelInd <- allModelFitted$Table$ModelIndex[1:10]
top10models <- allModelFitted$ModelList[modelInd,]


############# use this to save residual diagnostic plots, chainge model index as required
# width and hight are just a suggestion.
modelIndex = 1
png("testHoriz.png", width = 1280, height = 720)
plotResArima(gmo.ts.short, top10models[[modelIndex,"model"]], "Number of Apps.", F)
dev.off()

png("testVertical.png", width = 1280, height = 720)
plotResArima(gmo.ts.short, top10models[[modelIndex,"model"]], "Number of Apps.", vertical = T, kcalc = F)
dev.off()
############


for( i in seq(1, 10)){
  order <-  paste0(top10models[[i, "order"]], collapse = ",")
  seasonal <-  paste0(top10models[[i, "seasonal"]], collapse = ",")
  method <- top10models[[i, "method"]]
  fileString <- paste0("ARIMA(",order,")X(",seasonal,") ", method)
  
  png(paste0(fileString, "Horiz.png"),
      width = 1280, height = 720)
  print(plotResArima(gmo.ts.short, top10models[[i,"model"]], "Number of Apps.", F))
  dev.off()
  
  png(paste0(fileString, "Vert.png"),
      width = 1280, height = 720)
  print(plotResArima(gmo.ts.short, top10models[[i,"model"]], "Number of Apps.", T))
  dev.off()
  
  write.table(coefTable(top10models[[i,"model"]]), paste0(fileString, "Coefs.csv"),sep = ",", row.names = FALSE)
}

#overfitting 

# we fit 3 final models 

# 1. our best one (2,1,3)*(1,0,0)
# 2. Overfit 1 by increasing AR (3,1,3)*(1,0,0)
# 3. Overfit 1 by increasing AR (2,1,4)*(1,0,0)

Finalmodels = list(list("seasonal"= c(1,0,0),
                        "order" = list(    c(2,1,3),    c(3,1,3),    c(2,1,4))))

# Generate all models
FinalModelFitted =  modelListFit(gmo.ts.short, Finalmodels, lambda = 0.1)
#Produce all tables
FinalTable <- FinalModelFitted$Table
FinalTable

#original model is better on MASE, MAE, AIC,BIC and sig coefs 

modelInd_final <- FinalTable$ModelIndex[2:3]

overfitmodels <- FinalModelFitted$ModelList[modelInd_final,]

# printing coefs if needed
for( i in seq(1, 2)){
  order <-  paste0(overfitmodels[[i, "order"]], collapse = ",")
  seasonal <-  paste0(overfitmodels[[i, "seasonal"]], collapse = ",")
  method <- overfitmodels[[i, "method"]]
  fileString <- paste0("ARIMA(",order,")X(",seasonal,") ", method)
  
  png(paste0(fileString, "Horiz.png"),
      width = 1280, height = 720)
  print(plotResArima(gmo.ts.short, overfitmodels[[i,"model"]], "Number of Apps.", F))
  dev.off()
  
  png(paste0(fileString, "Vert.png"),
      width = 1280, height = 720)
  print(plotResArima(gmo.ts.short, overfitmodels[[i,"model"]], "Number of Apps.", T))
  dev.off()
  
  write.table(coefTable(overfitmodels[[i,"model"]]), paste0(fileString, "Coefs.csv"),sep = ",", row.names = FALSE)
}

# we can choose our best model ARIMA(2,1,3)*(1,0,0)_12

#Forecasting with the model

forecastPlot(gmo.ts.short,order = c(2,1,3), seasonal = c(1,0,0), lambda = 0.1, number_of_period = 10)


# Reference
# Ditiatkovski M., (2020), MATH1318 Time Series Analysis Assignment 1, multiple functions, RMIT 2020
# Zhang Y.C.(2020) MATH1318 Time Series Analysis Assignment 2, Timeseries Functions Chris, RMIT 2020

