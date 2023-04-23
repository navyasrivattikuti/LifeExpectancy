rm(list = ls())
#C:\\Users\\navya\\Documents\\Courses\\642 Statistics\\Working Folder\\Project Presentation
#C:\\Users\\LabUser\\Documents\\NS
setwd("C:\\Users\\LabUser\\Documents\\NS")

# Installing the packages and reading the LE_clean in R as a data frame
library(readxl)
Life_Exp <- read.csv("Life_Expectancy.csv")
LE <- data.frame(Life_Exp)

# Basic functions to understand data
head(LE)
nrow(LE)
ncol(LE)

#Renaming Variables
rownames(LE) <- paste0("row", 1:2938)
colnames(LE) <- c('Country', 'Year', 'Status', 'Life_exp', 'Adult_mortality', 'Infant_deaths', 'Alcohol', 'per_Expenditure', 'Hepatitis_B', 'Measels', 'BMI', 'under_5Deaths', 'Polio', 'Total_expenditure', 'Diphtheria', 'HIV_AIDS', 'GDP', 'Population', 'Thinness1', 'Thinness2', 'Income_COR', 'Schooling')
colnames(LE)

# creating new column continent
library(dplyr)
library(countrycode)

LE <- LE %>% 
  mutate(Continent = countrycode(sourcevar = Country, 
                              origin = "country.name",
                              destination = "continent"))
unique(LE$Continent)
unique(LE$Country[LE$Continent == "Asia"])
unique(LE$Country[LE$Continent == "Europe"])
unique(LE$Country[LE$Continent == "Africa"])
unique(LE$Country[LE$Continent == "Americas"])
unique(LE$Country[LE$Continent == "Oceania"])
head(LE)
LE_T <- LE[2:23] # Created a new data set removing country variables
ncol(LE_T) # Since country is removed, it still has same columns

# Data cleaning

glimpse(LE)

LE <- LE %>% 
  mutate(Status = as.factor(Status))

par(mfrow=c(2,7))
boxplot(LE$Life_exp,
        ylab = "Life Expectancy",
        main = "Life Expectancy",
        col= "#FF6666",
        outcol="#FF6666")
boxplot(LE$Adult_mortality,
        ylab = "Adult Mortality",
        main = "Adult Mortality",
        col= "#FF6666",
        outcol="#FF6666")
boxplot(LE$Alcohol,
        ylab = "Alcohol",
        main = "Alcohol",
        col= "#008080",
        outcol="#008080")
boxplot(LE$Hepatitis_B,
        ylab = "Hepatitis B",
        main = "Hepatitis B",
        col= "#FF6666",
        outcol="#FF6666")
boxplot(LE$BMI,
        ylab = "BMI",
        main = "BMI",
        col= "#008080",
        outcol="#008080")
boxplot(LE$Polio,
        ylab = "Polio",
        main = "Polio",
        col= "#FF6666",
        outcol="#FF6666")
boxplot(LE$Total_expenditure,
        ylab = "Total Expenditure",
        main = "Total Expenditure",
        col= "#FF6666",
        outcol="#FF6666")
boxplot(LE$Diphtheria,
        ylab = "Diphteria",
        main = "Diphteria",
        col= "#FF6666",
        outcol="#FF6666")
boxplot(LE$GDP,
        ylab = "GDP",
        main = "GDP",
        col= "#FF6666",
        outcol="#FF6666")
boxplot(LE$Population,
        ylab = "Population",
        main = "Population",
        col= "#FF6666",
        outcol="#FF6666")
boxplot(LE$Thinness1,
        ylab = "Thinness 1-19 years",
        main = "Thinness1",
        col= "#FF6666",
        outcol="#FF6666")
boxplot(LE$Thinness2,
        ylab = "Thinness 5-9 years",
        main = "Thinness2",
        col= "#FF6666",
        outcol="#FF6666")
boxplot(LE$Income_COR,
        ylab = "Income Composition",
        main = "Income Composition",
        col= "#008080",
        outcol="#008080")
boxplot(LE$Schooling,
        ylab = "Schooling",
        main = "Schooling",
        col= "#FF6666",
        outcol="#FF6666")

# These plots show none or few of the outliers in Alcohol, BMI, Income COR
# Replacing higher outliers with median values and lower outliers with mean values

LE_clean <- LE

#medians
LE_clean$Life_exp[is.na(LE_clean$Life_exp)] <- median(LE_clean$Life_exp, na.rm=TRUE)
LE_clean$Adult_mortality[is.na(LE_clean$Adult_mortality)] <- median(LE_clean$Adult_mortality, na.rm=TRUE)
LE_clean$Hepatitis_B[is.na(LE_clean$Hepatitis_B)] <- median(LE_clean$Hepatitis_B, na.rm=TRUE)
LE_clean$Polio[is.na(LE_clean$Polio)] <- median(LE_clean$Polio, na.rm=TRUE)
LE_clean$Diphtheria[is.na(LE_clean$Diphtheria)] <- median(LE_clean$Diphtheria, na.rm=TRUE)
LE_clean$Total_expenditure[is.na(LE_clean$Total_expenditure)] <- median(LE_clean$Total_expenditure, na.rm=TRUE)
LE_clean$GDP[is.na(LE_clean$GDP)] <- median(LE_clean$GDP, na.rm=TRUE)
LE_clean$Population[is.na(LE_clean$Population)] <- median(LE_clean$Population, na.rm=TRUE)
LE_clean$Thinness1[is.na(LE_clean$Thinness1)] <- median(LE_clean$Thinness1, na.rm=TRUE)
LE_clean$Thinness2[is.na(LE_clean$Thinness2)] <- median(LE_clean$Thinness2, na.rm=TRUE)
LE_clean$Schooling[is.na(LE_clean$Schooling)] <- median(LE_clean$Schooling, na.rm=TRUE)
#means
LE_clean$Alcohol[is.na(LE_clean$Alcohol)] <- mean(LE_clean$Alcohol, na.rm=TRUE)
LE_clean$BMI[is.na(LE_clean$BMI)] <- mean(LE_clean$BMI, na.rm=TRUE)
LE_clean$Income_COR[is.na(LE_clean$Income_COR)] <- mean(LE_clean$Income_COR, na.rm=TRUE)



anyNA(LE_clean)
#colSums(is.na(LE_clean))
summary(LE_clean)

LE_selected <- LE_clean %>% 
         select(-Country, -Year, -Continent)
summary(LE_selected)
glimpse(LE_selected)

library(moments)
skewness(LE_clean[4:20])


library(ggplot2)
library(tidyr)
library(purrr)
library(qqplotr)

df = pivot_longer(LE_clean, cols = everything())

my_qqplot <- function(.data, .title) {
  ggplot(data = .data, mapping = aes(sample = value)) +
    stat_qq_line() +
    stat_qq_point() +
    facet_wrap(~ name, scales = "free") +
    labs(x = "Theoretical Quantiles", y = "Sample Quantiles", title = .title)
}

qqplots <- df %>% 
  split(.$name) %>% 
  imap(my_qqplot)

# summarize the variables country, year and dependent variable 
summary(LE_clean[ , c(1, 2, 4)])
mean(LE_selected$Life_exp)


# Correlation plots
library(corrplot)
#?corrplot.mixed
par(mar = c(10, 10, 10, 10)) 
dev.off() 
corrplot.mixed(cor(LE_clean[ , 4:22], use="pairwise.complete.obs"),
               lower = "pie",
               upper = "number",
               diag = "l",
               t1.pos = "lt")
corrplot(cor(LE_clean[ , 4:22], use="pairwise.complete.obs"),
         t1.pos = 'td',
         method = "number",
         number.cex = 0.8, 
         type = "upper",
         order = "original",
         diag = FALSE)

corrplot(corr = cor(LE_clean[ , 4:22], use="pairwise.complete.obs"),
         method = "circle",
         type = "upper",
         tl.pos = "td",
         order = "original")

# We can observe that Adult_mortality, BMI, Income_COR, HIV/AIDs and Schooling are highly correlated with Life Expectancy. 
# Also we can remove facts under_5deaths and per_expenditure as they are highly correlated with infant deaths & Total_expenditure respectively
# This will ignore multi collinearity

LE_clean <- LE_clean %>% 
  select(-under_5Deaths, -per_Expenditure)
glimpse(LE_clean)

LE_selected <- LE_selected %>% 
  select(-under_5Deaths, -per_Expenditure)
glimpse(LE_selected)


# Understanding the dependent variable vs other factors which are highly correlated
library(ggplot2)
library(ggExtra)
library(gcookbook)
library(ggmosaic)
library(plotly)

p1 <- ggplot(LE_selected[which(LE_clean$Alcohol>0),], aes(x = Alcohol, y = Life_exp, color = factor(Status)))+
  geom_point() + theme_bw() + theme(legend.position = "bottom") + 
  labs(title = "Alcohol consumption vs Life Expectancy")
ggMarginal(p1, type = "histogram", groupColor = TRUE, groupFill = TRUE)

p2 <- ggplot(LE_selected[which(LE_clean$BMI>0),], aes(x = BMI, y = Life_exp, color = factor(Status)))+
  geom_point() + theme_bw() + theme(legend.position = "bottom") + 
  labs(title = "Body Mass Index vs Life Expectancy")
ggMarginal(p2, type = "histogram", groupColor = TRUE, groupFill = TRUE)

p3 <- ggplot(LE_selected[which(LE_clean$Adult_mortality>0),], aes(x = Adult_mortality, y = Life_exp, color = factor(Status)))+
  geom_point() + theme_bw() + theme(legend.position = "bottom") + 
  labs(title = "Adult Mortality vs Life Expectancy")
ggMarginal(p3, type = "histogram", groupColor = TRUE, groupFill = TRUE)

p4 <- ggplot(LE_selected[which(LE_clean$HIV_AIDS>0),], aes(x = log(HIV_AIDS), y = Life_exp, color = factor(Status)))+
  geom_point() + theme_bw() + theme(legend.position = "bottom") + 
  labs(title = "HIV/AIDS vs Life Expectancy")
ggMarginal(p4, type = "histogram", groupColor = TRUE, groupFill = TRUE)

p5 <- ggplot(LE_selected[which(LE_clean$Income_COR>0),], aes(x = Income_COR, y = Life_exp, color = factor(Status)))+
  geom_point() + theme_bw() + theme(legend.position = "bottom") + 
  labs(title = "Income Composition of Resources vs Life Expectancy")
ggMarginal(p1, type = "histogram", groupColor = TRUE, groupFill = TRUE)

p6 <- ggplot(LE_selected[which(LE_clean$Schooling>0),], aes(x = Schooling, y = Life_exp, color = factor(Status)))+
  geom_point() + theme_bw() + theme(legend.position = "bottom") + 
  labs(title = "Schooling vs Life Expectancy")
ggMarginal(p1, type = "histogram", groupColor = TRUE, groupFill = TRUE)



## Plot the data to see whether a linear regression function model the data adequately
library(lmtest)
library(sandwich)

# fit a simple linear model for cleaned data

ln_mod <- lm(Life_exp ~ Income_COR + Alcohol, data = LE_clean)
coeftest(ln_mod, vcov. = vcovHC, type = "HC1")
summary(ln_mod)$adj.r.squared


coeftest(ln_mod, vcov. = vcovHC, type = "HC1")
par(mfrow = c(2,2), mar = c(1, 6, 6, 1))
plot(ln_mod)
dev.off()
# plot the observations
plot(LE_clean$Income_COR, LE_clean$Life_exp,
     col = "steelblue",
     pch = 20,
     xlab = "Income_COR", 
     ylab = "Life Expectancy",
     cex.main = 0.9,
     main = "Life Expectancy vs Income_COR and a Linear OLS Regression Function")


# add the regression line to the plot
abline(ln_mod, 
       col = "red", 
       lwd = 2)

# fit a simple linear model for selected data

sn_mod <- lm(Life_exp ~ Income_COR + Alcohol, data = LE_selected)
coeftest(sn_mod, vcov. = vcovHC, type = "HC1")
summary(sn_mod)$adj.r.squared

par(mfrow = c(2,2), mar = c(1, 6, 6, 1))
plot(sn_mod)

# plot the observations
plot(LE_selected$Income_COR, LE_selected$Life_exp,
     col = "steelblue",
     pch = 20,
     xlab = "Income_COR", 
     ylab = "Life Expectancy",
     cex.main = 0.9,
     main = "Life Expectancy vs Alcohol Consumption and a Linear OLS Regression Function")

# add the regression line to the plot
abline(sn_mod, 
       col = "red", 
       lwd = 2)

# For all factors

L_m1 <- lm(Life_exp ~ Income_COR, data = LE_selected)


L_m2 <- lm(formula = Life_exp ~ Income_COR + Alcohol, data = LE_selected)


L_m3 <- lm(formula = Life_exp ~ Income_COR + Alcohol + Schooling , data = LE_selected)


L_m4 <- lm(formula = Life_exp ~ Income_COR + Alcohol + Schooling + 
             Adult_mortality, data = LE_selected)


L_m5 <- lm(formula = Life_exp ~ Income_COR + Alcohol + Schooling +  
             Adult_mortality + Infant_deaths, data = LE_selected)

L_m6 <- lm(formula = Life_exp ~ Income_COR + Alcohol + Schooling +  
             Adult_mortality + Infant_deaths + BMI, data = LE_selected)

L_m7 <- lm(formula = Life_exp ~ Income_COR + Alcohol + Schooling +  
             Adult_mortality + Infant_deaths + BMI + Total_expenditure + GDP, 
             data = LE_selected)

L_m8 <- lm(formula = Life_exp ~ Income_COR + Alcohol + Schooling +  
             Adult_mortality + Infant_deaths + BMI + Total_expenditure + GDP + 
             Thinness1 + Thinness2, data = LE_selected)

L_m9 <- lm(formula = Life_exp ~ Income_COR + Alcohol + Schooling +  
             Adult_mortality + BMI + Total_expenditure + GDP + 
             Thinness1 + Hepatitis_B + Measels + Polio + Diphtheria +
             HIV_AIDS, data = LE_selected)
summary(L_m9)
coeftest(L_m9, vcov. = vcovHC, type = "HC1")
par(mfrow = c(2,2), mar = c(1, 6, 6, 1))
plot(L_m9)
library(stargazer)

rob_se <- list(sqrt(diag(vcovHC(L_m1, type = "HC1"))),
               sqrt(diag(vcovHC(L_m2, type = "HC1"))),
               sqrt(diag(vcovHC(L_m3, type = "HC1"))),
               sqrt(diag(vcovHC(L_m4, type = "HC1"))),
               sqrt(diag(vcovHC(L_m5, type = "HC1"))),
               sqrt(diag(vcovHC(L_m6, type = "HC1"))),
               sqrt(diag(vcovHC(L_m7, type = "HC1"))),
               sqrt(diag(vcovHC(L_m8, type = "HC1"))),
               sqrt(diag(vcovHC(L_m9, type = "HC1"))))


stargazer(L_m1, L_m2, L_m3, 
          L_m4, L_m5, L_m6, L_m7, L_m8, L_m9,
          type="html",
          digits = 3,
          se = rob_se,
          title = "Linear Panel Regression Models of Life Expectancy",
          out="Linear_models.htm")

# Log model
K_m1 <- lm(formula = Life_exp ~ Income_COR + Alcohol + Schooling +  
             Adult_mortality + BMI + Total_expenditure + GDP + 
             Thinness1 + Hepatitis_B + Measels + Polio + Diphtheria +
             HIV_AIDS, data = LE_selected)

K_m2 <- lm(log(Life_exp) ~ Income_COR + Alcohol + Schooling +  
             Adult_mortality +  BMI + GDP + 
             Thinness1 +  Hepatitis_B + Measels + Polio + Diphtheria +
             HIV_AIDS, data = LE_selected)

K_m3 <- lm(Life_exp ~ Alcohol + Alcohol^2 + Income_COR + 
             Adult_mortality + BMI + Total_expenditure + GDP + 
              Thinness2 + Hepatitis_B + Measels + Polio + Diphtheria +
             HIV_AIDS, data = LE_selected)

K_m4 <- lm(formula = Life_exp ~ Income_COR + Income_COR^2 + Alcohol + Schooling +  
             Adult_mortality + BMI + Total_expenditure + GDP + 
             Thinness1 + Hepatitis_B + Measels + Polio + Diphtheria +
             HIV_AIDS, data = LE_selected)

K_m5 <- lm(formula = Life_exp ~ Income_COR + Income_COR^2 + Income_COR^3  + Alcohol + Schooling +  
             Adult_mortality + BMI + Total_expenditure + GDP + 
             Thinness1 + Hepatitis_B + Measels + Polio + Diphtheria +
             HIV_AIDS, data = LE_selected)

par(mfrow = c(2,2), mar = c(1, 6, 6, 1))
plot(K_m2)
par(mfrow = c(2,2), mar = c(1, 6, 6, 1))
plot(K_m3)
par(mfrow = c(2,2), mar = c(1, 6, 6, 1))
plot(K_m4)
par(mfrow = c(2,2), mar = c(1, 6, 6, 1))
plot(K_m5)

rob_se <- list(sqrt(diag(vcovHC(K_m1, type = "HC1"))),
               sqrt(diag(vcovHC(K_m2, type = "HC1"))),
               sqrt(diag(vcovHC(K_m3, type = "HC1"))),
               sqrt(diag(vcovHC(K_m4, type = "HC1"))),
               sqrt(diag(vcovHC(K_m5, type = "HC1"))))


stargazer(K_m1, K_m2, K_m3, K_m4, K_m5,
          type="html",
          digits = 3,
          se = rob_se,
          title = "Different Models of Life Expectancy",
          out="K_models.htm")

# subset the cleaned data
C2001 <- subset(LE_clean, Year == "2001")
#summary(C2000)

C2014 <- subset(LE_clean, Year == "2014")
#summary(C2014)
#summary(C2014$Alcohol)

# estimate simple regression models using 2000 and 2014 data
C2001_mod <- lm(Life_exp ~ Income_COR + Alcohol, data = C2001)
C2014_mod <- lm(Life_exp ~ Income_COR + Alcohol, data = C2014)

coeftest(C2001_mod, vcov. = vcovHC, type = "HC1")
coeftest(C2014_mod, vcov. = vcovHC, type = "HC1")
summary(C2001_mod)$adj.r.squared
summary(C2014_mod)$adj.r.squared

dev.off()

# plot the observations for LE2000
plot(C2001$Income_COR, C2001$Life_exp,
     col = "steelblue",
     pch = 20,
     xlab = "Income_COR", 
     ylab = "Life Expectancy",
     cex.main = 0.9,
     main = "Life Expectancy vs Income_COR in the year 2000")

# add the regression line to the plot LE2000
abline(C2001_mod, 
       col = "red", 
       lwd = 2)


# plot the observations for LE2014
plot(C2014$Income_COR, C2014$Life_exp,
     col = "steelblue",
     pch = 20,
     xlab = "Income_COR", 
     ylab = "Life Expectancy",
     cex.main = 0.9,
     main = "Life Expectancy vs Income_COR in the year 2014")

# add the regression line to the plot LE2014
abline(C2014_mod, 
       col = "red", 
       lwd = 2)

# OLS is not a good model to use for this because it clearly shows a positive coefficient
# which implies that increase in alcohol consumption increases life expectancy


# Panel data declaration
library(plm)

# PLM for cleaned data

C_mod <- plm(Life_exp ~ Income_COR + Alcohol, 
              data = LE_clean,
              index = c("Country", "Year"), 
              model = "within")

# obtain a summary based on clusterd standard errors 
coeftest(C_mod, vcov = vcovHC, type = "HC1")
summary(C_mod)$adj.r.squared

# estimate a combined time and entity fixed effects regression model
C2_mod <- plm(Life_exp ~ Income_COR + Alcohol, 
               data = LE_clean,
               index = c("Country", "Year"), 
               model = "within", 
               effect = "twoways")    ## individual + time

# obtain a summary based on clusterd standard errors 
coeftest(C2_mod, vcov = vcovHC, type = "HC1")
summary(C2_mod)$adj.r.squared

# including Income_COR
C3_mod <- plm(Life_exp ~ Income_COR + Alcohol, 
              data = LE_clean,
              index = c("Country", "Year"), 
              model = "within", 
              effect = "twoways")    ## individual + time

# summary based on clusterd standard errors 
coeftest(C3_mod, vcov = vcovHC, type = "HC1")

# Panel Data with Two Time Periods: "Before and After" Comparisons
# compute the differences 
diff_Life_exp <- C2014$Life_exp - C2001$Life_exp
diff_Alcohol <- C2014$Alcohol - C2001$Alcohol

# estimate a regression using differenced data
LE_diff_mod <- plm(diff_Life_exp ~ diff_Alcohol, index = "Country", model = "fd", data = LE_clean)
coeftest(LE_diff_mod, vcov = vcovHC, type = "HC1")

# Alcohol has negative coefficient but p>0.05 - it is not significant at 5% significance level 

# plot the differenced data
plot(diff_Alcohol, diff_Life_exp,
     col = "steelblue",
     pch = 20,
     xlab = "Difference in Alcohol Consumption", 
     ylab = "Difference in Life Expectancy",
     cex.main = 0.9,
     main = "Differences of Life Expectancy vs Alcohol Consumption in years 2000 & 2014")

# add the regression line to plot
abline(LE_diff_mod, 
       col = "red", 
       lwd = 2)

# Estimate all ten models using plm() for cleaned data
C_m1 <- lm(Life_exp ~ Income_COR + Alcohol, 
           data = LE_clean)

C_m2 <- plm(Life_exp ~ Income_COR + Alcohol + Country, 
            data = LE_clean)

C_m3 <- plm(Life_exp ~ Income_COR + Alcohol + Status + Continent,
            index = c("Country","Year"),
            model = "within",
            data = LE_clean)

C_m4 <- plm(Life_exp ~ Income_COR + Alcohol + Schooling,
            index = c("Country","Year"),
            model = "within",
            data = LE_clean)

C_m5 <- plm(Life_exp ~ Income_COR + Alcohol + Schooling + Adult_mortality,
            index = c("Country","Year"),
            model = "within",
            data = LE_clean)

C_m6 <- plm(Life_exp ~ Income_COR + Alcohol + Schooling + Adult_mortality + 
              Infant_deaths,
            index = c("Country","Year"),
            model = "within",
            data = LE_clean)

C_m7 <- plm(Life_exp ~ Income_COR + Alcohol + Schooling + Adult_mortality + 
              Infant_deaths + BMI,
            index = c("Country","Year"),
            model = "within",
            data = LE_clean)

C_m8 <- plm(Life_exp ~ Income_COR + Alcohol + Schooling + Adult_mortality + 
              Infant_deaths + Total_expenditure + GDP + Thinness1,
            index = c("Country","Year"),
            model = "within",
            data = LE_clean)

C_m9 <- plm(Life_exp ~ Income_COR + Alcohol + Schooling + Adult_mortality + 
              Infant_deaths + GDP +  Thinness2 + Hepatitis_B + Measels + Thinness1,
            index = c("Country","Year"),
            model = "within",
            data = LE_clean)

C_m10 <- plm(Life_exp ~ Income_COR + Alcohol + Schooling + Adult_mortality + 
              Infant_deaths + Hepatitis_B + Polio + Diphtheria 
              + HIV_AIDS,
            index = c("Country","Year"),
            model = "within",
            data = LE_clean)

coeftest(C_m10, vcov = vcovHC, type = "HC1")

par(mfrow = c(2,2), mar = c(1, 6, 6, 1))
plot(C_m10)

library(stargazer)

# gather clustered standard errors in a list
rob_se <- list(sqrt(diag(vcovHC(C_m1, type = "HC1"))),
               sqrt(diag(vcovHC(C_m2, type = "HC1"))),
               sqrt(diag(vcovHC(C_m3, type = "HC1"))),
               sqrt(diag(vcovHC(C_m4, type = "HC1"))),
               sqrt(diag(vcovHC(C_m5, type = "HC1"))),
               sqrt(diag(vcovHC(C_m6, type = "HC1"))),
               sqrt(diag(vcovHC(C_m7, type = "HC1"))),
               sqrt(diag(vcovHC(C_m8, type = "HC1"))),
               sqrt(diag(vcovHC(C_m9, type = "HC1"))),
               sqrt(diag(vcovHC(C_m10, type = "HC1"))))

# generate the table (html format)
# Regression models: in html format
stargazer(C_m1, C_m2, C_m3, C_m4, C_m5,
          C_m6, C_m7, C_m8, C_m9, C_m10,
          type="html",
          digits = 3,
          se = rob_se,
          title = "Linear Panel Regression Models of Life Expectancy",
          out="C_Panel_models.htm")

# To generate R^2 and adj. R^2 for panel linear models
R_m1 <- lm(Life_exp ~ Income_COR + Alcohol, 
           data = LE_clean)

R_m2 <- plm(Life_exp ~ Income_COR + Alcohol + Country, 
            model = "pool", 
            data = LE_clean)

R_m3 <- plm(Life_exp ~ Income_COR + Alcohol + Status + Continent,
            index = c("Country","Year"),
            model = "pool",
            data = LE_clean)

R_m4 <- plm(Life_exp ~ Income_COR + Alcohol + Schooling,
            index = c("Country","Year"),
            model = "pool",
            data = LE_clean)

R_m5 <- plm(Life_exp ~ Income_COR + Alcohol + Schooling + Adult_mortality,
            index = c("Country","Year"),
            model = "pool",
            data = LE_clean)

R_m6 <- plm(Life_exp ~ Income_COR + Alcohol + Schooling + Adult_mortality + 
              Infant_deaths,
            index = c("Country","Year"),
            model = "pool",
            data = LE_clean)

R_m7 <- plm(Life_exp ~ Income_COR + Alcohol + Schooling + Adult_mortality + 
              Infant_deaths + BMI,
            index = c("Country","Year"),
            model = "pool",
            data = LE_clean)

R_m8 <- plm(Life_exp ~ Income_COR + Alcohol + Schooling + Adult_mortality + 
              Infant_deaths + BMI + Total_expenditure + GDP,
            index = c("Country","Year"),
            model = "pool",
            data = LE_clean)

R_m9 <- plm(Life_exp ~ Income_COR + Alcohol + Schooling + Adult_mortality + 
              Infant_deaths + BMI + Total_expenditure + GDP + Thinness1 + Thinness2
            + Hepatitis_B,
            index = c("Country","Year"),
            model = "pool",
            data = LE_clean)

R_m10 <- plm(Life_exp ~ Income_COR + Alcohol + Schooling + Adult_mortality + 
               Infant_deaths + BMI + Hepatitis_B + Measels + Polio + Diphtheria + HIV_AIDS,
             index = c("Country","Year"),
             model = "pool",
             data = LE_clean)

stargazer(R_m1, R_m2, R_m3, R_m4, R_m5,
          R_m6, R_m7, R_m8, R_m9, R_m10,
          type="html",
          digits = 3,
          se = rob_se,
          title = "Linear Panel Regression Models of Life Expectancy",
          out="R_Panel_RSE.htm")


# Estimate all ten models using plm() for continent instead of country
M_m1 <- lm(Life_exp ~ Income_COR + Alcohol, 
           data = LE_clean)

M_m2 <- plm(Life_exp ~ Income_COR + Alcohol + Continent, 
            data = LE_clean)

M_m3 <- plm(Life_exp ~ Income_COR + Alcohol + Status + Continent,
            index = c("Continent","Year"),
            model = "within",
            data = LE_clean)

M_m4 <- plm(Life_exp ~ Income_COR + Alcohol + Schooling,
            index = c("Continent","Year"),
            model = "within",
            data = LE_clean)

M_m5 <- plm(Life_exp ~ Income_COR + Alcohol + Schooling + Adult_mortality,
            index = c("Continent","Year"),
            model = "within",
            data = LE_clean)

M_m6 <- plm(Life_exp ~ Income_COR + Alcohol + Schooling + Adult_mortality + 
              Infant_deaths,
            index = c("Continent","Year"),
            model = "within",
            data = LE_clean)

M_m7 <- plm(Life_exp ~ Income_COR + Alcohol + Schooling + Adult_mortality + 
              Infant_deaths + BMI,
            index = c("Continent","Year"),
            model = "within",
            data = LE_clean)

M_m8 <- plm(Life_exp ~ Income_COR + Alcohol + Schooling + Adult_mortality + 
              Infant_deaths + Total_expenditure + GDP + Thinness1,
            index = c("Continent","Year"),
            model = "within",
            data = LE_clean)

M_m9 <- plm(Life_exp ~ Income_COR + Alcohol + Schooling + Adult_mortality + 
              Infant_deaths + GDP +  Thinness2 + Hepatitis_B + Measels + Thinness1,
            index = c("Continent","Year"),
            model = "within",
            data = LE_clean)

M_m10 <- plm(Life_exp ~ Income_COR + Alcohol + Schooling + Adult_mortality + 
               Infant_deaths + Hepatitis_B + Polio + Diphtheria 
             + HIV_AIDS,
             index = c("Continent","Year"),
             model = "within",
             data = LE_clean)

# gather clustered standard errors in a list
rob_se <- list(sqrt(diag(vcovHC(M_m1, type = "HC1"))),
               sqrt(diag(vcovHC(M_m2, type = "HC1"))),
               sqrt(diag(vcovHC(M_m3, type = "HC1"))),
               sqrt(diag(vcovHC(M_m4, type = "HC1"))),
               sqrt(diag(vcovHC(M_m5, type = "HC1"))),
               sqrt(diag(vcovHC(M_m6, type = "HC1"))),
               sqrt(diag(vcovHC(M_m7, type = "HC1"))),
               sqrt(diag(vcovHC(M_m8, type = "HC1"))),
               sqrt(diag(vcovHC(M_m9, type = "HC1"))),
               sqrt(diag(vcovHC(M_m10, type = "HC1"))))

# generate the table (html format)
# Regression models: in html format
stargazer(M_m1, M_m2, M_m3, M_m4, M_m5,
          M_m6, M_m7, M_m8, M_m9, M_m10,
          type="html",
          digits = 3,
          se = rob_se,
          title = "Linear Panel Regression Models of Life Expectancy",
          out="M_Panel_models.htm")

# To generate R^2 and adj. R^2 for panel linear models
N_m1 <- lm(Life_exp ~ Income_COR + Alcohol, 
           data = LE_clean)

N_m2 <- plm(Life_exp ~ Income_COR + Alcohol + Continent, 
            model = "pool", 
            data = LE_clean)

N_m3 <- plm(Life_exp ~ Income_COR + Alcohol + Status + Continent,
            index = c("Continent","Year"),
            model = "pool",
            data = LE_clean)

N_m4 <- plm(Life_exp ~ Income_COR + Alcohol + Schooling,
            index = c("Continent","Year"),
            model = "pool",
            data = LE_clean)

N_m5 <- plm(Life_exp ~ Income_COR + Alcohol + Schooling + Adult_mortality,
            index = c("Continent","Year"),
            model = "pool",
            data = LE_clean)

N_m6 <- plm(Life_exp ~ Income_COR + Alcohol + Schooling + Adult_mortality + 
              Infant_deaths,
            index = c("Continent","Year"),
            model = "pool",
            data = LE_clean)

N_m7 <- plm(Life_exp ~ Income_COR + Alcohol + Schooling + Adult_mortality + 
              Infant_deaths + BMI,
            index = c("Continent","Year"),
            model = "pool",
            data = LE_clean)

N_m8 <- plm(Life_exp ~ Income_COR + Alcohol + Schooling + Adult_mortality + 
              Infant_deaths + BMI + Total_expenditure + GDP,
            index = c("Continent","Year"),
            model = "pool",
            data = LE_clean)

N_m9 <- plm(Life_exp ~ Income_COR + Alcohol + Schooling + Adult_mortality + 
              Infant_deaths + BMI + Total_expenditure + GDP + Thinness1 + Thinness2
            + Hepatitis_B,
            index = c("Continent","Year"),
            model = "pool",
            data = LE_clean)

N_m10 <- plm(Life_exp ~ Income_COR + Alcohol + Schooling + Adult_mortality + 
               Infant_deaths + BMI + Hepatitis_B + Measels + Polio + Diphtheria + HIV_AIDS,
             index = c("Continent","Year"),
             model = "pool",
             data = LE_clean)

stargazer(N_m1, N_m2, N_m3, N_m4, N_m5,
          N_m6, N_m7, N_m8, N_m9, N_m10,
          type="html",
          digits = 3,
          se = rob_se,
          title = "Linear Panel Regression Models of Life Expectancy (continent)",
          out="N_Panel_RSE.htm")

#F-Test to see if panel model is significant than linear model

pFtest(C_m10, L_m9)

# P value is less than , so we can say that panel fixed model fit is more significant than the linear model fit
