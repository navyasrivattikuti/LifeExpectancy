# LifeExpectancy

## People die everyday but a few geographies and periods in history have more people dying at younger ages than other geographies. What makes people die young or live longer? Our goal is to pick a few attributes and explore via analysis what changes or how do they impacted the life expectancy across geographies/time period.

## The Global life expectancy has seen an increase by 6 years FROM A MERE 66.8 to 73.4 years in the last 2 decades (Roughly from 2000 to 2019). Mostly Health care resources play a major role in affecting Life Expectancy. Wealthier countries have a higher average Life Expectancy than poorer countries. Most of the research showed that Life Expectance is majorly affected by the resources available in that country and depends on how they are utilized.

## The project relies on accuracy of data. The Global Health Observatory (GHO) data repository under World Health Organization (WHO) keeps track of the health status as well as many other related factors for all countries The data-sets are made available to public for the purpose of health data analysis. The data-set related to life expectancy, health factors for 193 countries has been collected from the same WHO data repository website and its corresponding economic data was collected from United Nation website. Among all categories of health-related factors only those critical factors were chosen which are more representative. It has been observed that in the past 15 years , there has been a huge development in health sector resulting in improvement of human mortality rates especially in the developing nations in comparison to the past 30 years. Therefore, in this project we have considered data from year 2000-2015 for 193 countries for further analysis. The individual data files have been merged together into a single data-set. On initial visual inspection of the data showed some missing values. As the data-sets were from WHO, we found no evident errors. Missing data was handled in R software by using Missmap command. The result indicated that most of the missing data was for population, Hepatitis B and GDP. The missing data were from less known countries like Vanuatu, Tonga, Togo, Cabo Verde etc. Finding all data for these countries was difficult and hence, it was decided that we exclude these countries from the final model data-set. The final merged file(final dataset) consists of 22 Columns and 2938 rows which meant 20 predicting variables. All predicting variables was then divided into several broad categories:​Immunization related factors, Mortality factors, Economical factors and Social factors.

## Hypothesis: 
## Alcohol Consumption is one of the biggest factors affecting life expectancy

## Research Question: 
## What changes are needed for a country to improve Life Expectancy

## Null Hypothesis, H0: There is no relationship between Life Expectancy and the independent variables  
## Alternate Hypothesis, Ha: There is a relationship between Life Expectancy and the independent variables

![image](https://user-images.githubusercontent.com/125774977/233864553-a532bcba-748a-4491-84cf-614a88042c15.png)

![image](https://user-images.githubusercontent.com/125774977/233864598-600a9a22-e005-4d21-a600-41e25f950549.png)

![image](https://user-images.githubusercontent.com/125774977/233864625-e8360134-3156-4f6c-92e4-afb254a0c4bd.png)




