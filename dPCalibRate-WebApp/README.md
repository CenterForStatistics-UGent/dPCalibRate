## dPCalibRate - Web Application
This folder contains:
1. server-side (server.R),
2. user interface (ui.R) and
3. R markdown report (report.Rmd)

files in use at http://antonov.ugent.be:3838/dPCalibRate.

These files may also be used to run the application locally from R or RStudio:
```{r eval=FALSE}
runApp("~/AdjustThisPath/dPCalibRate-WebApp")
```

or, alternatively, from your computer's console:
```{r eval=FALSE}
R -e "shiny::runApp("~/AdjustThisPath/dPCalibRate-WebApp")"
```

Make sure you have all of the following libraries installed:
1. shiny
2. shinyBS
3. shinydashboard
4. ggplot2
5. rmarkdown
6. lmtest
7. randtests
8. sandwich
9. cowplot