library(rhandsontable)
library(shinydashboard)
library(shinyBS)
library(shiny)
library(ggplot2)
library(rmarkdown)
library(lmtest)
library(randtests)
library(sandwich)
library(cowplot)
source("ci.cv.R")
numSamples<-20

server<-function(input, output,session)({

values <- reactiveValues(data=data.frame(Observed = c(rep(0.00,(numSamples))), Expected = c(rep(0.00,(numSamples)))))

adjust.data <- function(df_){
  	n <- nrow(df_)
  	if(n<input$nsample){
		df_[(n+1):input$nsample,1]<-0.0
  		df_[(n+1):input$nsample,2]<-0.0
  		row.names(df_)<-1:nrow(df_)
  	}else{
	  	df_<- df_[1:input$nsample,]
  	}
    df_$Observed<-as.numeric(df_$Observed)
    df_$Expected<-as.numeric(df_$Expected)	
    return(df_)
}

observe({
	if(!is.null(input$hot)){
		values$data <- hot_to_r(input$hot)
	}
})

observeEvent(input$nsample,({
	
	if(!is.null(input$hot)){
		values$data <- hot_to_r(input$hot)
	}
  	values$data<-adjust.data(values$data)
})
)

observeEvent(input$exampleLoad,{
	data.temp<-read.csv("testdata.csv")
	updateSliderInput(session,"nsample",value=nrow(data.temp))
	reorder.cols<-order(colnames(data.temp))   
	data.temp<-data.temp[,reorder.cols]
	data.temp<-data.temp[,c(2,1)]
	values$data<-data.temp
})

observeEvent(input$uploadModalButton, {
	inFile <- input$file1
	if (!is.null(inFile)){
		data.temp<-read.csv(inFile$datapath)
	} else {
    	output$uploadSuccess<-renderText("No data uploaded.")
		return()
	}
	if(ncol(data.temp)==2&&sum(sort(colnames(data.temp))==c("Expected","Observed"))==2){
    		updateSliderInput(session,"nsample",value=nrow(data.temp))   
		reorder.cols<-order(colnames(data.temp))
		data.temp<-data.temp[,reorder.cols]
		data.temp<-data.temp[,c(2,1)]
		values$data<-data.temp
		output$uploadSuccess<-renderText("Load succesful. Close this window to continue.")
    } else {
    	output$uploadSuccess<-renderText("Load failed. Wrong data format.")
    }
})


output$hot <- renderRHandsontable({
        rhandsontable(values$data,readOnly=F, useTypes= TRUE) %>%
      		hot_table(highlightCol = TRUE, highlightRow = FALSE) %>%
  			hot_col("Observed", format = "0.00") %>%
  			hot_col("Expected", format = "0.00")
})   

output$reportGenerate <- downloadHandler(
	filename = function() {
		paste(paste("report",strftime(strptime(date(),format="%a %b %d %H:%M:%S %Y"),format="%Y%m%d%H%M%S"), sep=''), sep = '.', 'pdf')
	},
	content = function(file) {
		src <- normalizePath('report.Rmd')
		df_ <- values$data
		# temporarily switch to the temp dir, in case you do not have write
		# permission to the current working directory
		owd <- setwd(tempdir())
		on.exit(setwd(owd))
		file.copy(src, 'report.Rmd',overwrite=TRUE)
		library(rmarkdown) 
		outReport <- rmarkdown::render('report.Rmd', pdf_document(),params=list(ls()))
		file.rename(outReport, file)
	}	
)

})
