#install.packages("shiny")
#install.packages("shinydashboard")
#install.packages("rhandsontable")
#install.packages("shinyBS")
#install.packages("rmarkdown")
#install.packages("ggplot2")
library(rhandsontable)
library(shinydashboard)
library(shinyBS)
library(shiny)
library(ggplot2)
library(rmarkdown)

numSamples<-20
ui<-shinyUI(dashboardPage(
    skin="green",
	dashboardHeader(title = "dPCalibRate"),
	dashboardSidebar(
	sidebarMenu(
    	menuItem("Calibrate", tabName = "data", icon = icon("line-chart")),
    	menuItem("Help", tabName = "help", icon = icon("question-circle"),
    		menuSubItem("Citation", tabName = "cite"),
   	 		menuSubItem("FAQ", tabName = "faq"))
    )

),
dashboardBody(
	tabItems(
    tabItem(tabName = "data",
		fluidRow(
			box(
				#HTML("<hr>"),
				title="Manual input data",
				width=6,
				height=150,
				sliderInput("nsample", "Number of samples:", 1, 100, numSamples)
			),
			box(
				title="Upload data",
				width=3,
				height=150,
				actionButton("uploadButton","Upload",icon=icon("upload")),
				bsModal("uploadData","Upload data","uploadButton",
					fileInput('file1', 'Choose file to upload',
                			accept = c(
                  		'text/csv',
                  		'text/comma-separated-values',
                  		'text/tab-separated-values',
                  		'text/plain',
                  		'.csv'
	                		)
      				),
					actionButton("uploadModalButton","Load"),
					uiOutput("uploadSuccess")
				)
			),
			box(
				title="Example data",
				width=3,
				height=150,
    			actionButton("exampleLoad", "Load example dataset")
			)

	   ),
	   fluidRow(
			box(
				title="Data",
				rHandsontableOutput("hot")
			)
		),
	   fluidRow(
			box(width=12,
	        	textInput("reportTitle","Title of the report",value="dPCR calibration report"),
	        	radioButtons('format', 'Document format', c('PDF'), inline = TRUE),
	        	downloadButton("reportGenerate","Generate report")
    	)      
      )
	),
	tabItem(tabName = "cite",
    	box(

		 	h4("Citation"),
		 	p("When using this application to analyse your dPCR data, please cite:"),
		 	p("Vynck, M. et al. (2017). Quality control of digital PCR assays and platforms.")		 	      
      	)
    ),
    tabItem(tabName = "faq",
    	box(

		 	h4("Questions? Feel free to contact me at Matthijs.Vynck@UGent.be.")#,
		 	#p("Feel free to contact me at Matthijs.Vynck@UGent.be.")
 	      
      	)
	)
  )
)
)
)
