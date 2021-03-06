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
    skin="red",
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
		 	p("Vynck, M. et al. (2017). Quality control of digital PCR assays and platforms. Analytical and Bioanalytical Chemistry."),
		 	p("Full text available at http://rdcu.be/uU5W."),
		 	p("More information on digital PCR data analysis on http://dpcr.ugent.be.")		 	      
      	)
    ),
    tabItem(tabName = "faq",
    	box(
			h4("Where can I find the underlying methodology?"),
			p("The full text of our paper is available at http://rdcu.be/uU5W."),
 	      	h4("I have problems uploading my data!"),
		 	p("The uploaded data must be in .csv format and contain two columns. The column headers should be named Observed and Expected. An example is available at http://users.ugent.be/~mvynck/calibration.csv.")	,
		 	h4("Questions? Feel free to contact me at Matthijs.Vynck@UGent.be.")
      	)
	)
  )
)
)
)
