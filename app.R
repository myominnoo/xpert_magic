

library(shiny)
library(bslib)
library(DT)
library(tidyverse)

source("rscript.R")

options(shiny.maxRequestSize = 30*1024^2)




# ui ----------------------------------------------------------------------


# sidebar -----------------------------------------------------------------

sidebar <-
	sidebarPanel(
	h3("Xpert Magic"),
	p("An Shiny App to convert raw outputs from GeneXpert machine into tidy tabular data for further analysis.",
		style="font-size:12px"),
	p("version 4.20231130", style = "font-size:10px;"),
	hr(),

	fileInput("upload", "Select raw data in .csv", accept = c(".csv", ".txt")),
	actionButton("process", "Process", icon = icon("microchip"), style = "float:right;"),
	br(),
	br(),
	p("- developed by", a("Myo Minn Oo", href="https://myominnoo.com"),
		style="font-size:12px;float:right;"),
	br()
)

# body --------------------------------------------------------------------

body <-
	mainPanel(
		tabsetPanel(
			type = "tabs",
			tabPanel(
				"Data",
				icon = icon("table"),
				br(),
				wellPanel(
					selectInput("vars", "Select columns to download:",
											choices = "", multiple = TRUE),
					downloadButton("download", "Save Data"),
				),
				br(),
				verbatimTextOutput("info"),
				dataTableOutput("data")
			),
			tabPanel(
				"About",
				icon = icon("info-circle"),
				verbatimTextOutput("about1"),
				verbatimTextOutput("about2")
			)
		)
)


# ui all ------------------------------------------------------------------

ui <-
	fluidPage(
	title = "Xpert Magic",
	theme = bs_theme(version = 5),
	br(),
	sidebarLayout(sidebar, body)
)


# server ------------------------------------------------------------------

server <- function(input, output, session) {

	env <- reactiveValues()

	observeEvent(input$process, {
		req(input$upload)
		n <- 5

		withProgress(
			message = 'Calculation in progress',
			detail = 'This may take a while...', value = 0, {
				raw <- readLines(input$upload$datapath, skipNul = TRUE)
				incProgress(1/n)
				df <-
					tryCatch(
						suppressWarnings(process_data(raw)),
						error = function(e) NULL
					)
				if (is.null(df)) {
					shinyalert::shinyalert("Oops!", "Something went wrong.", type = "error")
				}
				incProgress(1/n)
			}
		)

		env$data <- df
		print(df)
		str(df)
		vars <- names(df)
		print(vars)
		updateSelectInput(session = session, "vars", choices = vars,
											selected = vars)
	})

	# Quantification Table ----------------------------------------------------
	output$data <- renderDataTable({
		datatable(env$data,
							rownames = FALSE,
							editable = TRUE,
							filter = list(position = "top", clear = FALSE),
							options = list(scrollX = TRUE))
	})

	# download ----------------------------------------------------------------
	output$download <- downloadHandler(
		filename = function() paste0(
			"genexpert_data_", format(Sys.time(), "%d-%b-%Y_%H.%M.%S"), ".xlsx"
		),
		content = function(file) rio::export(env$data[, input$vars], file)
	)


	output$about1 <- renderPrint({
		sessionInfo()
	})

	output$about2 <- renderPrint({
		sessioninfo::package_info()
	})
}

# run ---------------------------------------------------------------------

shinyApp(ui, server)
