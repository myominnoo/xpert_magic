

library(shiny)
library(bslib)
library(DT)
library(tidyverse)
library(lubridate)


options(shiny.maxRequestSize = 30*1024^2)
# raw <- readLines("Xpert H 230522144713_2022.05.23_14.48.06JKHK.csv", skipNul = TRUE)
# raw <- readLines("Xpert Raw.csv", skipNul = TRUE)


# functions ---------------------------------------------------------------

process_raw <- function(raw, include, n = 5) {
	raw <- raw[raw != ""]
	
	## calculate position indexes for RESULT TABLE:
	## pri = position result index
	## prt = position result table
	pri <- c(grep("RESULT TABLE", raw), length(raw))
	prt <- mapply(`:`, pri[-length(pri)], pri[-1], 
								SIMPLIFY = TRUE)
	
	## create a character vector to remove irrelevant texts
	texts_included <- c("Patient ID|Sample Type,|Analyte Name|HPV 16|HPV 18_45|P3,|P4,|P5,|SAC,|Start Time,")
	
	## if prt is not a list, convert it to a list
	if (!is.list(prt))
		prt <- list(prt)
	# incProgress(1/n)
	
	## run a loop to process individual RESULT TABLE `irt`
	irt <- lapply(prt, function(x) {
		t <- raw[x]
		t <- t[grepl(texts_included, t)]
		t <- t[1:30]
		t <- t[!is.na(t)]
		
		if (length(t) >= 15) {
			pid <- unlist(strsplit(t[grepl("Patient ID,", t)], ","))[2]
			type <- unlist(strsplit(t[grepl("Sample Type,", t)], ","))[2]
			sttime <- unlist(strsplit(t[grepl("Start Time,", t)], ","))[2]
			sttime <- unlist(strsplit(sttime, " "))[1]
			
			
			analyte <- c("Analyte Name,", "HPV 16,", "HPV 18_45,",
									 "P3,", "P4,", "P5,", "SAC,")
			analyte <- paste0(analyte, collapse = "|")
			analyte <- t(do.call(rbind, strsplit(t[grepl(analyte, t)], ",")))
			
			pos <- any(analyte[4, 2:7] %in% "POS")
			
			ct_val <- which(analyte[, 1] == "Ct")[1]
			
			res <- cbind(start_time = sttime,
									 PID = pid)
			
			if (include == "type")
				res <- cbind(res, type = type)
			res <- cbind(
				res,
				result = ifelse(pos, "Positive", "Negative"),
				ct_hpv_16 = analyte[ct_val, which(analyte[1, ] == "HPV 16")[1]],
				ct_hpv_18_45 = analyte[ct_val, which(analyte[1, ] == "HPV 18_45")[1]],
				ct_p3 = analyte[ct_val, which(analyte[1, ] == "P3")[1]],
				ct_p4 = analyte[ct_val, which(analyte[1, ] == "P4")[1]],
				ct_p5 = analyte[ct_val, which(analyte[1, ] == "P5")[1]],
				ct_sac = analyte[ct_val, which(analyte[1, ] == "SAC")[1]]
			) 
			
			if (ncol(res) == 10) {
				if (grepl("TSH", pid)) data.frame(res) else NA
			} else {NA}
		} else {NA}
	})
	
	
	# incProgress(1/n)
	
	eligible <- (sapply(irt, is.null))
	cat("\n ***** Removing", sum(eligible), "non-HPV resuts *****\n")
	irt <- irt[!eligible]
	
	## combine all results into one table
	output <- do.call(rbind, irt)
	output <- output %>%
		mutate(across(ct_hpv_16:ct_sac, ~ as.numeric(.x)))
	# output[, -c(1:4)] <- apply(output[, -c(1:4)], 2, as.numeric)
	
	# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4355361/
	output$vl_hpv_16 <- output$ct_hpv_16 / (output$ct_sac / 2)
	output$vl_hpv_18_45 <- output$ct_hpv_18_45 / (output$ct_sac / 2)
	output$vl_hpv_p3 <- output$ct_p3 / (output$ct_sac / 2)
	output$vl_hpv_p4 <- output$ct_p4 / (output$ct_sac / 2)
	output$vl_hpv_p5 <- output$ct_p5 / (output$ct_sac / 2)
	names(output)
	output[is.na(output) | sapply(output, is.nan)] <- 0
	output[which(output$result == "Negative"), c(
		"vl_hpv_16", "vl_hpv_18_45", "vl_hpv_p3", "vl_hpv_p4", "vl_hpv_p5"
	)] <- 0
	
	output %>% 
		dplyr::mutate(start_time = lubridate::dmy(start_time)) %>% 
		dplyr::filter(!is.na(start_time))
}



# process_raw(raw, include = "type")



# ui ----------------------------------------------------------------------


# sidebar -----------------------------------------------------------------

sidebar <- sidebarPanel(
	h3("Xpert Magic"), 
	p("A tool to convert raw data into ct values and HPV results", 
		style="font-size:12px"), 
	p("version 2.20230118", style = "font-size:10px;"), 
	hr(), 
	
	fileInput("upload", "Select raw data in .csv", accept = ".csv"), 
	selectInput("include", "Include these data", 
							choices = c("Sample Type"="type", "None"="none"), selected = "type"), 
	actionButton("process", "Process", icon = icon("microchip")), 
	hr(), 
	hr(), 
	selectInput("vars", "Select columns to download:", 
							choices = "", multiple = TRUE), 
	downloadButton("download", "Save Data"), 
	br(), 
	br(), 
	p("- developed by Myo Minn Oo", style="font-size:12px;float:right;"), 
	br()
)

# body --------------------------------------------------------------------

body <- mainPanel(
	tabsetPanel(
		type = "pills", 
		tabPanel(
			"Quantification Table", 
			icon = icon("table"), 
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

ui <- fluidPage(
	title = "Xpert Magic", 
	theme = bs_theme(version = 5), 
	br(), 
	sidebarLayout(sidebar, body)
)


# server ------------------------------------------------------------------

server <- function(input, output, session) {
	
	env <- reactiveValues(data = NULL)
	
	observeEvent(input$process, {
		req(input$upload)
		n <- 5
		
		withProgress(
			message = 'Calculation in progress',
			detail = 'This may take a while...', value = 0, {
				df <- readLines(input$upload$datapath, skipNul = TRUE)
				incProgress(1/n)
				df <- suppressWarnings(process_raw(df, input$include, n))
				incProgress(1/n)
				vars <- names(df)
				updateSelectInput(session = session, "vars", choices = vars,
													selected = vars)
				incProgress(1/n)
			}
		)
		
		env$data <- df
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
			"hpv_genexpert_processed_at_", format(Sys.time(), "%d-%b-%Y_%H.%M.%S"), ".xlsx"
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