

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
	texts_included <- c("Patient ID|Sample Type,|Test Result,|Analyte Name|HPV 16|HPV 18_45|P3,|P4,|P5,|SAC,|Start Time,")
	
	## if prt is not a list, convert it to a list
	if (!is.list(prt))
		prt <- list(prt)
	
	## run a loop to process individual RESULT TABLE `irt`
	irt <- lapply(prt[1:100], function(x) {
		t <- raw[x]
		t <- t[grepl(texts_included, t)]
		t <- t[!is.na(t)]
		hpv <- any(grepl("Test Result,HPV", t))
		
		if (hpv) {
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
			bin <- unlist(strsplit(grep("Test Result,HPV", t, value = TRUE), split = ","))[2]
			bin <- unlist(strsplit(bin, split = "\\|"))
			bin <- gsub("HPV 16 |HPV 18_45 |OTHER HR HPV ", "", bin)
			
			if (include == "type") res <- cbind(res, type = type)
			res <- cbind(
				res,
				result = ifelse(pos, "Positive", "Negative"),
				hpv_16 = bin[1], 
				hpv_18_45 = bin[2], 
				hpv_hr = bin[3], 
				ct_hpv_16 = analyte[ct_val, which(analyte[1, ] == "HPV 16")[1]],
				ct_hpv_18_45 = analyte[ct_val, which(analyte[1, ] == "HPV 18_45")[1]],
				ct_p3 = analyte[ct_val, which(analyte[1, ] == "P3")[1]],
				ct_p4 = analyte[ct_val, which(analyte[1, ] == "P4")[1]],
				ct_p5 = analyte[ct_val, which(analyte[1, ] == "P5")[1]],
				ct_sac = analyte[ct_val, which(analyte[1, ] == "SAC")[1]]
			) 
			
			if (ncol(res) == 13) {
				if (grepl("TSH", pid)) data.frame(res) else NA
			} else {NA}
		} else {NA}
		
	})
	
	# incProgress(1/n)
	
	eligible <- (is.na(irt))
	cat("\n ***** Removing", sum(eligible), "non-HPV Gates Study resuts for HPV *****\n")
	irt <- irt[!eligible]
	
	## combine all results into one table
	output <- do.call(rbind, irt)
	if (is.null(output)) {
		hpv <- data.frame()
	} else {
		output <- output %>%
			mutate(across(ct_hpv_16:ct_sac, ~ as.numeric(.x)))
		
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
		
		hpv <- output %>% 
			dplyr::mutate(start_time = lubridate::dmy(start_time), 
										PID = gsub("\"", "", PID)) %>% 
			dplyr::filter(!is.na(start_time))
	}
	
	# incProgress(1/n)
	
	
	## run a loop to process individual RESULT TABLE `irt`
	irt <- lapply(prt, function(x) {
		t <- raw[x]
		t <- t[grepl(texts_included, t)]
		t <- t[!is.na(t)]
		ctng <- any(grepl("Test Result,CT", t))
		
		if (ctng) {
			
			pid <- grep("Patient ID,", t, value = TRUE)
			pid <- gsub("\"|Patient ID,", "", pid)
			pid <- paste(pid, collapse = ", ")
			
			type <- unlist(strsplit(t[grepl("Sample Type,", t)], ","))[2]
			sttime <- unlist(strsplit(t[grepl("Start Time,", t)], ","))[2]
			sttime <- unlist(strsplit(sttime, " "))[1]
			
			bin <- unlist(strsplit(grep("Test Result,CT", t, value = TRUE), split = ","))[2]
			bin <- unlist(strsplit(bin, split = "\\|"))
			bin <- gsub("CT |NG ", "", bin)
			bin <- ifelse(bin == "NOT DETECTED", "NEG", 
										ifelse(bin == "DETECTED", "POS", NA))
			
			res <- cbind(
				start_time = sttime, PID = pid, 
				type = type, 
				ct = bin[1], 
				ng = bin[2]
			)
			
			if (ncol(res) == 5) {
				if (grepl("TSH", pid)) data.frame(res) else NA
			} else {NA}
		} else {NA}
			
	})
	
	
	# incProgress(1/n)
	
	
	eligible <- (is.na(irt))
	cat("\n ***** Removing", sum(eligible), "non-HPV Gates Study resuts for CTNG *****\n")
	irt <- irt[!eligible]
	
	## combine all results into one table
	output <- do.call(rbind, irt)
	if (is.null(output)) {
		ctng <- data.frame()
	} else {
		ctng <- output %>% 
			dplyr::mutate(start_time = lubridate::dmy(start_time), 
										PID = gsub("\"", "", PID)) %>% 
			dplyr::filter(!is.na(start_time))
	}
	
	list(hpv = hpv, ctng = ctng)
}
# process_raw(raw, include = "type")



		
	



# process_raw(raw, include = "type")



# ui ----------------------------------------------------------------------


# sidebar -----------------------------------------------------------------

sidebar <- sidebarPanel(
	h3("Xpert Magic"), 
	p("A tool to convert raw outputs of HPV & CTNG results from GeneXpert machine into tidy tabular data for further analysis.", 
		style="font-size:12px"), 
	p("version 3.20230206", style = "font-size:10px;"), 
	hr(), 
	
	fileInput("upload", "Select raw data in .csv", accept = ".csv"), 
	selectInput("include", "Include these data", 
							choices = c("Sample Type"="type", "None"="none"), selected = "type"), 
	actionButton("process", "Process", icon = icon("microchip")), 
	br(), 
	br(), 
	p("- developed by Myo Minn Oo", style="font-size:12px;float:right;"), 
	br()
)

# body --------------------------------------------------------------------

body <- mainPanel(
	tabsetPanel(
		type = "tabs", 
		tabPanel(
			"HPV", 
			icon = icon("table"), 
			br(), 
			wellPanel(
				selectInput("hpv_vars", "Select columns to download:", 
										choices = "", multiple = TRUE), 
				downloadButton("hpv_download", "Save HPV Data"), 
			), 
			br(), 
			verbatimTextOutput("hpv_info"), 
			dataTableOutput("hpv_data")
		), 
		tabPanel(
			"CTNG", 
			icon = icon("table"), 
			br(), 
			wellPanel(
				selectInput("ctng_vars", "Select columns to download:", 
										choices = "", multiple = TRUE), 
				downloadButton("ctng_download", "Save CTNG Data"), 
			), 
			br(), 
			verbatimTextOutput("ctng_info"), 
			dataTableOutput("ctng_data")
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
	
	env <- reactiveValues(hpv = NULL, ctng = NULL)
	
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
				hpv_vars <- names(df[["hpv"]])
				updateSelectInput(session = session, "hpv_vars", choices = hpv_vars,
													selected = hpv_vars)
				ctng_vars <- names(df[["ctng"]])
				updateSelectInput(session = session, "ctng_vars", choices = ctng_vars,
													selected = ctng_vars)
				incProgress(1/n)
			}
		)
		
		env$hpv <- df[["hpv"]]
		env$ctng <- df[["ctng"]]
	})

	# Quantification Table ----------------------------------------------------
	output$hpv_data <- renderDataTable({
		datatable(env$hpv,
							rownames = FALSE,
							editable = TRUE,
							filter = list(position = "top", clear = FALSE),
							options = list(scrollX = TRUE))
	})
	output$ctng_data <- renderDataTable({
		datatable(env$ctng,
							rownames = FALSE,
							editable = TRUE,
							filter = list(position = "top", clear = FALSE),
							options = list(scrollX = TRUE))
	})
	
	# download ----------------------------------------------------------------
	output$hpv_download <- downloadHandler(
		filename = function() paste0(
			"hpv_genexpert_processed_at_", format(Sys.time(), "%d-%b-%Y_%H.%M.%S"), ".xlsx"
		),
		content = function(file) rio::export(env$hpv[, input$hpv_vars], file)
	)
	output$ctng_download <- downloadHandler(
		filename = function() paste0(
			"ctng_genexpert_processed_at_", format(Sys.time(), "%d-%b-%Y_%H.%M.%S"), ".xlsx"
		),
		content = function(file) rio::export(env$ctng[, input$ctng_vars], file)
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