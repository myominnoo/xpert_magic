

library(shiny)
library(bslib)
library(DT)
library(tidyverse)


options(shiny.maxRequestSize = 30*1024^2)

# functions ---------------------------------------------------------------

process_raw <- function(raw, include) {
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
 
  	## run a loop to process individual RESULT TABLE `irt`
  	irt <- lapply(prt, function(x) {
  		t <- raw[x]
  		t <- t[grepl(texts_included, t)]
  		t <- t[1:30]
  		t <- t[!is.na(t)]
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
    res <- cbind(res,
                 result = ifelse(pos, "Positive", "Negative"),
                 ct_hpv_16 = analyte[ct_val, which(analyte[1, ] == "HPV 16")[1]],
                 ct_hpv_18_45 = analyte[ct_val, which(analyte[1, ] == "HPV 18_45")[1]],
                 ct_p3 = analyte[ct_val, which(analyte[1, ] == "P3")[1]],
                 ct_p4 = analyte[ct_val, which(analyte[1, ] == "P4")[1]],
                 ct_p5 = analyte[ct_val, which(analyte[1, ] == "P5")[1]],
                 ct_sac = analyte[ct_val, which(analyte[1, ] == "SAC")[1]])
    if (ncol(res) %in% 9:10)
      data.frame(res)
  })

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
  output
}


# ui ----------------------------------------------------------------------


# sidebar -----------------------------------------------------------------

sidebar <- sidebarPanel(
  h3("Xpert Magic"), 
  p("A tool to convert raw data into ct values and HPV results", 
    style="font-size:12px"), 
  p("version 1.20221118", style = "font-size:10px;"), 
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
  
  df <- eventReactive(input$process, {
    req(input$upload)
    
    
    # Create a Progress object
    progress <- shiny::Progress$new()
    # Make sure it closes when we exit this reactive, even if there's an error
    on.exit(progress$close())

    progress$set(message = "Processing", value = 0)
    
    raw <- readLines(input$upload$datapath, skipNul = TRUE)
    str(input$include)
    raw <- suppressWarnings(process_raw(raw, input$include))
    vars <- names(raw)
    updateSelectInput(session = session, "vars", choices = vars, 
                      selected = vars) 
    
    # Number of times we'll go through the loop
    n <- 5
    for (i in 1:n) {
      # Each time through the loop, add another row of data. This is
      # a stand-in for a long-running computation.
      raw <- readLines(input$upload$datapath, skipNul = TRUE)
      raw <- suppressWarnings(process_raw(raw, input$include))
      vars <- names(raw)
      updateSelectInput(session = session, "vars", choices = vars,
                        selected = vars)
      # Increment the progress bar, and update the detail text.
      progress$inc(1/n, detail = paste("", round(i/n*100, 0), "% ... "))
      # Pause for 0.1 seconds to simulate a long computation.
      Sys.sleep(0.1)
    }
    
    raw
  })
  
  d <- reactiveValues(data = NULL, download = NULL)
  observe({d$data <- d$download <- df()})
  
# Quantification Table ----------------------------------------------------
  output$data <- renderDataTable({
    datatable(d$data,
              rownames = FALSE,
              editable = TRUE,
              filter = list(position = "top", clear = FALSE),
              options = list(scrollX = TRUE))
  })
  

# download ----------------------------------------------------------------
  output$download <- downloadHandler(
    filename = function() paste0(
      "hpv_genexpert_processed_at_", format(Sys.time(), "%d-%b-%Y_%H.%M.%S"), ".csv"
    ), 
    content = function(file) write.csv(
      d$download[, input$vars], file, row.names = FALSE
    )
  )
}

# run ---------------------------------------------------------------------

shinyApp(ui, server)
