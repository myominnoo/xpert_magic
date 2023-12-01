


get_meta <- function(t, cols, vars) {
	## extract meta information about the test run
	l <- sapply(cols$cols, function(col) {
		d <- grep(pattern = paste0("^", col), x = t, value = TRUE,
							useBytes = TRUE) |>
			stringr::str_split(col, simplify = TRUE)
		d[d != ""] |>
			paste(collapse = " ")
	})
	l

	idx <- which(t == "Analyte Result"):which(t == "Detail")
	d <- t[idx]
	d <- d[!(d %in% c("Analyte Result", "Detail"))] |>
		stringr::str_split(",", simplify = TRUE) |>
		data.frame(check.names = FALSE)
	names(d) <- d[1, ]
	d <- d[-1, ]
	d

	rbind(l) |>
		data.frame(check.names = FALSE) |>
		cbind(d) |>
		setNames(c(gsub(",", "", cols$cols), vars$vars))
}


process_data <- function(raw) {
	# tryCatch({

		## remove empty vectors
		raw <- raw[which(raw != "")]

		## get start and end positions based on the value RESULT TABLE:
		start <- which(raw == "RESULT TABLE")
		# end <- which(raw == "Run History")
		end <- which(raw == "Melt Peaks")
		## generate serial position number from start to end
		if (length(start) == 1) {
			start_end <- start:ifelse(end < start, length(raw), end)
		} else {
			start_end <- mapply(`:`, start, end, SIMPLIFY = FALSE)
		}

		cols <- readxl::read_excel("cols.xlsx", sheet = "cols")
		vars <- readxl::read_excel("cols.xlsx", sheet = "vars")

		## get meta information
		if (!is.list(start_end)) {
			start_end <- list(start_end)
		}

		f <- lapply(start_end, function(se) {
			get_meta(raw[se], cols, vars)
		})

		do.call(rbind, f) |>
			data.frame(check.names = FALSE)
			# tidyr::unnest(cols = dplyr::everything())
	# }, error = function(e) message("Something unexpected happens: ", e))
}

## read data as line text: this will make separate lines to character vectors

# raw <- readLines("data/Xpert H 230522144713_2022.05.23_14.48.064.csv", skipNul = TRUE)
# process_data(raw = raw) |>
# 	str()
#
# raw <- readLines("data/Xpert H 230522144713_2022.05.23_14.48.06JKHK.csv", skipNul = TRUE)
# process_data(raw = raw) |>
# 	str()
#
# raw <- readLines("data/Xpert Raw.csv", skipNul = TRUE)
# process_data(raw = raw) |>
# 	str()
#
# raw <- readLines("data/873180744V87.0_2023.11.29_17.36.05.csv", skipNul = TRUE)
# process_data(raw = raw) |>
# 	str()


