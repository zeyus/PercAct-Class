"Utilities for working with Motion capture data
from QTM.
"

library(stringi)
library(assertthat)
library(stats)
library(XML)
library(dplyr)
library(tidyr)
# adapted from
# https://github.com/zeyus/QTM_Bela_Sonification/blob/main/scripts/step_01_import_raw_data.R

#' This function processes a QTM TSV file and reads it into a dataframe.
#'
#'   You will get both the trajectory data (dataframe) as well as metadata.
#'
#'   Args:
#'     `data_file` (character): The path to the QTM TSV file.
#'     `condition_regex` (character): A regex to extract the condition
#'        from the file name.
#'     `null_value` (character): The value in the files that represents a
#'       missing measurement.
#'    `only_data` (bool): If TRUE, only the data will be returned.
#'   `long_format` (bool): If TRUE, the data will be returned in long format.
#'
#'   Returns:
#'     * A dataframe with the trajectory data. If `only_data` is TRUE
#'       otherwise, if `only_data` is FALSE
#'     * A list with two elements:
#'       * `data`: The trajectory data as a dataframe.
#'       * `metadata`: A list with the metadata.
#'           The `metadata` contains the following elements:
#'             * `condition`: The condition of the data.
#'             * `frequency`: The frequency (Hz) of the data.
#'             * `marker_count`: The number of markers.
#'             * `frame_count`: The number of frames.
#'             * `marker_names`: The names of the markers.
process_qtm_tsv <- function(
    data_file,
    condition_regex = "group[0-9]+_([^_\\.]+).*",
    null_value = "NA",
    only_data = FALSE,
    long_format = FALSE) {
  message(paste("Processing", data_file))

  # get the condition
  cond <- stri_match_last_regex(data_file, condition_regex)[2]
  # make condition lowercase
  cond <- tolower(cond)
  message(paste("Condition:", cond))

  # read in the data
  dat <- readLines(data_file)

  # now get other relevant metadata
  # get frequency
  freq <- stri_match_first_regex(dat, "^FREQUENCY.*")
  freq <- freq[!is.na(freq)]

  # get marker count
  marker_count <- stri_extract_first_regex(dat, "^NO_OF_MARKERS.*")
  marker_count <- marker_count[!is.na(marker_count)]
  marker_count_value <- as.integer(stri_split_fixed(marker_count, "\t")[[1]][2])

  # get frame count
  frame_count <- stri_extract_first_regex(dat, "^NO_OF_FRAMES.*")
  frame_count <- frame_count[!is.na(frame_count)]
  frame_count_value <- as.integer(stri_split_fixed(frame_count, "\t")[[1]][2])

  # get marker names
  marker_names <- stri_extract_first_regex(dat, "^MARKER_NAMES.*")
  marker_names <- marker_names[!is.na(marker_names)]
  marker_names_values <- stri_split_fixed(
    marker_names,
    "\t"
  )[[1]][1:marker_count_value + 1]
  message("File information:")
  message(paste(freq, "Hz", "frequency"))
  message(paste(marker_count_value, "markers"))
  message(paste(frame_count_value, "frames"))

  # now just keep the tracking information
  dat <- stri_extract_first_regex(dat, "^[0-9]+\t.*")
  dat <- dat[!is.na(dat)]

  # now ensure the number of frames is correct
  # this is the number of lines in the data
  assertthat::assert_that(
    length(dat) == frame_count_value,
    msg = paste(
      "Number of frames is not correct, found:",
      length(dat),
      "expected:",
      frame_count_value
    )
  )

  # ensure the number of markers is correct
  # this is 3 columns per marker, plus index and time
  num_found_markers <- length(stri_split_fixed(dat[[1]], "\t")[[1]])
  assertthat::assert_that(
    num_found_markers == marker_count_value * 3 + 2,
    msg = paste(
      "Number of markers is not correct, found:",
      num_found_markers,
      "expected:",
      marker_count_value * 3 + 2
    )
  )

  message(paste("File has", length(dat), "frames"))
  message(paste("File has", marker_count_value, "markers"))

  message("Creating data frame...")

  col_names <- c(
    "index",
    "elapsed_time",
    paste0(
      rep(marker_names_values, each = 3),
      c("_x", "_y", "_z")
    )
  )


  # now we need to create a data frame

  # split each line by tab
  dat <- stri_split_fixed(dat, "\t", simplify = TRUE)
  # set the column names
  colnames(dat) <- col_names
  # and then convert to a data frame
  dat <- as_tibble(dat)
  # now we need to replace the null values with NA
  dat[dat == null_value] <- NA
  # and convert to numeric
  dat <- mutate_all(dat, as.numeric)

  if (long_format) {
    # let's make the data long format so we can
    # easily group by subject or marker
    dat <- dat %>%
      pivot_longer(
        cols = contains("_x") | contains("_y") | contains("_z"),
        cols_vary = "slowest",
        names_to = "marker",
        values_to = "value"
      ) %>%
      mutate(
        subject = stri_replace_first_regex(marker, "^([AB])_.*", "$1"),
        axis = stri_extract_last_regex(marker, "[xyz]$"),
        marker = stri_replace_first_regex(
          marker,
          "^[AB]_([a-zA-Z_]+)_[xyz]$",
          "$1"
        )
      )
    # move axes to columns
    dat <- dat %>%
      pivot_wider(
        names_from = axis,
        values_from = value
      )
    dat$marker <- factor(dat$marker)
  }

  metadata <- list(
    condition = cond,
    frequency = freq,
    marker_count = marker_count,
    frame_count = frame_count,
    marker_names = marker_names
  )

  # return the data and metadata
  if (only_data) {
    return(dat)
  }
  return(list(data = dat, metadata = metadata))
}


#' This function does a linear gap fill for a vector.
#'
#'   Only columns with at least 2 non-NA values will be gap filled.
#'
#'   Args:
#'     x (vector): The vector to gap fill.
#'
#'   Returns:
#'     The gap filled vector.
gap_fill_linear <- function(x) {
  # get the indices of the NA values
  na_indices <- which(is.na(x))
  # get the indices of the non-NA values
  non_na_indices <- which(!is.na(x))

  if (length(na_indices) == 0) {
    # if there are no NA values, just return the vector
    return(x)
  }

  if (length(non_na_indices) < 2) {
    # if there are less than 2 non-NA values, we can't do a linear interpolation
    return(x)
  }
  # get the values of the non-NA indices
  non_na_values <- x[non_na_indices]
  # get the values of the NA indices
  na_values <- x[na_indices]
  # now we can do a linear interpolation
  na_values <- approx(
    x = non_na_indices,
    y = non_na_values,
    xout = na_indices,
    method = "linear"
  )$y
  # now we can replace the NA values with the interpolated values
  x[na_indices] <- na_values
  # and return the vector
  return(x)
}


#' This function reads a QTM QML label file.
#'
#'  Args:
#'   xmlpath (str): The path to the QML file.
#'
#' Returns:
#'  A dataframe with the trajectory names and colors.
read_trajectory_labels <- function(xmlpath) {
  # load the XML file
  xmlfile <- xmlParse(xmlpath)

  # get the labels, which are in the following format:
  # <QTM_Label_List_Ver_1.00>
  #     <Trajectories>
  #         <Trajectory>
  #             <Name>A_head_top</Name>
  #             <Color R="0" G="147" B="0"/>
  #         </Trajectory>
  #     </Trajectories>
  # </QTM_Label_List_Ver_1.00>

  # get the trajectory names
  traj_names <- xpathSApply(xmlfile, "//Trajectory/Name", xmlValue)

  # get the trajectory colors
  traj_colors <- xpathSApply(xmlfile, "//Trajectory/Color", xmlAttrs)

  # convert the colors to hex
  traj_colors <- rgb(
    as.numeric(traj_colors[1, ]),
    as.numeric(traj_colors[2, ]),
    as.numeric(traj_colors[3, ]),
    alpha = 255,
    maxColorValue = 255
  )


  # combine the names and colors into a data frame
  traj_labels <- data.frame(
    traj_names,
    traj_colors,
    stringsAsFactors = FALSE
  )

  return(traj_labels)
}
