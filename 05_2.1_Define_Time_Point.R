
# Define time point --------------------------------------------------

define_time_point <- function(items_cluster, sequence_start_time, sequence_end_time, months_per_period = 12) {
  # Generate breaks based on the specified number of months
  breaks <- seq(as.Date(sequence_start_time), as.Date(sequence_end_time), by = paste(months_per_period, "months"))
  
  # Determine the number of periods and create labels
  period_count <- length(breaks) - 1
  period_labels <- if (months_per_period == 12) {
    paste("Year", 1:period_count)  # Yearly labels
  } else {
    paste("Period", 1:period_count)  # Other intervals
  }
  
  items_cluster %>%
    mutate(
      # Convert time data
      Date = as.Date(round_started_at),
      # Create Time Period using dynamic slicing with labels
      Time_Period = cut(Date, breaks = breaks, labels = period_labels, include.lowest = TRUE))%>%
    arrange(USER_ID, round_started_at) %>%
    group_by(USER_ID, Time_Period) %>%
    mutate(
          # Assign item numbers within each time period for each user
          Item_at_time_point = 1:length(Time_Period)
          ) %>%
   ungroup()  
}

# Example:
# Monthly slicing
monthly_points <- define_time_point(items_C5, "2021-08-01", "2023-08-01", 1)
# Yearly slicing
yearly_points <- define_time_point(items_C5, "2021-08-01", "2023-08-01", 12)
# Quarterly slicing (every three months)
quarterly_points <- define_time_point(items_C5, "2021-08-01", "2023-08-01", 3)


