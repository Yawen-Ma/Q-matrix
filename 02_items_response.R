
# Responses ---------------------------------------------------------------

# Ignore the games and levels, just pick up any 8 sequence from the same skill families
# 8 * 4 = 32

items_selection <- function(data) {
  result <- data %>%
    group_by(USER_ID, game_id, LEVEL) %>%
    # any 2 levels from each game
    arrange(USER_ID, game_id, LEVEL, ROUND_STARTED_AT_TIME) %>%
    # for the same level, keep the first engagement
    filter(row_number() == 1) %>%
    ungroup() %>%
    group_by(USER_ID, skill_family) %>%
    filter(row_number() %in% c(1:8)) %>%
    ungroup() %>%
    arrange(USER_ID, game_id, LEVEL, ROUND_STARTED_AT_TIME) %>%
    mutate(item = row_number())  
  
  return(result)
}

# For example, items_C5 should include 8 responses from each skill family, 
# while not everyone has engaged that much, thus they might not have 32 responses in total
# on average, they have engaged with 19.08 items, ranges from (1, 32) 
# with 1st quantile: 15 items; 3rd quantile: 22 items
items_C5 <- items_selection(C5)

items_C5 %>% group_by(USER_ID) %>% reframe(n()) %>% summary()

# Subset data who has engaged 32 items in total ---------------------------

subset <- function(data) {
  # Initial subsetting to users with exactly 32 entries
  sub_user <- data %>%
    group_by(USER_ID) %>%
    summarise(n = n()) %>%
    # find out users who have engaged all 32 items
    filter(n == 32) %>%
    pull(USER_ID) # Extracting USER_IDs with exactly 32 records
  
  # Filtering the main data based on the extracted USER_IDs
  sub_data <- data %>%
    filter(USER_ID %in% sub_user)
  
  # Transformations, mutating, and pivoting
  result <- sub_data %>%
    group_by(USER_ID, skill_family) %>%
    mutate(sequence = row_number(),  # Assumes data is already grouped and ordered as needed
           response_accuracy = as.integer(is_level_mastered),
           item = paste(skill_family, sequence, sep = "_")) %>%
    ungroup() %>%
    # Take out response accuracy at this stage for generating q-matrix
    select(USER_ID, item, response_accuracy) %>%
    # Restructure the dataset, user by items (32 columns for each item)
    pivot_wider(names_from = item, values_from = response_accuracy, id_cols = USER_ID)
  
  return(result)
}

# For example, the code can be used to filter responses from participants who engaged with 32 items
C5_subset <- subset(items_C5)
# Change the order of the columns so that the skill families following PA, ED, V, MC
C5_subset <- C5_subset[, c(1, # user_id
                                 2,5,6,7,8,10,11,12, # PA items 1:8
                                 3,4,9,17,18,21,32,33, # ED items 1:8
                                 13,14,22,23,24,25,26,27, # V items 1:8
                                 15,16,19,20,28,29,30,31)] # MC items 1:8

n_distinct(C5_subset$USER_ID) 

# 194 users 
