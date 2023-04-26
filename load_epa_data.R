library(dplyr)


load_epa_data = function(snapshot_id=1) {
  # Load EPA data
  epa.data = readRDS('data/df_data_12list.RDS')
  # Single snapshot
  epa.df = epa.data[[snapshot_id]]
  epa.df$case_cntl = NULL
  epa.df$year = NULL
  names(epa.df) = c('station_id', 'x', 'y', 'pm2_5')
  epa.df$station_id = NULL
  
  # Filter to remove duplicates
  epa.df = epa.df[!duplicated(epa.df[, c('x', 'y')]), ]
  
  # Convert to mean zero and standard deviation 1
  epa.df$pm2_5 = epa.df$pm2_5 %>% scale
  
  return(epa.df)
}


train_test_split = function(df, train_size=0.7) {
  # Train-test split data
  train_size = as.integer(train_size * nrow(epa.df))
  epa.train.idx = sample(1:nrow(epa.df), size=train_size)
  epa.test.idx = setdiff(1:nrow(epa.df), epa.train.idx)
  epa.train.df = epa.df[epa.train.idx, ]
  epa.test.df = epa.df[epa.test.idx, ]
  
  # Create X and y matrices
  X.train = epa.train.df[, names(epa.train.df) != 'pm2_5']
  X.test = epa.test.df[, names(epa.test.df) != 'pm2_5']
  y.train = epa.train.df$pm2_5
  y.test = epa.test.df$pm2_5
  
  # Package into training and testing data frames
  df.train = cbind(X.train, y.train)
  df.test = cbind(X.test, y.test)
  names(df.train) = c('x', 'y', 'pm2_5')
  names(df.test) = c('x', 'y', 'pm2_5') 
  
  return(list("X.train"=X.train,
              "y.train"=y.train,
              "X.test"=X.test,
              "y.test"=y.test,
              "df.train"=df.train,
              "df.test"=df.test))
}