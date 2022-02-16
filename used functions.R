library(rpart)
library(imputeTS)
library(FNN)
library(class)
library(kknn)
library(ggplot2)
library(caret)
#Noise and imputation methods

noise <- function(clean_data,target,noise_level){
  
  to_noise <- as.integer(noise_level*nrow(clean_data))
  noised <- 0
  real_index <- c(1:nrow(clean_data))
  while(noised < to_noise){
    center <- sample(real_index,1)
    frame  <- sample(c(0:7),1)
    noise_id <- c((center-as.integer(frame/2)):(center+as.integer(frame/2)))
    noise_id <- noise_id[which(noise_id %in% real_index)]
    clean_data[noise_id,target] <- NA
    noised <- sum(is.na(clean_data))
  }
  return(list(noise_id = 0,noised_data = clean_data ))
}





mean_imputation <- function(noised_data,reference,target){
  
  only_noised <- noised_data[which(is.na(noised_data[,target])),]
  imputed <- na_mean(noised_data)
  imputed$uncertainty <- 1
  
  mean <- round(mean(noised_data[complete.cases(noised_data),target]))
  
  mean_vec <- rep(mean,nrow(noised_data[complete.cases(noised_data),]))
  #err <- (imputed[which(rownames(imputed) %in% rownames(only_noised)),target] - reference[which(rownames(reference) %in% rownames(only_noised)),target])^2
  err_count <- (mean_vec - noised_data[complete.cases(noised_data),target])^2
  rmse_norm <- sqrt(mean(err_count))/sqrt(max(err_count))
  mean_un <- 1-rmse_norm
  imputed[which(rownames(imputed) %in% rownames(only_noised)),"uncertainty"] <- mean_un
  return(list(imputed_data = imputed,error = 0))
}

#TIME-Eknn
learning_imputation <- function(noised_data,historic_length,historic_variables,reference,to_remove,target){
  
  noised_data$date <- as.double(noised_data$date)
  base <- noised_data[complete.cases(noised_data),]
  
  only_noised <- noised_data[which(is.na(noised_data[,target])),]
  
  
  f_eknn <- EkNNval_certain_time_regression(base[,-which(names(base) %in% c(target_name,"uncertainty"))],
                                            base[,target_name],
                                            base[,-which(names(base) %in% c(target_name,"uncertainty"))],
                                            10,
                                            1,
                                            base[,target_name],0,max(base[,target_name]))
  
  
  
  
  
  f_eknn_noised <- EkNNval_certain_time_regression(base[,-which(names(base) %in% c(target_name,"uncertainty"))],
                                                   base[,target_name],
                                                   only_noised[,-which(names(base) %in% c(target_name,"uncertainty"))],
                                                   10,
                                                   1,
                                                   NULL,
                                                   0,max(base[,target_name]))
  
  noised_data$uncertainty <- 1
  for (i in 1:nrow(only_noised) ) {
    noised_data[which(rownames(noised_data) %in% rownames(only_noised)[i]),target_name] <- round(f_eknn_noised$ypred[i])
    noised_data[which(rownames(noised_data) %in% rownames(only_noised)[i]),"uncertainty"] <-  f_eknn_noised$m[i,(f_eknn_noised$ypred[i] + 1)]
    
  }
  
  return(list(imputed_data = noised_data, error = f_eknn$err))
  
}

#LOCF
mean_imputation <- function(noised_data,target,reference,beta) {
  
  idx <- c(1:nrow(noised_data))
  noised_data <- cbind(noised_data,idx)
  err <- data.frame(err = 0,dist = 0)
  err <- err[-1,]
  noised_data$uncertainty <- 1
  labelled_data <- noised_data[complete.cases(noised_data),]
  for (i in 1:nrow(noised_data)) {
    if(is.na(noised_data[i,target])){
      past <- labelled_data[which(labelled_data$idx < noised_data[i,"idx"]),]
      past <- past[order(past$idx,decreasing = T),]
      past <- past[1,]
      past_distance <- noised_data[i,"idx"] - past[1,"idx"]
      total_distance <- noised_data[i,"idx"] - past[1,"idx"]
      #past_coef <- 1 - (past_distance/total_distance)
      noised_data[i,target] <- past[1,target]
      #noised_data[i,target] <- round(sum(past_coef*past[,target])) 
      noised_data[i,"uncertainty"] <- exp(-beta*min(past_distance))
      gap <- (noised_data[i,target] - reference[i,target])^2
      err <- rbind(err,cbind(gap,dist))
      #gap <- (noised_data[i,target] - reference[i,target])^2
      #dist <- futur[1,"idx"] - past[1,"idx"]
      #err <- rbind(err,cbind(gap,dist))
    }
  }
  return(list(imputed_data = noised_data[complete.cases(noised_data),-which(names(noised_data) %in% "idx")],error = err$gap))
}

#CMA
weighted_evidential_imputation <- function(noised_data,target,reference,beta) {
  
  idx <- c(1:nrow(noised_data))
  noised_data <- cbind(noised_data,idx)
  err <- data.frame(err = 0,dist = 0)
  err <- err[-1,]
  noised_data$uncertainty <- 1
  labelled_data <- noised_data[complete.cases(noised_data),]
  for (i in 1:nrow(noised_data)) {
    if(is.na(noised_data[i,target])){
      past <- labelled_data[which(labelled_data$idx < noised_data[i,"idx"]),]
      futur <- labelled_data[which(labelled_data$idx > noised_data[i,"idx"]),]
      past <- past[order(past$idx,decreasing = T),]
      futur <- futur[order(futur$idx),]
      past <- past[1:5,]
      futur <- futur[1:5,]
      #total_distance <- futur[nrow(futur),"idx"]-past[nrow(past),"idx"]
      past_distance <- noised_data[i,"idx"] - past[,"idx"]
      futur_distance <- futur[,"idx"] - noised_data[i,"idx"]
      total_distance <- max(past_distance) + max(futur_distance) 
      past_coef <-  1 - (past_distance/total_distance)
      futur_coef <- 1 - (futur_distance/total_distance)
      norm <- sum(past_coef) + sum(futur_coef)
      past_coef <- past_coef/norm
      futur_coef <- futur_coef/norm
      #print(futur_coef)
      #print(past_coef)
      #past_coef <- rev(past_coef)
      #futur_coef <-rev(futur_coef)
      #print(total_distance)
      #print(futur_distance)
      #print(past_distance)
      #print(sum(past_coef))
      #print(sum(futur_coef))
      print(sum(past_coef)+sum(futur_coef))
      noised_data[i,target] <- round(sum(past_coef*past[,target]) + sum(futur_coef*futur[,target]))
      distance <- c(futur_distance,past_distance)
      noised_data[i,"uncertainty"] <- exp(-beta*min(distance))
      gap <- (noised_data[i,target] - reference[i,target])^2
      dist <- futur[1,"idx"] - past[1,"idx"]
      err <- rbind(err,cbind(gap,dist))
      #gap <- (noised_data[i,target] - reference[i,target])^2
      #dist <- futur[1,"idx"] - past[1,"idx"]
      #err <- rbind(err,cbind(gap,dist))
    }
  }
  return(list(imputed_data = noised_data[complete.cases(noised_data),-which(names(noised_data) %in% "idx")],error = err$gap))
}


error_temporal <- function(labelled_data,target){
  err <- c()
  for (i in 2:(nrow(labelled_data)-1)) {
    fut <- labelled_data[i+1,target]
    past <- labelled_data[i-1,target]
    total_distance <- labelled_data[i+1,"idx"]-labelled_data[i-1,"idx"]
    past_distance <- labelled_data[i,"idx"] - labelled_data[i-1,"idx"]
    futur_distance <- labelled_data[i+1,"idx"] - labelled_data[i,"idx"]
    past_coef <- 1 - (past_distance/total_distance)
    futur_coef <- 1 - (futur_distance/total_distance)
    expectation <- round(past_coef*(labelled_data[i-1,target]) +  futur_coef*(labelled_data[i+1,target]))
    error <- (labelled_data[i,target] - expectation)^2
    err <- c(err,error)
  }
  
  return(err)
}


history_creation <- function(data,target,historic_length) {
  history_column_index <- c()
  for (i in 1:historic_length) {
    data[,ncol(data)+1] <- NA
    colnames(data)[ncol(data)] <- paste(paste(colnames(data)[target],"lags_",sep = "_"),i,sep = "_")
    history_column_index <- c(history_column_index,ncol(data))
  }
  return(list(data = data,history_column_index = history_column_index))
}


lags_fill <- function(data,filling_target,target,history_column_index,historic_length,prediction_horizon) {
  for (i in (prediction_horizon+historic_length):nrow(data)) {
    for (j in 1:historic_length) {
      idx <- (j + prediction_horizon) - 1
      data[i,history_column_index[j]] <- data[i - idx,filling_target]
      #if(filling_target == target){
      #  print(data[i,"uncertainty"])
      #  print(data[i-idx,"uncertainty"])
      #  data[i,"uncertainty"] <- data[i,"uncertainty"] * data[i - idx,"uncertainty"]
      #}
      
    }
  }
  return(data[complete.cases(data),])
}

history <- function(data,historic_variables,historic_length,target,prediction_horizon){
  if(historic_length == 0){return(data)}
  for(to_fill in historic_variables){
    data <- history_creation(data,to_fill,historic_length)
    history_column_index <- data$history_column_index
    data <- lags_fill(data$data,to_fill,target,history_column_index,historic_length,prediction_horizon)
  }
  return(data)
}

history_creation <- function(data,target,historic_length) {
  history_column_index <- c()
  for (i in 1:historic_length) {
    data[,ncol(data)+1] <- NA
    colnames(data)[ncol(data)] <- paste(paste(colnames(data)[target],"lags_",sep = "_"),i,sep = "_")
    history_column_index <- c(history_column_index,ncol(data))
  }
  return(list(data = data,history_column_index = history_column_index))
}


lags_fill2 <- function(data,filling_target,target,history_column_index,historic_length) {
  for (i in (historic_length+1):nrow(data)) {
    for (j in 1:historic_length) {
      data[i,history_column_index[j]] <- data[i - j,filling_target]
      
    }
  }
  return(data[complete.cases(data),])
}

history2 <- function(data,historic_variables,historic_length,target){
  if(historic_length == 0){return(data)}
  for(to_fill in historic_variables){
    data <- history_creation(data,to_fill,historic_length)
    history_column_index <- data$history_column_index
    data <- lags_fill2(data$data,to_fill,target,history_column_index,historic_length)
  }
  return(data)
}


#Evidential K-nearest neighbors

#EkNN regression Conjunctive fusion
EkNNval_uncertain_regression <- function(xtrain,ytrain,y_certainty,xtst,K,ytst=NULL,minset,maxset){
  
  add <- max(ytrain)*0.15
  add <- as.integer(add)
  maxset <- max(ytrain)+add
  omega <- seq(minset,maxset)
  #Train/test samples to matrix - Labels to integer
  xtst<-as.matrix(xtst)
  
  xtrain<-as.matrix(xtrain)
  
  #If NULL initialize parameters Alpha gamma
  
  #Napp Train size
  Napp<-nrow(xtrain)
  
  ytrain <- ytrain + 1
  #Class number
  M<- length(omega)
  #Test sample size
  N<-nrow(xtst)
  
  #nneighbors for test samples from train samples
  knn<-get.knnx(xtrain, xtst, k=K)
  
  #Index matrix
  is<-t(knn$nn.index)
  
  #Distance matrix
  ds<-t(knn$nn.dist)
  
  
  #Certainty matrix
  labels_mass <- t(knn$nn.dist)
  
  for (i in 1:nrow(labels_mass)) {
    for (j in 1:ncol(labels_mass)) {
      labels_mass[i,j] <- y_certainty[is[i,j]]
      
    }
    
  }
  #Mass matrix
  m = rbind(matrix(0,M,N),rep(1,N))
  #Mass for each neighbor, combination with dempster rule
  for(i in 1:N){
    for(j in 1:K){
      m1 <- rep(0,M+1)
      m1[ytrain[is[j,i]]] <- (exp(-ds[j,i]) * 0.999) * labels_mass[j,i]
      
      m1[M+1] <- 1 - m1[ytrain[is[j,i]]]
      m[1:M,i] <- m1[1:M]*m[1:M,i] + m1[1:M]*m[M+1,i] + m[1:M,i]*m1[M+1]
      m[M+1,i] <- m1[M+1] * m[M+1,i]
      m<-m/matrix(colSums(m),M+1,N,byrow=TRUE)
      
    }
  }
  
  #Mass matrix
  m<-t(m)
  #Pignistic Transformation
  for (i in 1:nrow(m)) {
    
    m[i,1:(ncol(m)-1)] <- m[i,1:(ncol(m)-1)] + ( m[i,ncol(m)] * 1/(ncol(m)-1))
    
  }
  
  ypred <- rep(0,N)
  
  #Pignistic expectation
  for (i in 1:N) {
    sum <- omega * m[i,1:M]
    ypred[i] <- sum(sum)
  }
  
  
  
  if(!is.null(ytst)) err<-sqrt(mean((ypred - ytst)^2)) else err<-NULL
  
  return(list(m=m[,-ncol(m)],ypred=ypred,err=err))
  
}


EkNNval_certain_time_regression <- function(xtrain,ytrain,xtst,K,date,ytst=NULL,minset,maxset){
  
  omega <- seq(minset:maxset)
  
  #Train/test samples to matrix - Labels to integer
  xtst<-as.matrix(xtst)
  
  xtrain<-as.matrix(xtrain)
  
  ytrain <- ytrain+1
  
  #If NULL initialize parameters Alpha gamma
  
  #Napp Train size
  Napp<-nrow(xtrain)
  
  #Class number
  M <- length(omega)
  
  #Test sample size
  N<-nrow(xtst)
  
  #nneighbors for test samples from train samples
  knn<-get.knnx(xtrain[,date], xtst[,date], k=K)
  
  #Index matrix
  is<-t(knn$nn.index)
  
  #Distance matrix
  ds<-t(knn$nn.dist)
  
  #for (i in 1:N) {
  #  ds[,i] <- scale(ds[,i],scale = T,center = F)
  
  #}
  
  #Certainty matrix
  labels_mass <- t(knn$nn.dist)
  
  #Mass matrix
  m = rbind(matrix(0,M,N),rep(1,N))
  
  #Mass for each neighbor, combination with dempster rule
  for(i in 1:N){
    for(j in 1:K){
      m1 <- rep(0,M+1)
      m1[ytrain[is[j,i]]] <- 0.95 * exp(-0.01*ds[j,i])
      
      
      m1[M+1] <- 1 - m1[ytrain[is[j,i]]]
      m[1:M,i] <- m1[1:M]*m[1:M,i] + m1[1:M]*m[M+1,i] + m[1:M,i]*m1[M+1]
      m[M+1,i] <- m1[M+1] * m[M+1,i]
      m<-m/matrix(colSums(m),M+1,N,byrow=TRUE)
    }
  }
  
  #Mass matrix
  m<-t(m)
  mf <- m
  #Pignistic Transformation
  for (i in 1:nrow(m)) {
    
    m[i,1:(ncol(m)-1)] <- m[i,1:(ncol(m)-1)] + ( m[i,ncol(m)] * 1/(ncol(m)-1))
    
  }
  
  ypred <- rep(0,N)
  
  #Pignistic expectation
  for (i in 1:N) {
    sum <- omega * m[i,1:M]
    ypred[i] <- sum(sum)
  }
  
  if(!is.null(ytst)) err<-sqrt(mean((ypred - ytst)^2)) else err<-NULL
  if(!is.null(ytst)) m_err<-sqrt(max((ypred - ytst)^2)) else m_err<-NULL
  
  return(list(m=mf,m_err = m_err,ypred=ypred,err=err))
  
}


EkNNval_certain_regression <- function(xtrain,ytrain,y_certainty,xtst,K,ytst=NULL,minset,maxset){
  
  add <- as.integer(max(ytrain))*0.15
  maxset <- max(ytrain)+add
  omega <- seq(minset:maxset)
  #Train/test samples to matrix - Labels to integer
  xtst<-as.matrix(xtst)
  
  xtrain<-as.matrix(xtrain)
  
  ytrain <- ytrain+1
  
  #If NULL initialize parameters Alpha gamma
  
  #Napp Train size
  Napp<-nrow(xtrain)
  
  #Class number
  M <- length(omega)
  
  #Test sample size
  N<-nrow(xtst)
  
  #nneighbors for test samples from train samples
  knn<-get.knnx(xtrain, xtst, k=K)
  
  #Index matrix
  is<-t(knn$nn.index)
  
  #Distance matrix
  ds<-t(knn$nn.dist)
  
  
  #Certainty matrix
  labels_mass <- t(knn$nn.dist)
  
  for (i in 1:nrow(labels_mass)) {
    for (j in 1:ncol(labels_mass)) {
      labels_mass[i,j] <- y_certainty[is[i,j]]
      
    }
    
  }
  #Mass matrix
  m = rbind(matrix(0,M,N),rep(1,N))
  #Mass for each neighbor, combination with dempster rule
  for(i in 1:N){
    for(j in 1:K){
      m1 <- rep(0,M+1)
      m1[ytrain[is[j,i]]] <- exp(-ds[j,i]) * 0.999
      
      m1[M+1] <- 1 - m1[ytrain[is[j,i]]]
      m[1:M,i] <- m1[1:M]*m[1:M,i] + m1[1:M]*m[M+1,i] + m[1:M,i]*m1[M+1]
      m[M+1,i] <- m1[M+1] * m[M+1,i]
      m<-m/matrix(colSums(m),M+1,N,byrow=TRUE)
    }
  }
  
  #Mass matrix
  m<-t(m)
  #Pignistic Transformation
  for (i in 1:nrow(m)) {
    
    m[i,1:(ncol(m)-1)] <- m[i,1:(ncol(m)-1)] + ( m[i,ncol(m)] * 1/(ncol(m)-1))
    
  }
  
  ypred <- rep(0,N)
  
  #Pignistic expectation
  for (i in 1:N) {
    sum <- omega * m[i,1:M]
    ypred[i] <- sum(sum)
  }
  
  if(!is.null(ytst)) err<-sqrt(mean((ypred - ytst)^2)) else err<-NULL
  
  return(list(m=m[,-ncol(m)],ypred=ypred,err=err))
  
}




#Evaluation

evaluation_naive <- function(imputed_data,target_name,reference,start_day,prediction_horizon,n_neighbor){
  #print("NAIVE")
  naive_acc <- c()
  for (i in start_day:(nrow(imputed_data)-prediction_horizon)) {
    train_set <- imputed_data[1:i,]
    test_set <- imputed_data[i+prediction_horizon,]
    test_set[,target_name] <- reference[i+prediction_horizon,target_name]
    prediction <- mean(train_set[(nrow(train_set)-14):nrow(train_set),target_name])
    p_vector <- rep(prediction,nrow(test_set))
    acc <- sqrt(mean( (p_vector - test_set[,target_name])^2 ))
    naive_acc <- c(naive_acc,prediction)
  }
  return(naive_acc)
}



evaluation_eknn <- function(imputed_data,target_name,date,reference,start_day,prediction_horizon,n_neighbor){
  #print("EKNN REGRESSION")
  eknn_certain <- c()
  for (i in start_day:(nrow(imputed_data)-prediction_horizon)) {
    train_set <- imputed_data[1:i,]
    #set <- variable_selection_eknn(train_set,n_neighbor,target_name,prediction_horizon)
    #train_set <- train_set[,set]
    #train_set[,-which(names(train_set) == c("uncertainty",target_name) )] <- scale(train_set[,-which(names(train_set) == c("uncertainty",target_name) )],scale = T,center = T)
    
    test_set <- imputed_data[i+prediction_horizon,]
    test_set[,target_name] <- reference[i+prediction_horizon,target_name]
    #test_set <- reference[(i+1):(i+prediction_horizon),set]
    #test_set[,-which(names(test_set) == target_name )] <- scale(test_set[,-which(names(test_set) == target_name )],scale = T,center = T)
    #print(sum(is.na(train_set)))
    #print(sum(is.na(test_set)))
    
    #print(i)
    f_eknn <- EkNNval_uncertain_regression(train_set[,-which(names(train_set) %in% c(target_name,"uncertainty","date"))],
                                           train_set[,target_name],
                                           train_set[,"uncertainty"],
                                           test_set[,-which(names(test_set) %in% c(target_name,"uncertainty","date"))],
                                           n_neighbor,
                                           test_set[,target_name],0,5000)
    
    eknn_certain <- c(eknn_certain,f_eknn$ypred)
    nas <- sum(is.na(f_eknn$ypred))
    if(nas > 0){
      print("ISSUE EKNN")
      data <- NULL
    }
    
  }
  
  return(eknn_certain)
}


evaluation_eknnc <- function(imputed_data,target_name,date,reference,start_day,prediction_horizon,n_neighbor){
  #print("EKNN CERTAIN REGRESSION")
  eknn_certain <- c()
  print("début")
  for (i in start_day:(nrow(imputed_data)-prediction_horizon)) {
    train_set <- imputed_data[1:i,]
    #set <- variable_selection_eknn(train_set,n_neighbor,target_name,prediction_horizon)
    #train_set <- train_set[,set]
    #train_set[,-which(names(train_set) == c("uncertainty",target_name) )] <- scale(train_set[,-which(names(train_set) == c("uncertainty",target_name) )],scale = T,center = T)
    
    test_set <- imputed_data[i+prediction_horizon,]
    test_set[,target_name] <- reference[i+prediction_horizon,target_name]
    #test_set <- reference[(i+1):(i+prediction_horizon),set]
    #test_set[,-which(names(test_set) == target_name )] <- scale(test_set[,-which(names(test_set) == target_name )],scale = T,center = T)
    #print(sum(is.na(train_set)))
    #print(sum(is.na(test_set)))
    
    #print(i)
    f_eknn <- EkNNval_certain_regression(train_set[,-which(names(train_set) %in% c(target_name,"uncertainty","date"))],
                                         train_set[,target_name],
                                         train_set[,"uncertainty"],
                                         test_set[,-which(names(test_set) %in% c(target_name,"uncertainty","date"))],
                                         n_neighbor,
                                         test_set[,target_name],0,5000)
    
    eknn_certain <- c(eknn_certain,f_eknn$ypred)
    nas <- sum(is.na(f_eknn$ypred))
    if(nas > 0){
      print("ISSUE EKNN")
      data <- NULL
    }
    
  }
  print("fini")
  return(eknn_certain)
}






#Experiment Pipeline
pipeline_regression <- function(file,n_neighbor,prediction_horizon,start_day,to_remove,
                                target,target_name,noise_level,imp,historic_length,historic_variables){
  #celui l?
  oid <- read.csv(file)
  
  smooth <- 7
  
  r <- 50
  
  clean_data <- oid[which(oid$location == "France"),]
  
  clean_data <- clean_data[-c(1:150),]
  
  clean_data <- clean_data[,c(4,6,9)]
  
  clean_data <- clean_data[complete.cases(clean_data),]
  
  clean_data[which(clean_data$new_deaths < 0),"new_deaths"] <- abs(clean_data[which(clean_data$new_deaths < 0),"new_deaths"])
  
  clean_data[which(clean_data$new_cases < 0),"new_cases"] <- abs(clean_data[which(clean_data$new_cases < 0),"new_cases"])
  
  clean_data[which(is.na(clean_data$new_deaths)),"new_deaths"] <- 0
  
  clean_data[which(is.na(clean_data$new_cases)),"new_cases"] <- 0
  
  
  
  
  result_matrix_temporal <- data.frame(matrix(nrow = 1,ncol = 4))
  colnames(result_matrix_temporal) <- c("noise_level","EKNNC","EKNN"
                                        ,"Naive")
  
  result_matrix_temporal[1,1] <- noise_level
  
  
  for (t in 1:r) {
    reference <- clean_data
    
    noised_data <- noise(clean_data,target,noise_level)$noised_data
    for (i in 1:5) {
      noised_data[i,target_name] <- reference[i,target_name] 
      
    }
    for (i in (nrow(noised_data)-4):nrow(noised_data)) {
      noised_data[i,target_name] <- reference[i,target_name] 
    }
    
    if(imp == "cma"){
      imputed_by_temporal <- weighted_evidential_imputation(noised_data,target,reference,0.05)$imputed_data
    }
    if(imp == "locf"){
      imputed_by_temporal <- mean_imputation(noised_data,target,reference,0.1)$imputed_data
    }
    if(noise_level != 0 & imp == "teknn"){imputed_by_temporal <- learning_imputation(noised_data,historic_length,historic_variables,clean_data,to_remove,target)$imputed_data 
    }else {  imputed_by_temporal <- clean_data
    imputed_by_temporal$uncertainty <- 1
    }
    
    original <- ncol(imputed_by_temporal)
    imputed_by_temporal$s_deaths <- NA 
    print("ok")
    
    for (i in smooth:nrow(imputed_by_temporal)) {
      deaths <- imputed_by_temporal[(i-(smooth-1)):i,target]
      imputed_by_temporal[i,(original+1)] <- median(deaths)
      
    }
    
    
    for (i in 1:nrow(imputed_by_temporal)) {
      if(is.na(imputed_by_temporal[i,(original+1)])){
        imputed_by_temporal[i,(original+1)] <- imputed_by_temporal[i,target]
      }  
    }
    
    imputed_by_temporal[,target] <- imputed_by_temporal[,(original+1)]
    imputed_by_temporal <- imputed_by_temporal[,1:original]
    imputed_by_temporal[,target] <- as.integer(imputed_by_temporal[,target])
    #print(sum(is.na(clean_data)))
    
    imputed_by_temporal$s_cases <- NA
    for (i in smooth:nrow(imputed_by_temporal)) {
      cases <- imputed_by_temporal[(i-(smooth-1)) :i,"new_cases"] 
      imputed_by_temporal[i,(original+1)] <- median(cases)
    }
    
    
    for (i in 1:nrow(imputed_by_temporal)) {
      if(is.na(imputed_by_temporal[i,(original+1)])){
        imputed_by_temporal[i,(original+1)] <- imputed_by_temporal[i,"new_cases"]
      }  
    }
    
    imputed_by_temporal[,"new_cases"] <- imputed_by_temporal[,(original+1)]
    imputed_by_temporal <- imputed_by_temporal[,1:original]
    imputed_by_temporal[,"new_cases"] <- as.integer(imputed_by_temporal[,"new_cases"])

    
    reference$uncertainty <- 1
    
    original <- ncol(reference)
    reference$s_deaths <- NA 
    
    
    for (i in smooth:nrow(reference)) {
      deaths <- reference[(i-(smooth-1)):i,target]
      reference[i,(original+1)] <- median(deaths)
      
    }
    
    
    for (i in 1:nrow(reference)) {
      if(is.na(reference[i,(original+1)])){
        reference[i,(original+1)] <-reference[i,target]
      }  
    }
    
    reference[,target] <- reference[,(original+1)]
    reference <- reference[,1:original]
    reference[,target] <- as.integer(reference[,target])
    #print(sum(is.na(clean_data)))
    
    reference$s_cases <- NA
    for (i in smooth:nrow(reference)) {
      cases <- reference[(i-(smooth-1)) :i,"new_cases"] 
      reference[i,(original+1)] <- median(cases)
    }
    
    
    for (i in 1:nrow(reference)) {
      if(is.na(reference[i,(original+1)])){
        reference[i,(original+1)] <- reference[i,"new_cases"]
      }  
    }
    
    reference[,"new_cases"] <- reference[,(original+1)]
    reference <- reference[,1:original]
    reference[,"new_cases"] <- as.integer(reference[,"new_cases"])
    
    
    ##
    
  
    
    #date
    
    imputed_by_temporal$date <- as.double(imputed_by_temporal$date)
    true_date <- reference$date
    reference$date <- as.double(reference$date)
    
    
    #noised_data$date <- as.double(noised_data$date)
    
    #KNN
    imputed_by_temporal <- history(imputed_by_temporal,historic_variables,historic_length,target,prediction_horizon)
    reference <- history(reference,historic_variables,historic_length,target,prediction_horizon)
    ###
    imputed_by_temporal <-  imputed_by_temporal[,-2]
    reference <-  reference[,-2]
    #print(colnames(imputed_by_temporal))
    #print(colnames(reference))
    #max normalization
    for (i in 4:ncol(imputed_by_temporal)) {
      imputed_by_temporal[,i] <- imputed_by_temporal[,i]/max(imputed_by_temporal[,i])
    }
    
    for (i in 4:ncol(reference)) {
      reference[,i] <- reference[,i]/max(reference[,i])
    }
    
    
    print("ici")

    reality <- c()
    for (i in start_day:(nrow(reference)-prediction_horizon)) {
      
      reality <- c(reality,reference[i+prediction_horizon,target_name])
      date <- c(date,as.character(reference[i+prediction_horizon,"date"]))
      
    }
    
    print("ici")
    if(t==1){imputres <- imputed_by_temporal$new_deaths }else{
      deaths_result <- imputed_by_temporal$new_deaths
      imputres <- imputres + imputed_by_temporal$new_deaths 
    }
    print("ici")
    if(t==1){
      naive_result_imputation_by_temporal <- evaluation_naive(imputed_by_temporal[,-to_remove],target_name,reference[,-to_remove],start_day,prediction_horizon,n_neighbor)
      result_matrix_temporal[1,4] <- sqrt(median((naive_result_imputation_by_temporal - reality)^2))
    }
    
    else{
      naive_result <- evaluation_naive(imputed_by_temporal[,-to_remove],target_name,reference[,-to_remove],start_day,prediction_horizon,n_neighbor)
      result_matrix_temporal[1,4] <- result_matrix_temporal[1,4] + sqrt(median((naive_result - reality)^2))
      naive_result_imputation_by_temporal <- naive_result + naive_result_imputation_by_temporal
    }
    
    print("ici")
    if(t==1){
      eknn_result_imputation_by_temporal <- evaluation_eknn(imputed_by_temporal,target_name,to_remove,reference,start_day,prediction_horizon,n_neighbor)
      result_matrix_temporal[1,3] <- sqrt(median((eknn_result_imputation_by_temporal - reality)^2))
    }
  
    else{
      eknn_result <-  evaluation_eknn(imputed_by_temporal,target_name,to_remove,reference,start_day,prediction_horizon,n_neighbor)
      result_matrix_temporal[1,3] <- result_matrix_temporal[1,3]  + sqrt(median((eknn_result - reality)^2))
      #print(eknn_result_imputation_by_temporal)
      eknn_result_imputation_by_temporal <- eknn_result + eknn_result_imputation_by_temporal
    }
    print("ici")
    if(t ==1){
      eknnc_result_imputation_by_temporal <- evaluation_eknnc(imputed_by_temporal,target_name,to_remove,reference,start_day,prediction_horizon,n_neighbor)
      result_matrix_temporal[1,2] <- sqrt(median((eknnc_result_imputation_by_temporal - reality)^2))
    }
    
    else{
      eknnc_result <-  evaluation_eknnc(imputed_by_temporal,target_name,to_remove,reference,start_day,prediction_horizon,n_neighbor)
      result_matrix_temporal[1,2] <- result_matrix_temporal[1,2] + sqrt(median((eknnc_result - reality)^2))
      eknnc_result_imputation_by_temporal <- eknnc_result + eknnc_result_imputation_by_temporal
    }
    
    
  }
  print("ici")
  start <- length(true_date) - dim(reference)[1]
  start <- start + 1
  true_date <- true_date[start:length(true_date)]
  reference$date <- true_date
  print("MAX")
  print(max(imputed_by_temporal$uncertainty))
  print("MEAN")
  print(mean(imputed_by_temporal$uncertainty))
  
  date <- c()
  reality <- c()
  
  for (i in start_day:(nrow(reference)-prediction_horizon)) {
    
    reality <- c(reality,reference[i+prediction_horizon,target_name])
    date <- c(date,as.character(reference[i+prediction_horizon,"date"]))
    
  }
  
  
  print("ici ?")
  
  #print(as.integer(eknn_result_imputation_by_temporal))
  imputres <- imputres/r
  eknnc_result_imputation_by_temporal <- eknnc_result_imputation_by_temporal/r
  naive_result_imputation_by_temporal <- naive_result_imputation_by_temporal/r
  eknn_result_imputation_by_temporal <- eknn_result_imputation_by_temporal/r
  #dff_eknn <-  reality - eknn_result_imputation_by_temporal
  #dff_eknnc <-   reality - eknnc_result_imputation_by_temporal 
  #df_diff <- cbind(date,dff_eknn)
  #df_diff <- cbind(df_diff,dff_eknnc)
  #df_diff <- as.data.frame(df_diff)
  #colnames(df_diff) <- c("date","eknn","eknnc")
  
  #titre <- "Difference between predictions and reality : temporal : K = :"
  #titre <- "Difference between predictions and reality : mean : K = :"
  i_df <- cbind(imputed_by_temporal$date,imputres)
  i_df <- as.data.frame(i_df)
  i_df <- cbind(i_df,reference$new_deaths)
  colnames(i_df) <- c("date","deaths","real")
  
  #ggplot(i_df, aes(date)) + 
  #  geom_line(aes(y = deaths, colour = "imputed deaths" )) +
  #  geom_line(aes(y = real, colour = "real data" )) +
  #  ggtitle(paste("imputation_temp",paste(n_neighbor,paste(": history =",historic_length)))) +
  #  theme(panel.background = element_rect(fill = 'white'),
  #        plot.background=element_rect(fill = "white"),
  #        panel.background = element_rect(fill = 'gray20'),
  #        plot.background=element_rect(fill = "gray20"),legend.position="bottom",plot.title = element_text(size = 20, face = "bold"),
  #        axis.text = element_text(size = 16,face = "bold",colour = "white")
  #        ,axis.title=element_text(size=16,face="bold",colour = "white"),
  #        text = element_text(size=16,face = "bold",colour = "white"),legend.direction = "horizontal",
  #        legend.background = element_rect(fill = "gray20", color = NA),
  #        legend.key = element_rect(color = "gray20", fill = "gray20"),
  #        legend.title = element_text(color = "white"),
  #        legend.text = element_text(color = "white",size = 16),
  #        legend.key.size = unit(2.5,"line"),
  #        panel.grid.major = element_blank(), 
  #        panel.grid.minor = element_blank(),
  #        panel.border = element_rect(color = "white", fill = NA )) + 
  #  xlab("Date") + ylab("Prediction") 
  #ggsave(paste("imp-temporal",paste(paste(paste(n_neighbor,noise_level,sep = " "),historic_length,sep = " "),".png",sep = ""),sep = ""),device = "png",height = 7,width = 10)
  
  
  #print("MOYENNES")
  #print("KNN")
  #print("EKNN")
  #print(as.integer(eknn_result_imputation_by_temporal))
  #print(sum(is.na(clean_data)))
  #print(length(eknn_result_imputation_by_temporal))
  #print(length(knn_result_imputation_by_temporal))
  #print(length(date))
  #print(length(reality))
  print("RMSE")
  print(sqrt(mean((eknn_result_imputation_by_temporal - reality)^2)))
  print(sqrt(mean((eknnc_result_imputation_by_temporal - reality)^2)))
  print(sqrt(mean((naive_result_imputation_by_temporal - reality)^2)))
  
  result_matrix_temporal[1,3] <- result_matrix_temporal[1,3]/r
  
  result_matrix_temporal[1,2] <- result_matrix_temporal[1,2]/r
  
  result_matrix_temporal[1,4] <- result_matrix_temporal[1,4]/r
  
  
  df_temporal <- cbind(cbind(eknn_result_imputation_by_temporal,naive_result_imputation_by_temporal),
                       eknnc_result_imputation_by_temporal)
  
  
  
  df_temporal <- as.data.frame(df_temporal)
  df_temporal <- cbind(date,df_temporal)
  df_temporal <- cbind(df_temporal,reality)
  colnames(df_temporal) <- c("date","eknn","naive","eknnc","reality")
  if(imp == "locf"){
    titre <- "LOCF imputation : K = :"
  }
  if(imp == "cma"){
    titre <- "CMA imputation : K = :"
  }
  if(imp == "teknn"){
    titre <- "Learning imputation : K = :"
  }
  df_temporal$date <- as.Date(df_temporal$date)
  
  
  ggplot(df_temporal, aes(date)) + 
    geom_line(aes(y = eknnc, colour = "EKNN" ),size = 1) +
    geom_line(aes(y = eknn, colour = "EKNN uncertain labels" ),size = 1) + 
    geom_line(aes(y = naive, colour = "Baseline" ),size = 1) +
    geom_line(aes(y = reality, colour = "Reality" ),size = 1) +
    scale_x_date(date_labels = "%b %Y") + theme_minimal() +
    theme(panel.background = element_rect(fill = 'white'),
          plot.background=element_rect(fill = "white"),
          legend.position="top",plot.title = element_text(size = 20, face = "bold"),
          axis.text = element_text(size = 13)
          ,axis.title=element_text(size=20,face="bold"),
          legend.direction = "horizontal",
          legend.text = element_text(size = 20),
          legend.key.size = unit(1.5,"line"),
          legend.title = element_blank(), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) + 
    xlab("Date") + ylab("Number of deaths") 
  ggsave(paste(paste("pipe-",imp,sep = ""),paste(paste(paste(n_neighbor,noise_level,sep = " "),historic_length,sep = " "),".png",sep = ""),sep = ""),device = "png",height = 7,width = 10)
  write.csv(df_temporal, file = paste(paste("pipe-",imp,sep = ""),paste(paste(paste(n_neighbor,noise_level,sep = " "),historic_length,sep = " "),".csv",sep = ""),sep = ""))
  
  return(list(temporal_imputation = result_matrix_temporal))  
  
}











#Imputation evaluation
pipeline_imputation <- function(file,n_neighbor,prediction_horizon,start_day,to_remove,
                                target,target_name,noise_level,historic_length,historic_variables){
  oid <- read.csv(file)
  
  smooth <- 7
  
  r <- 50
  
  clean_data <- oid[which(oid$location == "France"),]
  
  clean_data <- clean_data[-c(1:150),]
  
  clean_data <- clean_data[,c(4,6,9)]
  
  clean_data <- clean_data[complete.cases(clean_data),]
  
  clean_data[which(clean_data$new_deaths < 0),"new_deaths"] <- abs(clean_data[which(clean_data$new_deaths < 0),"new_deaths"])
  
  clean_data[which(clean_data$new_cases < 0),"new_cases"] <- abs(clean_data[which(clean_data$new_cases < 0),"new_cases"])
  
  clean_data[which(is.na(clean_data$new_deaths)),"new_deaths"] <- 0
  
  clean_data[which(is.na(clean_data$new_cases)),"new_cases"] <- 0
  
  #original <- ncol(clean_data)
  #clean_data$s_deaths <- NA 
  
  
  #for (i in smooth:nrow(clean_data)) {
  #  deaths <- clean_data[(i-(smooth-1)):i,target]
  #  clean_data[i,4] <- median(deaths)
  #  
  #}
  
  
  #for (i in 1:nrow(clean_data)) {
  #  if(is.na(clean_data[i,4])){
  #    clean_data[i,4] <- clean_data[i,3]
  #  }  
  #}
  
  #clean_data[,target] <- clean_data[,(original+1)]
  #clean_data <- clean_data[,1:original]
  #clean_data[,target] <- as.integer(clean_data[,target])
  #print(sum(is.na(clean_data)))
  
  #clean_data$s_cases <- NA
  #for (i in smooth:nrow(clean_data)) {
  #  cases <- clean_data[(i-(smooth-1)) :i,"new_cases"] 
  #  clean_data[i,4] <- median(cases)
  #}
  
  
  #for (i in 1:nrow(clean_data)) {
  #  if(is.na(clean_data[i,4])){
  #    clean_data[i,4] <- clean_data[i,"new_cases"]
  #  }  
  #}
  
  #clean_data[,"new_cases"] <- clean_data[,(original+1)]
  #clean_data <- clean_data[,1:original]
  #clean_data[,"new_cases"] <- as.integer(clean_data[,"new_cases"])
  
  
  
  for (t in 1:r) {
    reference <- clean_data
    
    noised_data <- noise(clean_data,target,noise_level)$noised_data
    
    for (i in 1:1) {
      noised_data[i,target_name] <- reference[i,target_name] 
      
    }
    for (i in nrow(noised_data):nrow(noised_data)) {
      noised_data[i,target_name] <- reference[i,target_name] 
    }
    
    
    imputed_by_temporal <- weighted_evidential_imputation(noised_data,target,reference,0.05)$imputed_data
    imputed_by_mean <- mean_imputation(noised_data,target,reference,0.05)$imputed_data
    if(noise_level != 0){  imputed_by_learning <- learning_imputation(noised_data,historic_length,historic_variables,clean_data,to_remove,target)$imputed_data 
    }else{  imputed_by_learning <- clean_data
    imputed_by_learning$uncertainty <- 1
    }
    
    
    #cma
    original <- ncol(imputed_by_temporal)
    imputed_by_temporal$s_deaths <- NA 
    
    
    for (i in smooth:nrow(imputed_by_temporal)) {
      deaths <- imputed_by_temporal[(i-(smooth-1)):i,target]
      imputed_by_temporal[i,(original+1)] <- median(deaths)
      
    }
    
    
    for (i in 1:nrow(imputed_by_temporal)) {
      if(is.na(imputed_by_temporal[i,(original+1)])){
        imputed_by_temporal[i,(original+1)] <- imputed_by_temporal[i,target]
      }  
    }
    
    imputed_by_temporal[,target] <- imputed_by_temporal[,(original+1)]
    imputed_by_temporal <- imputed_by_temporal[,1:original]
    imputed_by_temporal[,target] <- as.integer(imputed_by_temporal[,target])
    #print(sum(is.na(clean_data)))
    
    imputed_by_temporal$s_cases <- NA
    for (i in smooth:nrow(imputed_by_temporal)) {
      cases <- imputed_by_temporal[(i-(smooth-1)) :i,"new_cases"] 
      imputed_by_temporal[i,(original+1)] <- median(cases)
    }
    
    
    for (i in 1:nrow(imputed_by_temporal)) {
      if(is.na(imputed_by_temporal[i,(original+1)])){
        imputed_by_temporal[i,(original+1)] <- imputed_by_temporal[i,"new_cases"]
      }  
    }
    
    imputed_by_temporal[,"new_cases"] <- imputed_by_temporal[,(original+1)]
    imputed_by_temporal <- imputed_by_temporal[,1:original]
    imputed_by_temporal[,"new_cases"] <- as.integer(imputed_by_temporal[,"new_cases"])
    
    #locf
    original <- ncol(imputed_by_mean)
    imputed_by_mean$s_deaths <- NA 
    
    
    for (i in smooth:nrow(imputed_by_mean)) {
      deaths <- imputed_by_mean[(i-(smooth-1)):i,target]
      imputed_by_mean[i,(original+1)] <- median(deaths)
      
    }
    
    
    for (i in 1:nrow(imputed_by_mean)) {
      if(is.na(imputed_by_mean[i,(original+1)])){
        imputed_by_mean[i,(original+1)] <- imputed_by_mean[i,target]
      }  
    }
    
    imputed_by_mean[,target] <- imputed_by_mean[,(original+1)]
    imputed_by_mean <- imputed_by_mean[,1:original]
    imputed_by_mean[,target] <- as.integer(imputed_by_mean[,target])
    #print(sum(is.na(clean_data)))
    
    imputed_by_mean$s_cases <- NA
    for (i in smooth:nrow(imputed_by_mean)) {
      cases <- imputed_by_mean[(i-(smooth-1)) :i,"new_cases"] 
      imputed_by_mean[i,(original+1)] <- median(cases)
    }
    
    
    for (i in 1:nrow(imputed_by_mean)) {
      if(is.na(imputed_by_mean[i,(original+1)])){
        imputed_by_mean[i,(original+1)] <- imputed_by_mean[i,"new_cases"]
      }  
    }
    
    imputed_by_mean[,"new_cases"] <- imputed_by_mean[,(original+1)]
    imputed_by_mean <- imputed_by_mean[,1:original]
    imputed_by_mean[,"new_cases"] <- as.integer(imputed_by_mean[,"new_cases"])
    

    #Teknn
    original <- ncol(imputed_by_learning)
    imputed_by_learning$s_deaths <- NA 
    
    
    for (i in smooth:nrow(imputed_by_learning)) {
      deaths <- imputed_by_learning[(i-(smooth-1)):i,target]
      imputed_by_learning[i,(original+1)] <- median(deaths)
      
    }
    
    
    for (i in 1:nrow(imputed_by_learning)) {
      if(is.na(imputed_by_learning[i,(original+1)])){
        imputed_by_learning[i,(original+1)] <- imputed_by_learning[i,target]
      }  
    }
    
    imputed_by_learning[,target] <- imputed_by_learning[,(original+1)]
    imputed_by_learning <- imputed_by_learning[,1:original]
    imputed_by_learning[,target] <- as.integer(imputed_by_learning[,target])
    #print(sum(is.na(clean_data)))
    
    imputed_by_learning$s_cases <- NA
    for (i in smooth:nrow(imputed_by_learning)) {
      cases <- imputed_by_learning[(i-(smooth-1)) :i,"new_cases"] 
      imputed_by_learning[i,(original+1)] <- median(cases)
    }
    
    
    for (i in 1:nrow(imputed_by_learning)) {
      if(is.na(imputed_by_learning[i,(original+1)])){
        imputed_by_learning[i,(original+1)] <- imputed_by_learning[i,"new_cases"]
      }  
    }
    
    imputed_by_learning[,"new_cases"] <- imputed_by_learning[,(original+1)]
    imputed_by_learning <- imputed_by_learning[,1:original]
    imputed_by_learning[,"new_cases"] <- as.integer(imputed_by_learning[,"new_cases"])
    
    
    #ref
    original <- ncol(reference)
    reference$s_deaths <- NA 
    
    
    for (i in smooth:nrow(reference)) {
      deaths <- reference[(i-(smooth-1)):i,target]
      reference[i,(original+1)] <- median(deaths)
      
    }
    
    
    for (i in 1:nrow(reference)) {
      if(is.na(reference[i,(original+1)])){
        reference[i,(original+1)] <- reference[i,target]
      }  
    }
    
    reference[,target] <- reference[,(original+1)]
    reference <- reference[,1:original]
    reference[,target] <- as.integer(reference[,target])
    #print(sum(is.na(clean_data)))
    
    reference$s_cases <- NA
    for (i in smooth:nrow(reference)) {
      cases <- reference[(i-(smooth-1)) :i,"new_cases"] 
      reference[i,(original+1)] <- median(cases)
    }
    
    
    for (i in 1:nrow(reference)) {
      if(is.na(reference[i,(original+1)])){
        reference[i,(original+1)] <- reference[i,"new_cases"]
      }  
    }
    
    reference[,"new_cases"] <- reference[,(original+1)]
    reference <- reference[,1:original]
    reference[,"new_cases"] <- as.integer(reference[,"new_cases"])
    
    
    ##
    
    
    #date
    
    
    #noised_data$date <- as.double(noised_data$date)
    
    #KNN


    
    
    
    #temporal
    #if(t==1){knn_result_imputation_by_temporal <- evaluation_wknn(imputed_by_temporal[,-to_remove],target_name,reference[,-to_remove],start_day,prediction_horizon,5)}else{
    #  knn_result <- evaluation_wknn(imputed_by_temporal[,-to_remove],target_name,reference[,-to_remove],start_day,prediction_horizon,5)
    #  knn_result_imputation_by_temporal <- knn_result + knn_result_imputation_by_temporal
    #}
    
    if(t==1)
      {imputres <-  imputed_by_temporal$new_deaths
       imputresm <- imputed_by_mean$new_deaths
       imputresl <- imputed_by_learning$new_deaths}else{
       deaths_result <- imputed_by_temporal$new_deaths
       imputres <- imputres + imputed_by_temporal$new_deaths
       deaths_resultm <- imputed_by_mean$new_deaths
       imputresm <- imputresm + imputed_by_mean$new_deaths
       deaths_resultl <- imputed_by_learning$new_deaths
       imputresl <- imputresl + imputed_by_learning$new_deaths 
    }
    
    

    
  }

  #print(as.integer(eknn_result_imputation_by_temporal))
  imputres <- imputres/r
  imputresl <- imputresl/r 
  imputresm <- imputresm/r 
  #dff_eknn <-  reality - eknn_result_imputation_by_temporal
  #dff_eknnc <-   reality - eknnc_result_imputation_by_temporal 
  #df_diff <- cbind(date,dff_eknn)
  #df_diff <- cbind(df_diff,dff_eknnc)
  #df_diff <- as.data.frame(df_diff)
  #colnames(df_diff) <- c("date","eknn","eknnc")
  
  i_df <- cbind(imputres,imputresm)
  i_df <- cbind(i_df,imputresl)
  i_df <- cbind(i_df,reference$new_deaths)
  i_df <- as.data.frame(i_df)
  i_df <- cbind(noised_data$date,i_df)
  colnames(i_df) <- c("date","temporal","mean","learning","real")
  i_df$date <- as.Date(i_df$date)
  ggplot(i_df, aes(x = date)) + 
    geom_line(aes(y = temporal, colour = "CMA" ),size = 1) +
    geom_line(aes(y = mean, colour = "LOCF" ),size = 1) +
    geom_line(aes(y = learning, colour = "TEKNN" ),size = 1) +
    geom_line(aes(y = real, colour = "Reality" ),size = 1) +
    scale_x_date(date_labels = "%b %Y") + theme_minimal() +
    theme(legend.position = "top",
          legend.direction = "horizontal",
          legend.text = element_text(size = 20),
      legend.title = element_blank(),
      axis.text.x = element_text(size = 13),
      axis.title.x = element_text(size = 30),
      axis.title.y = element_text(size = 30)) + 
    xlab("Date") + ylab("Number of deaths") 
  ggsave(paste("imp",paste(noise_level,".png",sep = ""),sep = ""),device = "png",height = 7,width = 10)
  write.csv(i_df, file = paste("imput",paste(noise_level,".csv",sep = "")))
  

  
  return(list(temporal_imputation = 0)) 
  }
  
