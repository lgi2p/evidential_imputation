#setwd("Projet-COVID")
setwd("/Users/nsc/Downloads/Imputation_credibiliste-main/Evidential_imputation")
source("used functions.R")

set.seed(123)
n_neighbor <- c(1,10,20)
prediction_horizon <- 7
start_day <- 21
to_remove <- 1
target <- 3
target_name <- "new_deaths"
noise_level <- c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7)
historic_length <- 4
historic_variables <- c(2,3)
imp <- "cma"
file <- "C:/Users/Rayan/Documents/Projet-COVID/owid-covid-data.csv"
file <- "C:/Users/Rayan/Documents/COVIDREZ/OLD 2/bilan/owid-covid-data.csv"
file <- "owid-covid-data.csv"


result_matrix_temporal <- data.frame(matrix(nrow = 0,ncol = 4))
colnames(result_matrix_temporal) <- c("noise_level","EKNNC","EKNN","Naive")


for(l in  historic_length){
  path <- getwd()
  dir.create(as.character(l))
  setwd(as.character(l))
  print(l)
  for (j in n_neighbor) {
    print(j)
    for (i in noise_level) {
      print(i)
      result <- pipeline_regression(file,j,prediction_horizon,
                                    start_day,to_remove,target,target_name,i,imp,l,historic_variables)
      result_matrix_temporal <- rbind(result_matrix_temporal,result$temporal_imputation)
      }
    
    ggplot(result_matrix_temporal, aes(noise_level)) + 
      geom_line(aes(y = EKNNC ,colour = "EKNN" ),size = 2) +
      geom_line(aes(y = EKNN ,colour = "EKNN uncertain labels" ),size = 2) +
      #geom_line(aes(y = Naive, colour = "Baseline" ),size = 2) +
      theme(panel.background = element_rect(fill = 'white'),
            plot.background=element_rect(fill = "white"),
            legend.position="top",plot.title = element_text(size = 20, face = "bold"),
            axis.text = element_text(size = 20,face = "bold"),
            axis.title=element_text(size=20,face="bold"),
            legend.direction = "horizontal",
            legend.text = element_text(size = 20),
            legend.key.size = unit(1.5,"line"),
            legend.title = element_blank(), 
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank()) + 
     xlab("Noise level") + ylab("RMedSE") 
    ggsave(paste(imp,paste(paste(j,l,sep = " "),".png",sep = ""),sep = ""),height = 7,width = 10)
    write.csv(result_matrix_temporal, file = paste(imp,paste(paste(j,l,sep = " "),".csv",sep = ""),sep = "") )
    result_matrix_temporal <- data.frame(matrix(nrow = 0,ncol = 4))
    colnames(result_matrix_temporal) <- c("noise_level","EKNNC","EKNN","Naive")
    
    
  }
  setwd(path)
}





