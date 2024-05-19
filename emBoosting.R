# I will be using EM to soft cluster, and provide hard cluster initial values
# priming the world
 
# libraries
library(janitor)
library(dplyr)
library(MASS)
library(cluster)
library(dbscan)
library(mixtools)
library(factoextra)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(insight)

set.seed(1971)

seed <- 1971

# function for fitting hard cluster to EM algorithm
clus_init_EM <- function(clus_method, cluster_size, dataset) {
  
  set.seed(seed)
  
  df_temp <- dataset
  
  if (clus_method == "kmeans") {
    mvn_cluster <- kmeans(df_temp, centers = cluster_size)
    print(min(table(mvn_cluster$cluster)))
    
    if (!is.infinite(min(table(mvn_cluster$cluster))) & min(table(mvn_cluster$cluster)) > 40) {
      df_temp$cluster <- mvn_cluster$cluster
    } else {
      print("Clusters had an infinite subclass")
      return(NULL)
    }
    
  } else if (clus_method == "hmedian") {
    mvn_cluster <- pam(df_temp, cluster_size, metric = "euclidean", stand = FALSE)
    
    print(min(table(mvn_cluster$cluster)))
    
    if (!is.infinite(min(table(mvn_cluster$cluster))) & min(table(mvn_cluster$cluster)) > 40) {
      df_temp$cluster <- mvn_cluster$cluster
    } else {
      print("Clusters had Infinite subclass")
      return(NULL)
    }
    
  } else if (clus_method == "hclust") {
    mvn_cluster <- hclust(dist(df_temp))
    
    print(min(table(cutree(tree = mvn_cluster, k = cluster_size))))
    
    if (!is.infinite(min(table(cutree(tree = mvn_cluster, k = cluster_size)))) & min(table(cutree(tree = mvn_cluster, k = cluster_size))) > 40) {
      df_temp$cluster <- cutree(tree = mvn_cluster, k = cluster_size)
    } else {
      print("Clusters had Infinite subclass")
      return(NULL)
    } 
  } else if (clus_method == "dbscan") {
      
    mvn_cluster <- dbscan(df_temp, eps = 5)
    
    print(min(table(mvn_cluster$cluster)))
    
    if (!is.infinite(min(table(mvn_cluster$cluster))) & min(table(mvn_cluster$cluster)) > 40) {
        df_temp$cluster <- mvn_cluster$cluster
    } else {
      print("Clusters had Infinite subclass")
      return(NULL)
    }
  }
  
  # gets the estimated mean
  mean_list <- by(df_temp[,-length(names(df_temp))], df_temp$cluster, colMeans)
  
  sigma_list <- by(df_temp[,-length(names(df_temp))], df_temp$cluster, cov)
  
  cluster_lambda <- table(df_temp$cluster) / length(df_temp$cluster)
  
  # applies EM for mixed normal
  mvn_EM <- mvnormalmixEM(df_temp[,-length(names(df_temp))], arbmean = TRUE,
                          mu = mean_list, sigma = sigma_list, lambda = cluster_lambda, epsilon = 1e-02)

  return(mvn_EM)
}

## actual
df_nba <- read.table("Seasons_Stats.csv", sep = ",", header = TRUE)

df_nba_player <- df_nba %>%
  dplyr::filter(Year %in% c("2016", "2017"))

df_nba_player_filter <- df_nba_player %>%
  dplyr::select(Player, Pos, X3P., X2P., TRB., AST., STL., BLK., TOV.)

df_nba_player_filter <- df_nba_player_filter %>%
  drop_na()

df_nba_player_filter$Pos <- ifelse(df_nba_player_filter$Pos == "PF-C", "PF", df_nba_player_filter$Pos)

gghistogram(data = df_nba_player_filter, x = "X3P.", color = "Pos", bins = 60, xlab = "3-Point %", ylab = "Count", title = "Histogram of 3-Point Accuracy")

df_nba_player_filter %>%
  dplyr::group_by(Pos) %>%
  summarise(`3P%` = mean(X3P.),
            `2P%` = mean(X2P.),
            TRB = mean(TRB.),
            AST = mean(AST.),
            STL = mean(STL.),
            BLK = mean(BLK.),
            TOV = mean(TOV.)) %>%
  export_table(df, format = "html")

df_nba_numeric <- df_nba_player_filter %>%
  dplyr::select_if(is.numeric)

mvn_clusters_nba <- list()

for (i in cluster_methods) {
  
  print(i)
  
  mvn_clusters_temp <- clus_init_EM(i, 5, df_nba_numeric)
  
  mvn_clusters_nba[[i]] <- mvn_clusters_temp
  
}

# creates random values for random EM
# lambda
lambda_base <- runif(5)
lambda_base <- lambda_base / sum(lambda_base)

# mu
mu_base <- matrix(runif(5*5, 0, 21), nrow = 5, ncol = 5)

# sigma
sigma_base <- list(matrix(runif(25, 0, 21), nrow = 5, ncol = 5, byrow = TRUE),
              matrix(runif(25, 0, 21), nrow = 5, ncol = 5, byrow = TRUE),
              matrix(runif(25, 0, 21), nrow = 5, ncol = 5, byrow = TRUE),
              matrix(runif(25, 0, 21), nrow = 5, ncol = 5, byrow = TRUE), 
              matrix(runif(25, 0, 21), nrow = 5, ncol = 5, byrow = TRUE))

mvn_base <- mvnormalmixEM(df_nba_numeric, lambda_base, k = 5, epsilon = 1e-02)

v_mean_act <- colMeans(df_nba_numeric, na.rm = TRUE)

v_actual_lambda <- table(df_nba_player_filter$Pos)
v_actual_lambda <- v_actual_lambda / sum(v_actual_lambda)
v_actual_lambda <- sort(v_actual_lambda, decreasing = TRUE)

df_actual_mean <- df_nba_player_filter %>%
    dplyr::group_by(Pos) %>%
    summarise(
      X3P. = mean(X3P.),
      X2P. = mean(X2P.),
      TRB. = mean(TRB.),
      AST. = mean(AST.),
      STL. = mean(STL.),
      BLK. = mean(BLK.),
      TOV. = mean(TOV.))

df_actual_var <- by(df_nba_player_filter[-c(1,2)], df_nba_player_filter$Pos,cov)

# making a base df
df_lambda <- as.data.frame(as.list(v_actual_lambda))

df_lambda$Method <- "actual"

# completely random EM estimates

df_base_lamdba <- sort(mvn_base$lambda, decreasing = TRUE)

df_base_lamdba <- append(df_base_lamdba, "random")

names(df_base_lamdba) <- c("SF", "PG", "SG", "PF", "C", "Method")

df_lambda <- rbind(df_lambda, df_base_lamdba)

# cluster methods
for (i in names(mvn_clusters_nba)) {
  em_temp <- mvn_clusters_nba[[i]]
  
  v_lambda_temp <- sort(em_temp$lambda, decreasing = TRUE)
  
  v_lambda_temp <- append(v_lambda_temp, i)
  
  names(v_lambda_temp) <- c("SF", "PG", "SG", "PF", "C", "Method")
  
  df_lambda <- rbind(df_lambda, v_lambda_temp)
}

library(ggthemes)
library(reshape2)

df_reshape <- df_lambda
df_reshape[-length(names(df_reshape))] <- sapply(df_reshape[-length(names(df_reshape))], as.numeric)

v_actual <- df_reshape[df_reshape$Method == "actual", -length(names(df_reshape))]

df_reshape <- melt(df_lambda, measure.vars = c("SF", "PG", "SG", "PF", "C"))
df_reshape$value <- as.numeric(df_reshape$value)

df_reshape <- df_reshape %>%
  dplyr::filter(Method == "actual") %>%
  dplyr::select(variable, value) %>%
  dplyr::rename(actual = value) %>%
  dplyr::right_join(df_reshape, by = "variable") %>%
  dplyr::mutate(Error = abs(value - actual))

df_reshape$variable <- as.character(df_reshape$variable)

df_reshape <- df_reshape[order(df_reshape$variable),]

# plot
ggplot(df_reshape, aes(variable, round(value, digits = 3), class = Method, color = Method, size = Error^1.5)) +
  geom_point(alpha = 0.80) + 
  labs(title= "Accuracy of Lambda", 
       subtitle="By each Priming Method",
       x= "Player Position",
       y= "Lambda Estimates") +
  guides(size = guide_legend(title="Absolute Error"))

# actual mean
df_actual_mean <- df_nba_player_filter %>%
  dplyr::group_by(Pos) %>%
  summarise(
    X3P. = mean(X3P.),
    X2P. = mean(X2P.),
    TRB. = mean(TRB.),
    AST. = mean(AST.),
    STL. = mean(STL.),
    BLK. = mean(BLK.),
    TOV. = mean(TOV.))

column_names <- df_actual_mean$Pos

df_actual_mean <- as.data.frame(t(df_actual_mean[,names(df_actual_mean)[-1]]))

names(df_actual_mean) <- column_names

df_actual_mean <- df_actual_mean[,names(v_actual)]

df_actual_mean$Method <- "actual"

# completely random EM estimates mean
df_base_mean <- mvn_base$mu

df_base_mean <- data.frame(df_base_mean)

names(df_base_mean) <- mvn_base$lambda

df_base_mean <- df_base_mean[,as.character(sort(as.numeric(names(df_base_mean)), decreasing = TRUE))]

df_base_mean$Method <- "random"

names(df_base_mean) <- c("SF", "PG", "SG", "PF", "C", "Method")

df_mean <- rbind(df_actual_mean, df_base_mean)

# cluster methods
for (i in names(mvn_clusters_nba)) {
  
  em_temp <- mvn_clusters_nba[[i]]
  
  df_mean_temp <- data.frame(em_temp$mu)
  
  names(df_mean_temp) <- em_temp$lambda
  
  df_mean_temp <- df_mean_temp[,as.character(sort(as.numeric(names(df_mean_temp)), decreasing = TRUE))]
  
  df_mean_temp$Method <- i
  
  names(df_mean_temp) <- c("SF", "PG", "SG", "PF", "C", "Method")
  
  df_mean <- rbind(df_mean, df_mean_temp)
}

# hard coded... Sorry
df_mean$metrics <- rep(rownames(df_mean)[1:7], 4)

df_mean_actual <- df_mean %>%
  melt() %>%
  dplyr::filter(Method == "actual")

df_mean_actual <- df_mean_actual[,-1]

names(df_mean_actual) <- c("metrics", "variable", "actual value")

df_mean <- df_mean %>%
  melt()

df_mean_join <- df_mean %>% 
  dplyr::inner_join(df_mean_actual, by = c("metrics", "variable"))

df_mean_join$Error <- abs(df_mean_join$value - df_mean_join$`actual value`)

df_mean_join <- df_mean_join[,-5]

df_mean_plot <- df_mean_join %>%
  dplyr::group_by(Method, variable) %>%
  summarise(value = mean(value),
            MAE = mean(Error))

df_mean_plot$variable <- as.character(df_mean_plot$variable)

df_mean_plot <- df_mean_plot[order(df_mean_plot$variable),]

# plot
ggplot(df_mean_plot, aes(variable, round(value, digits = 3), class = Method, color = Method, size = MAE^1.5)) +
  geom_point(alpha = 0.80) + 
  labs(title= "Accuracy of Means", 
       subtitle="For each Priming Method",
       x= "Player Position",
       y= "Avg. Mean Estimates") +
  guides(size = guide_legend(title="Mean Absolute Error"))

# covariance
df_position <- melt(v_actual_lambda)

names(df_position) <- c("Position", "lambda")

df_position <- df_position$Position

# actual cov
cov.list <- lapply(unique(df_nba_player_filter$Pos),
                   function(x) cov(df_nba_player_filter[df_nba_player_filter$Pos==x,-c(1,2)],
                                  use="na.or.complete"))

position <- data.frame(Pos = unique(df_nba_player_filter$Pos), L1 = c(1,2,3,4,5))

df_actual_cov <- cov.list %>%
  melt() %>%
  as.data.frame() %>%
  dplyr::left_join(position, by = "L1") %>%
  dplyr::group_by(Pos) %>%
  dplyr::summarise(value = mean(value))

column_names <- df_actual_cov$Pos

df_actual_cov$Method <- "actual"

# completely random EM estimates mean
df_base_sigma <- mvn_base$sigma %>%
  melt() %>%
  group_by(L1) %>%
  summarise(value = mean(value))

position <- data.frame(Lambda = mvn_base$lambda, L1 = c(1,2,3,4,5))

position <- position[order(position$Lambda),]

position$Pos <- df_position

position <- position[,c(2,3)]

df_base_sigma <- df_base_sigma %>%
  dplyr::left_join(position, by = "L1")

df_base_sigma <- df_base_sigma[, -1]

df_base_sigma$Method <- "random"

df_sigma <- rbind(df_actual_cov, df_base_sigma)

# cluster methods
for (i in names(mvn_clusters_nba)) {
  
  em_temp <- mvn_clusters_nba[[i]]
  
  df_sigma_temp <- em_temp$sigma %>%
    melt() %>%
    group_by(L1) %>%
    summarise(value = mean(value))
  
  position <- data.frame(Lambda = mvn_clusters_nba[[i]]$lambda, L1 = c(1,2,3,4,5))
  
  position <- position[order(position$Lambda),]
  
  position$Pos <- df_position
  
  position <- position[,c(2,3)]
  
  df_sigma_temp <- df_sigma_temp %>%
    dplyr::left_join(position, by = "L1")
  
  df_sigma_temp <- df_sigma_temp[, -1]
  
  df_sigma_temp$Method <- i
  
  df_sigma <- rbind(df_sigma, df_sigma_temp)
}

df_sigma_plot <- df_sigma %>%
  dplyr::filter(Method == "actual") %>%
  dplyr::select(Pos, value) %>%
  dplyr::rename(actual = value) %>%
  dplyr::right_join(df_sigma, by = "Pos") %>%
  dplyr::mutate(Error = abs(actual - value)) %>%
  dplyr::select(Pos, value, Method, Error)

df_sigma_plot <- df_sigma_plot %>% arrange(Pos)

# plot
ggplot(df_sigma_plot, aes(Pos, round(value, digits = 3), class = Method, color = Method, size = Error^1.5)) +
  geom_point(alpha = 0.80) + 
  labs(title= "Accuracy of Covariance", 
       subtitle="For each Priming Method",
       x= "Player Position",
       y= "Avg. Covariance Estimates") +
  guides(size = guide_legend(title="Mean Absolute Error"))

#iterations plot
v_iterations <- length(mvn_base$all.loglik)

df_iter <- data.frame(Method = "random", iterations = v_iterations)

for (i in names(mvn_clusters_nba)) {
  
  em_temp <- mvn_clusters_nba[[i]]
  
  df_iter <- rbind(df_iter, data.frame(Method = i, iterations = length(em_temp$all.loglik)))
}

ggplot(data = df_iter, aes(x = reorder(Method, -iterations), y = iterations, fill = rainbow(3))) +
  geom_bar(stat="identity") +
  ylab("Number of Iterations") +
  xlab("Method") +
  labs(title = "Bar plot of Iterations", subtitle = "By Method") +
  theme(legend.position = "none")  
