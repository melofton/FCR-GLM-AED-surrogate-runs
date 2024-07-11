# Title: Assemble Surrogate Dataset
# Author: Mary Lofton
# Date: 11JUL24

# Purpose: Create a dataset for Bobby's students to run with for surrogate development

# Notes:

# Data frame format: Parameter 1, parameter 2, temp/weather deviation, datetime, 
# variable, prediction, observation

# load packages ----
source("./Scripts/install.R")

library(tidyverse)
library(lubridate)
library(lhs)
library(glmtools)
library(GLM3r)

# set working directory; you can change this to be any calibration folder ----
setwd("./GLM-AED") 

# create matrix of parameter values using maximin space-filling design ----
phyto_groups <- c("hot","cold","Nfixer")
param_names <- c("pd%R_growth", "pd%w_p", "pd%I_K", "pd%R_resp", "pd%KePHY", "pd%K_N", "pd%K_P")

# Latin hypercube function
# n runs for m factors
# so start with 1000 runs for 7*3 = 21 factors (parameters)

mylhs <- function(n, m)
{
  ## generate the Latin hypercube 
  l <- (-(n - 1)/2):((n - 1)/2)
  L <- matrix(NA, nrow=n, ncol=m)
  for(j in 1:m) L[,j] <- sample(l, n)
  
  ## draw the random uniforms and turn the hypercube into a sample
  U <- matrix(runif(n*m), ncol=m)
  X <- (L + (n - 1)/2 + U)/n
  colnames(X) <- paste0("x", 1:m)
  
  ## return the design and the grid it lives on for visualization
  return(list(X=X, g=c((l + (n - 1)/2)/n,1)))
}

Dlist <- mylhs(n = 1000, m = 21)

# data wrangling to get parameter values in correct range
scale_R_growth <- function(x, na.rm = FALSE) x*3 + 0.5
scale_w_p <- function(x, na.rm = FALSE){
  
  for(i in 1:length(x)){
  if(x[i] == 0.5){x[i] <- 0} else if(x[i] < 0.5){x[i] <- x[i]*-2} else {x[i] <- (x[i]-0.5)*2}
  }
  
  return(x)
  
}
scale_I_K <- function(x, na.rm = FALSE) x*490 + 10
scale_R_resp <- function(x, na.rm = FALSE) x*0.24 + 0.01
scale_KePHY <- function(x, na.rm = FALSE) x*0.01
scale_K_N <- function(x, na.rm = FALSE) x*4.9 + 0.1
scale_K_P <- function(x, na.rm = FALSE) x*0.49 + 0.01



param_values <- tibble(data.frame(Dlist$X)) %>%
  mutate_at(c("x1","x2","x3"), scale_R_growth) %>%
  mutate_at(c("x4","x5","x6"), scale_w_p) %>%
  mutate_at(c("x7","x8","x9"), scale_I_K) %>%
  mutate_at(c("x10","x11","x12"), scale_R_resp) %>%
  mutate_at(c("x13","x14","x15"), scale_KePHY) %>%
  mutate_at(c("x16","x17","x18"), scale_K_N) %>%
  mutate_at(c("x19","x20","x21"), scale_K_P) 
colnames(param_values) <- c("R_growth_hot","R_growth_cold","R_growth_Nfixer",
                            "w_p_hot","w_p_cold","w_p_Nfixer",
                            "I_K_hot","I_K_cold","I_K_Nfixer",
                            "R_resp_hot","R_resp_cold","R_resp_Nfixer",
                            "KePHY_hot","KePHY_cold","KePHY_Nfixer",
                            "K_N_hot","K_N_cold","K_N_Nfixer",
                            "K_P_hot","K_P_cold","K_P_Nfixer")


# set nml filepath
nml_file <- file.path('./aed/aed2_phyto_pars_24MAY24_MEL.nml')

# set file location of output
nc_file <- file.path('./output/output.nc') 

# save starting version of nml in environment so you can reset after
start_nml <- glmtools::read_nml(nml_file = nml_file)

# for-loop to run GLM using different parameter values

  
  for(j in 1:length(unlist(param_values[,1]))){
    
  # read in nml
  nml <- glmtools::read_nml(nml_file = nml_file)
  
  # replace parameter value as desired
  curr_R_growth <- unname(unlist(param_values[j,c(1:3)]))
  curr_w_p <- unname(unlist(param_values[j,c(4:6)]))
  curr_I_K <- unname(unlist(param_values[j,c(7:9)]))
  curr_R_resp <- unname(unlist(param_values[j,c(10:12)]))
  curr_KePHY <- unname(unlist(param_values[j,c(13:15)]))
  curr_K_N <- unname(unlist(param_values[j,c(16:18)]))
  curr_K_P <- unname(unlist(param_values[j,c(19:21)]))
  
  # set nml parameter values
  new_nml <- glmtools::set_nml(nml, arg_name = param_names[1], arg_val = curr_R_growth)
  new_nml1 <- glmtools::set_nml(new_nml, arg_name = param_names[2], arg_val = curr_w_p)
  new_nml2 <- glmtools::set_nml(new_nml1, arg_name = param_names[3], arg_val = curr_I_K)
  new_nml3 <- glmtools::set_nml(new_nml2, arg_name = param_names[4], arg_val = curr_R_resp)
  new_nml4 <- glmtools::set_nml(new_nml3, arg_name = param_names[5], arg_val = curr_KePHY)
  new_nml5 <- glmtools::set_nml(new_nml4, arg_name = param_names[6], arg_val = curr_K_N)
  new_nml6 <- glmtools::set_nml(new_nml5, arg_name = param_names[7], arg_val = curr_K_P)
  
  # create path to write permuted nml to file
  write_path <- nml_file
  
  # write permuted nml to file
  glmtools::write_nml(new_nml6, file = write_path)
  
  # run GLM-AED using GLM3r
  GLM3r::run_glm()

  # pull variable of interest from model output
  var <- glmtools::get_var(nc_file, var_name = "PHY_tchla", reference="surface", z_out=1.6) %>%
    filter(hour(DateTime) == 12)
  
  # pull parameters from model output
  R_growth <- new_nml6$phyto_data$`pd%R_growth`
  w_p <- new_nml6$phyto_data$`pd%w_p`
  I_K <- new_nml6$phyto_data$`pd%I_K`
  R_resp <- new_nml6$phyto_data$`pd%R_resp`
  KePHY <- new_nml6$phyto_data$`pd%KePHY`
  K_N <- new_nml6$phyto_data$`pd%K_N`
  K_P <- new_nml6$phyto_data$`pd%K_P`
  
  # assemble dataframe for that model run
  temp <- data.frame(R_growth_hot = R_growth[1],
                     R_growth_cold = R_growth[2],
                     R_growth_Nfixer = R_growth[3],
                     w_p_hot = w_p[1],
                     w_p_cold = w_p[2],
                     w_p_Nfixer = w_p[3],
                     I_K_hot = I_K[1],
                     I_K_cold = I_K[2],
                     I_K_Nfixer = I_K[3],
                     R_resp_hot = R_resp[1],
                     R_resp_cold = R_resp[2],
                     R_resp_Nfixer = R_resp[3],
                     KePHY_hot = KePHY[1],
                     KePHY_cold = KePHY[2],
                     KePHY_Nfixer = KePHY[3],
                     K_N_hot = K_N[1],
                     K_N_cold = K_N[2],
                     K_N_Nfixer = K_N[3],
                     K_P_hot = K_P[1],
                     K_P_cold = K_P[2],
                     K_P_Nfixer = K_P[3],
                     deviation = 0,
                     datetime = var$DateTime,
                     variable = "PHY_tchla_1.6",
                     prediction = var$PHY_tchla_1.6)

  # make sure you reset nml
  glmtools::write_nml(start_nml, file = nml_file)
  
  # bind to other model runs
  if(j == 1){
    final <- temp
  } else {
    final <- bind_rows(final, temp)
  }

    }

write.csv(final, file = "./model_runs_surrogate.csv",row.names = FALSE)



# # Code needed for Parul runs eventually
# # read in observation dataset
# obs <- read_csv("/home/rstudio/RProjects/FCR-GLM-metrics/observations/EXOChla.csv") %>%
#   filter(Depth == 1.6) %>%
#   rename(observed = PHY_TCHLA) %>%
#   select(-Depth)
# 
# # write final dataset to file
# final1 <- final %>%
#   filter(datetime %in% obs$DateTime)
# write.csv(final1, file = "./collated_model_scenarios_EXOdates.csv",row.names = FALSE)
# 
# # plot parameter space
# parms <- final1 %>%
#   select(R_growth_cyano, w_p_cyano) %>%
#   distinct()
# 
# ggplot(data = parms, aes(x = R_growth_cyano, y = w_p_cyano))+
#   geom_point(size = 3)+
#   theme_bw()
# 
# # plot output
# plot_data <- final1 %>%
#   group_by(R_growth_cyano, w_p_cyano) %>%
#   mutate(model_run = cur_group_id()) %>%
#   rename(DateTime = datetime) %>%
#   left_join(., obs, by = "DateTime")
# 
# p <- ggplot(data = plot_data)+
#   geom_line(aes(x = DateTime, y = prediction, group = as.factor(model_run), color = as.factor(R_growth_cyano), linetype = as.factor(w_p_cyano)))+
#   geom_point(aes(x = DateTime, y = observed), size = 0.5)+
#   xlab("")+
#   ylab("chlorophyll-a (ug/L) at 1.6 m")+
#   labs(color = "Cyano growth rate", linetype = "Cyano sinking rate")+
#   theme_classic()+
#   theme(legend.position = "bottom")
# p
# ggsave(p, filename = "./figures/example_param_scenarios.png", units = "in",
#        dev = "png", height = 3, width = 7)
# 
# # look at f factors
# glmtools::sim_vars(file = nc_file)
# 
# f_factor_names <- c("PHY_green_fI","PHY_green_fNit","PHY_green_fPho","PHY_green_fSal",
#                "PHY_green_fSil","PHY_green_fT","PHY_cyano_fI","PHY_cyano_fNit","PHY_cyano_fPho","PHY_cyano_fSal",
#                "PHY_cyano_fSil","PHY_cyano_fT","PHY_diatom_fI","PHY_diatom_fNit","PHY_diatom_fPho","PHY_diatom_fSal",
#                "PHY_diatom_fSil","PHY_diatom_fT")
# 
# 
# 
# plot_factors <- f_factors %>%
#   pivot_longer(PHY_green_fI_1.6:PHY_diatom_fT_1.6, names_to = "var_name", values_to = "f_factor") %>%
#   separate(var_name, c("PHY","group","factor_name","depth1","depth2")) %>%
#   select(-c("PHY","depth1","depth2")) %>%
#   filter(DateTime >= "2018-08-06") 
# 
# factor_plot <- ggplot(data = plot_factors, aes(x = DateTime, y = f_factor, group = factor_name, color = factor_name))+
#   geom_line()+
#   facet_grid(rows = vars(group))+
#   theme_bw()
# factor_plot