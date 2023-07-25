

#library("odin.dust")

#library("patchwork")
#library('mcstate')
#library(didehpc)
#library(pkgdepends)

#library("coda")
#library(binom)

#library(bayesplot)
#library(reshape2)
#library(ggpubr)
#library(patchwork)
#library(RColorBrewer)
install.packages("luz")
## actually used libraries..
library(odin)
library(ggplot2)
library(sifter)
library(zoo)
library(dplyr)
library(torch)
library(luz)
library(tidyverse)

??load_file
##Set default theme
theme_set(theme_minimal()+
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  strip.background = element_blank(),
                  panel.border = element_rect(colour = "black",fill=NA)))

#Required functions
source('sim/data_gen.R')
source('shared/run_pmcmc.R')
source('shared/plot_particle_filter.R')
source('shared/addCIs.R')
source('shared/model_parameters.R')
source('shared/equilibrium-init-create-stripped.R')

##Generate simulated data##

n_sims<-20
volatility<-0.8
init_EIR<-20
duration=20*365
out_step=30

generate_sim_compendium<-function(n_sims,volatility,init_EIR,duration,out_step){
  sims_compendium<-data.frame(
    run=numeric(),t=numeric(),prev_true=numeric(),EIR_true=numeric(),vol_true=numeric(),inc_true=numeric()
  )
  for(i in 1:n_sims){
    sims_compendium<-sims_compendium%>%add_row(
      cbind(
        data.frame(run=i),
          data_gen(EIR_volatility = volatility, init_EIR = init_EIR,time = duration,out_step = out_step)
      )
    )
  }
  return(sims_compendium)
}

#generate_torch_ds_lags<-function(sim_compendium,fit_var,pred_var,lag_no){
  lags=1:lag_no
  map_lag <- lags %>% map(~partial(lag, n = .x))
  fit_tensor<-torch_tensor(
    as.matrix(
      sim_compendium%>%group_by(run)%>%
        mutate(across(.cols = {{fit_var}}, .fns = map_lag, .names = "{.col}_lag{lags}"))%>%
        filter(step>lag_no)%>%
        ungroup()%>%
        select(contains("lag"))
    )
  )
  pred_tensor<-torch_tensor(
    as.matrix(sim_compendium%>%
                group_by(run)%>%
                filter(step>lag_no)%>%
                ungroup()%>%
                select(pred_var)))
  return(tensor_dataset(fit_tensor, pred_tensor))
}


generate_torch_ds_lags<-function(sim_compendium,fit_var,pred_var,t_var,lag_no){
  lags=1:lag_no
  map_lag <- lags %>% map(~partial(lag, n = .x))
  fit_tensor<-torch_tensor(
    as.matrix(
      sim_compendium%>%group_by(run)%>%
        arrange(run,!!sym(paste0(t_var)))%>%
        mutate(across(.cols = {{fit_var}}, .fns = map_lag, .names = "{.col}_forward{lags}"))%>%
        filter(!!sym(t_var)>lag_no)%>%
        ungroup()%>%
        select(contains("lag"))
    )
  )
  pred_tensor<-torch_tensor(
    as.matrix(sim_compendium%>%
                group_by(run)%>%
                arrange(run,!!sym(t_var))%>%
                filter(!!sym(t_var)>lag_no)%>%
                ungroup()%>%
                select(pred_var)))
  return(tensor_dataset(fit_tensor, pred_tensor))
}

predict_prev_inf_lags<-generate_preds_valid_lag(model=net,fit_var="prev_true",pred_var = "inc_true_1000",t_var="step_back",
                                                epochs=50, sims_compendium_train=sims_compendium_train,sims_compendium_test = sims_compendium_test,sims_compendium_valid = sims_compendium_valid
)
compare<-sims_compendium_test%>%
  group_by(run)%>%
  filter(step_back>24)%>%
  select(t,inc_true_1000)
compare<-sims_compendium_test%>%
  group_by(run)%>%
  arrange(run,step_back)%>%
  filter(step_back>24)%>%
  select(t,inc_true_1000)
compare$pred<-predict_prev_inf_lags$predictions


ggplot(compare,aes(x=t,y=inc_true_1000))+
  geom_line()+geom_line(aes(y=pred),col="red")+
  facet_wrap(~run)+theme(strip.text.x = element_text(size=0))

length(compare$t)

fit_var = "prev_true"
pred_var = "inc_true"
t_var="step_back"
sims_compendium_train%>%group_by(run)%>%
  arrange(run,!!sym(t_var))%>%
  mutate(across(.cols = {{fit_var}}, .fns = map_lag, .names = "{.col}_lag{lags}"))%>%
  filter(!!sym(t_var)>lag_no)%>%
  ungroup()

sims_compendium_train%>%
  group_by(run)%>%
  arrange(run,!!sym(t_var))%>%
  filter(!!sym(t_var)>lag_no)%>%
  ungroup()%>%
  select(pred_var)

sims_compendium_train

sims_compendium_train%>%group_by(run)%>%
  arrange(run,!!sym(t_var))%>%
  mutate(across(.cols = {{fit_var}}, .fns = map_lag, .names = "{.col}_lag{lags}"))%>%
  filter(!!sym(t_var)>lag_no)%>%
  ungroup()

sims_compendium_train
check_torch<-generate_torch_ds_lags(sim_compendium = sims_compendium_train,fit_var = "prev_true",pred_var = "inc_true_1000",t_var="step_back",lag_no=24)
check_torch$tensors[[1]]

sims_compendium_train$inc_true_1000
generate_preds_valid_lag<-function(model,fit_var,pred_var,t_var,sims_compendium_train,sims_compendium_test,sims_compendium_valid,
                                   d_hidden=100,epochs=50,loss=nn_mse_loss(),optimizer=optim_adam,lag_no=24){
  
  train_ds<-generate_torch_ds_lags(sim_compendium = sims_compendium_train,fit_var = fit_var,pred_var = pred_var,t_var=t_var,lag_no=lag_no)
  test_ds<-generate_torch_ds_lags(sim_compendium = sims_compendium_test,fit_var = fit_var,pred_var = pred_var,t_var=t_var,lag_no=lag_no)
  valid_ds<-generate_torch_ds_lags(sim_compendium = sims_compendium_valid,fit_var = fit_var,pred_var = pred_var,t_var=t_var,lag_no=lag_no)
  train_dl<-dataloader(train_ds, batch_size = 100, shuffle = TRUE)
  test_dl<-dataloader(test_ds, batch_size = 100, shuffle = TRUE)
  valid_dl<-dataloader(valid_ds, batch_size = 100, shuffle = TRUE)
  
  # output dimensionality (number of predicted features)
  d_out <- ncol(train_ds$tensors[[2]])
  d_in <-  ncol(train_ds$tensors[[1]])
  n<-nrow(train_ds$tensors[[1]])
  
  fitted <- model %>%
    setup(loss = loss, 
          optimizer = optimizer,
          metrics = list(luz_metric_mae())) %>%
    set_hparams(
      d_in = d_in,
      d_hidden = d_hidden, d_out = d_out
    ) %>%
    fit(train_ds, epochs = epochs, valid_data = valid_dl)
  
  pred<-fitted %>% predict(test_ds)
  # sims_compendium_test$pred<-as.vector(t(as.matrix(pred)))
  return(list(fit_model=fitted,
              predictions=as.vector(t(as.matrix(pred))))
  )
}



net <- nn_module(
  initialize = function(d_in, d_hidden, d_out) {
    self$net <- nn_sequential(
      nn_linear(d_in, d_hidden),
      nn_relu(),
      nn_linear(d_hidden, d_hidden),
      nn_relu(),
      nn_linear(d_hidden, d_out),
        )
  },
  forward = function(x) {
    self$net(x)
  }
)

sims_compendium_train<-generate_sim_compendium(n_sims,volatility,init_EIR,duration,out_step)
sims_compendium_train<-sims_compendium_train%>%
  mutate(inc_true_1000=inc_true*1000,
         prev_true_zeroed=replace(prev_true,prev_true<0,1e-12),
         log_odds_prev=log(prev_true_zeroed/(1-prev_true_zeroed)))

sims_compendium_test<-generate_sim_compendium(n_sims,volatility,init_EIR,duration,out_step)
sims_compendium_test<-sims_compendium_test%>%
  mutate(inc_true_1000=inc_true*1000,
         prev_true_zeroed=replace(prev_true,prev_true<0,1e-12),
         log_odds_prev=log(prev_true_zeroed/(1-prev_true_zeroed)))

sims_compendium_valid<-generate_sim_compendium(n_sims,volatility,init_EIR,duration,out_step)
sims_compendium_valid<-sims_compendium_valid%>%
  mutate(inc_true_1000=inc_true*1000,
         prev_true_zeroed=replace(prev_true,prev_true<0,1e-12),
         log_odds_prev=log(prev_true_zeroed/(1-prev_true_zeroed)))

sims_compendium_test$step<-sims_compendium_test$t/30
sims_compendium_train$step<-sims_compendium_train$t/30
sims_compendium_valid$step<-sims_compendium_valid$t/30
sims_compendium_train$step_back<-(max(sims_compendium_train$t)-sims_compendium_train$t)/30+1
sims_compendium_test$step_back<-(max(sims_compendium_test$t)-sims_compendium_test$t)/30+1
sims_compendium_valid$step_back<-(max(sims_compendium_valid$t)-sims_compendium_valid$t)/30+1
sims_compendium_train$step_back<-(max(sims_compendium_train$t)-sims_compendium_train$t)/30+1


max(predict_prev_inf_lags$predictions)



compare


########## after here isn't relevant#####

check<-data.frame(t=1:10,y=rep(1:2,5),x=rnorm(10))
order_data<-function(data,arg_var){
  return(data%>%arrange(y,!!sym(arg_var))%>%group_by(y))
}
order_data(check,"x")



generate_torch_ds<-function(sim_compendium,fit_var,pred_var){
  fit_tensor<-torch_tensor(
    as.matrix(sim_compendium%>%
                select(t,run,fit_var)%>%
                pivot_wider(names_from=t,values_from = fit_var)%>%
                select(-run)
              )
    )
  pred_tensor<-torch_tensor(
    as.matrix(sim_compendium%>%
                select(t,run,pred_var)%>%
                pivot_wider(names_from=t,values_from = pred_var)%>%
                select(-run)
    )
  )
  return(
    tensor_dataset(fit_tensor, pred_tensor)
  )
}

generate_lagged_matrix <- function(y, lag_range) {
  n <- length(y)
  lagged_matrix <- matrix(NA, nrow = n - max(lag_range) + 1, ncol = length(lag_range))
  
  for (i in 1:length(lag_range)) {
    lag <- lag_range[i]
    lagged_matrix[, i] <- y[(lag + 1):(n - max(lag_range) + lag)]
  }
  
  return(lagged_matrix)
}

set.seed(42)  # For reproducibility
t <- 1:100
y <- cumsum(runif(length(t)))

# Generate the lagged matrix for t = 10:100 with lags 1 to 10
lag_range <- 10
lagged_matrix <- generate_lagged_matrix(y[10:100], lag_range)

# Print the lagged matrix
print(lagged_matrix)



predict_inc_prev_lags$predictions


show_train<-sims_compendium_train%>%
  group_by(run)%>%
  filter(step>24)%>%
  select(t,prev_true)

ggplot(show_train,aes(x=t,y=prev_true))+
  geom_line()+
  facet_wrap(~run)+theme(strip.text.x = element_text(size=0))


compare$pred<-predict_inc_prev_lags$predictions

predict_inc_prev_lag


sims_compendium_test$step=sims_compendium_test$t/30

  
  sims_compendium_test
  fit_tensor<-torch_tensor(
    as.matrix(sim_compendium%>%
                filter(t>lag_no)%>%
                select(t,run,fit_var)%>%
                pivot_wider(names_from=t,values_from = fit_var)%>%
                select(contains(la))
    )
  )
  pred_tensor<-torch_tensor(
    as.matrix(sim_compendium%>%
                select(t,run,pred_var)%>%
                pivot_wider(names_from=t,values_from = pred_var)%>%
                select(-run)
    )
  )
  return(
    tensor_dataset(fit_tensor, pred_tensor)
  )
}


generate_torch_ds<-function(sim_compendium,fit_var,pred_var,lag_no){
  fit_tensor<-torch_tensor(
    as.matrix(sim_compendium%>%
                select(t,run,fit_var)%>%
                pivot_wider(names_from=t,values_from = fit_var)%>%
                select(-run)
    )
  )
  pred_tensor<-torch_tensor(
    as.matrix(sim_compendium%>%
                select(t,run,pred_var)%>%
                pivot_wider(names_from=t,values_from = pred_var)%>%
                select(-run)
    )
  )
  return(
    tensor_dataset(fit_tensor, pred_tensor)
  )
}

generate_preds<-function(model,fit_var,pred_var,sims_compendium_train,sims_compendium_test,
                         d_hidden=10000,epochs=400,loss=nn_mse_loss(),optimizer=optim_adam){
  
  train_ds<-generate_torch_ds(sim_compendium = sims_compendium_train,fit_var = fit_var,pred_var = pred_var)
  test_ds<-generate_torch_ds(sim_compendium = sims_compendium_test,fit_var = fit_var,pred_var = pred_var)
  
  # output dimensionality (number of predicted features)
  d_out <- ncol(train_ds$tensors[[2]])
  d_in <-  ncol(train_ds$tensors[[1]])
  n<-nrow(train_ds$tensors[[1]])
  
  fitted <- model %>%
  setup(loss = loss, optimizer = optimizer) %>%
  set_hparams(
    d_in = d_in,
    d_hidden = d_hidden, d_out = d_out
  ) %>%
  fit(train_ds, epochs = epochs)

  pred<-fitted %>% predict(test_ds)
  sims_compendium_test$pred<-as.vector(t(as.matrix(pred)))
  return(list(fit_model=fitted,
              predictions_df=sims_compendium_test))
  }
  
net <- nn_module(
  initialize = function(d_in, d_hidden, d_out) {
    self$net <- nn_sequential(
      nn_linear(d_in, d_hidden),
      nn_relu(),
      nn_linear(d_hidden, d_out)
    )
  },
  forward = function(x) {
    self$net(x)
  }
)

check_ds<-generate_torch_ds(sim_compendium = sims_compendium_train,fit_var="inc_true",pred_var = "prev_true")
check_ds$tensors[1]
generate_preds_valid<-function(model,fit_var,pred_var,sims_compendium_train,sims_compendium_test,sims_compendium_valid,
                         d_hidden=10000,epochs=400,loss=nn_mse_loss(),optimizer=optim_adam){
  
  train_ds<-generate_torch_ds(sim_compendium = sims_compendium_train,fit_var = fit_var,pred_var = pred_var)
  test_ds<-generate_torch_ds(sim_compendium = sims_compendium_test,fit_var = fit_var,pred_var = pred_var)
  valid_ds<-generate_torch_ds(sim_compendium = sims_compendium_valid,fit_var = fit_var,pred_var = pred_var)
  train_dl<-dataloader(train_ds, batch_size = 100, shuffle = TRUE)
  test_dl<-dataloader(test_ds, batch_size = 100, shuffle = TRUE)
  valid_dl<-dataloader(valid_ds, batch_size = 100, shuffle = TRUE)
  
  # output dimensionality (number of predicted features)
  d_out <- ncol(train_ds$tensors[[2]])
  d_in <-  ncol(train_ds$tensors[[1]])
  n<-nrow(train_ds$tensors[[1]])
  
  fitted <- model %>%
    setup(loss = loss, 
          optimizer = optimizer,
          metrics = list(luz_metric_mae())) %>%
    set_hparams(
      d_in = d_in,
      d_hidden = d_hidden, d_out = d_out
    ) %>%
    fit(train_ds, epochs = epochs, valid_data = valid_dl)
  
  pred<-fitted %>% predict(test_ds)
  sims_compendium_test$pred<-as.vector(t(as.matrix(pred)))
  return(list(fit_model=fitted,
              predictions_df=sims_compendium_test))
}


predict_inc_prev$fit_model$
sims_compendium_train$EIR_true
predict_inc_prev<-generate_preds(model=net,fit_var="inc_true",pred_var = "prev_true",
                                 sims_compendium_train=sims_compendium_train,sims_compendium_test = sims_compendium_test
                                 )

predict_inc_prev<-generate_preds_valid(model=net,fit_var="inc_true",pred_var = "prev_true",
                                 sims_compendium_train=sims_compendium_train,sims_compendium_test = sims_compendium_test,sims_compendium_valid = sims_compendium_valid
)
log_odds_prev
predict_inc_prev<-generate_preds_valid(model=net,fit_var="inc_true",pred_var = "prev_true",
                                       sims_compendium_train=sims_compendium_train,sims_compendium_test = sims_compendium_test,sims_compendium_valid = sims_compendium_valid
)

sims_compendium_train$log_odds_prev
ggplot(predict_inc_prev$predictions_df,aes(x=t,y=prev_true))+
  geom_line()+geom_line(aes(y=pred),col="red")+
  facet_wrap(~run)+theme(strip.text.x = element_text(size=0))


check$tensors[1]
ds$tensors[1]
pivot_wider(sim_compendium,cols=)

x_mat_test<-sapply(1:max(sims_compendium_test$run),function(i){
  sims_compendium_test$inc_true[which(sims_compendium_test$run==i)]
}
)

y_mat_test<-sapply(1:max(sims_compendium_test$run),function(i){
  sims_compendium_test$prev_true[which(sims_compendium_test$run==i)]
}
)
x_test<-torch_tensor(t(x_mat_test))
y_test<-torch_tensor(t(y_mat_test))


dl_test <- dataloader(ds_test, batch_size = 100, shuffle = TRUE)

n=100
sort(sample.int(
  n = n,
  size = n
))
n_timesteps=12
rep(1:(max(sims_compendium$step)-n_timesteps),length(unique(sims_compendium$run)))
sims_compendium$step<-sims_compendium$t/30
names(sims_compendium$inc_true)
mean_inc<-mean(sims_compendium$inc_true)
sd_inc<-sd(sims_compendium$inc_true)

train_data<-demand_dataset(sims_compendium,"EIR_true","inc_true",mean_inc,sd_inc,24)
train_data$x

y_var<-"inc_true"
vect<-as.double()
check_data$starts
class(vect)

data<-rep(1:10,5)
n=5

# input dimensionality (number of input features)
d_in <- 3
# number of observations in training set
n <- 1000
x
x <- torch_randn(n, d_in)
coefs <- c(0.2, -1.3, -0.5)
y <- x$matmul(coefs)$unsqueeze(2) + torch_randn(n, 1)



x_mat<-sapply(1:max(sims_compendium$run),function(i){
  sims_compendium$inc_true[which(sims_compendium$run==i)]
}
)

y_mat<-sapply(1:max(sims_compendium$run),function(i){
  sims_compendium$prev_true[which(sims_compendium$run==i)]
}
)
x<-torch_tensor(t(x_mat))
y<-torch_tensor(t(y_mat))

ds <- tensor_dataset(x, y)

dl <- dataloader(ds, batch_size = 100, shuffle = TRUE)





ds$tensors
pred<-fitted %>% predict(ds_test)

plot(ds_test$tensors[[2]][2,])
lines(pred[2,])
lines(pred[1,])
(pred_vect)
sims_compendium_test$inc_true
as.numeric(torch_flatten(t(pred)))



windows(height=10,width=15)
ggplot(sims_compendium_test,aes(x=t,y=prev_true))+
  geom_line()+geom_line(aes(y=pred),col="red")+
  facet_wrap(~run)+theme(strip.text.x = element_text(size=0))


sims_compendium_test
pred_vect<-pred$view()
t$t()$reshape(9)

predict()
fitted <- net %>%
  setup(loss = nn_mse_loss(), optimizer = optim_adam) %>%
  set_hparams(
    d_in = d_in,
    d_hidden = d_hidden, d_out = d_out
  ) %>%
  fit(dl, epochs = 200)

plot(fitted)
evaluate(fitted, dl_test)
preds <- predict(fitted, dl_test)
plot(preds[1,])
lines(dl_test$dataset[1][[2]])

fitted %>% predict(dl_test)

y_mat_test[,1]
x <- torch_randn(n, d_in)
coefs <- c(0.2, -1.3, -0.5)
y <- x$matmul(coefs)$unsqueeze(2) + torch_randn(n, 1)



ds <- tensor_dataset(x, y)
ds[1]

unlist(lapply(1:(max(sims_compendium$step)-n_timesteps),function(i){
  which(sims_compendium$step==i)
}
))
## set up the data ######################
demand_dataset <- dataset(
  name = "demand_dataset",
  initialize = function(sims_compendium,x_var,y_var,train_mean,train_sd,n_timesteps) {
    self$n_timesteps <- n_timesteps
    self$x <- torch_tensor(as.matrix(sims_compendium%>%select(x_var)))
    self$y<-torch_tensor((as.matrix(sims_compendium%>%select(y_var)) - train_mean) / train_sd)
    self$starts <- unlist(lapply(1:(max(sims_compendium$step)-n_timesteps),function(i){
      which(sims_compendium$step==i)
    }
    ))
  },
  .getitem = function(i) {
    start <- self$starts[i]
    end <- start + self$n_timesteps - 1
    
    list(
      x = self$x[start:end],
      y = self$y[end + 1]
    )
  },
  .length = function() {
    length(self$starts)
  }
)

model <- nn_module(
  initialize = function(input_size,
                        hidden_size,
                        dropout = 0.2,
                        num_layers = 1,
                        rec_dropout = 0) {
    self$num_layers <- num_layers
    
    self$rnn <- nn_lstm(
      input_size = input_size,
      hidden_size = hidden_size,
      num_layers = num_layers,
      dropout = rec_dropout,
      batch_first = TRUE
    )
    
    self$dropout <- nn_dropout(dropout)
    self$output <- nn_linear(hidden_size, 1)
  },
  forward = function(x) {
    (x %>%
       # these two are equivalent
       # (1)
       # take output tensor,restrict to last time step
       self$rnn())[[1]][, dim(x)[2], ] %>%
      # (2)
      # from list of state tensors,take the first,
      # and pick the final layer
      # self$rnn())[[2]][[1]][self$num_layers, , ] %>%
      self$dropout() %>%
      self$output()
  }
)

input_size <- 1
hidden_size <- 32
num_layers <- 2
rec_dropout <- 0.2

model <- model %>%
  setup(optimizer = optim_adam, loss = nn_mse_loss()) %>%
  set_hparams(
    input_size = input_size,
    hidden_size = hidden_size,
    num_layers = num_layers,
    rec_dropout = rec_dropout
  )

rates_and_losses <- model %>% 
  lr_finder(train_data, start_lr = 1e-3, end_lr = 1)
rates_and_losses %>% plot()


ds[1]

ggplot(sims_compendium,aes(x=t,y=log(EIR_true),by=as.factor(run)))+
  geom_line()

ggplot(sims_compendium,aes(x=t,y=prev_true,by=as.factor(run)))+
  geom_line()


plot(data_sim_comptest$EIR_true,type = 'l')

##Three previously run simulated data sets are saved in the folder
## 'anc_pmcmc/sim/sim_datasets
data_sim_comptest3 <- readRDS('sim/sim_datasets/data_sim3.RDS')

##Test run_pmcmc function##
test_run <- run_pmcmc(data = data_sim_comptest3,
                      n_particles = 10,
                      proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
                      max_EIR=1000,
                      max_steps = 1e7,
                      atol = 1e-5,
                      rtol = 1e-6,
                      n_steps = 10,
                      n_threads = 2,
                      lag_rates = 10)
plot_particle_filter(test_run$history,true_history=data_sim_comptest3,times=data_sim_comptest3$t)

##Set up cluster##
root <- "T:/jth/contexts" ##Edit to your contexts path
sources <- c("sim/run_pmcmc.R",
             "shared/model_parameters.R","shared/equilibrium-init-create-stripped.R")


ctx <- context::context_save("T:/jth/contexts", sources = sources,
                             packages = c('statmod','coda'),
                             package_sources = conan::conan_sources(c("mrc-ide/odin.dust",'mrc-ide/mcstate')))

##Set up and run on 32 node core##
config_32 <- didehpc::didehpc_config(cores = 32, parallel = TRUE)
obj_32 <- didehpc::queue_didehpc(ctx,config = config_32)

obj_32$cluster_load(TRUE)
obj_32$config
obj_32$login()

run_32_200 <- obj_32$enqueue(run_pmcmc(data = data_sim_comptest3,
                                      n_particles = 200,
                                      proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
                                      max_EIR=1000,
                                      max_steps = 1e7,
                                      atol = 1e-5,
                                      rtol = 1e-6,
                                      n_steps = 1000,
                                      n_threads = 32))
run_32_200$status()
run_32_200$log()

##Produce some diagnostics
1 - coda::rejectionRate(as.mcmc(result_32_200$mcmc)) ##Acceptance rate
coda::effectiveSize(as.mcmc(result_32_200$mcmc)) ##ESS
cov(result_32_200$pars) ##Covariance
summary(as.mcmc(result_32_200$mcmc)) ##Summarize mcmc run
plot(as.mcmc(result_32_200$mcmc)) ##Plot traces and distributions

##Plot particle filter results vs simulated data
plot_particle_filter(result_32_200$history,true_history=data_sim_comptest3,times=data_sim_comptest3$t)

##Below creates various figures to show prevalence, EIR, and incidence trajectories
##Used for ASTMH presentation
prev_cis <- addCIs(data_sim_comptest3,data_sim_comptest3$positive,data_sim_comptest3$tested)
##Remove 5% burn-in and reformat prevalence output
prev_history <- data.frame(t(result_32_200$history['prev', 51:1000, -1]))
long_prev_history <- prev_history%>%
  mutate(time=c(1:nrow(prev_history))*30)%>%
  melt(id='time')

windows(7,5)
prev_plot <- ggplot(long_prev_history)+
  geom_line(aes(x=time,y=value,group=variable),col = "#A6CEE3",alpha=0.2)+
  geom_point(data=prev_cis,aes(x=t,y=mean),pch = 19,
             col = "#666666")+
  geom_errorbar(data=prev_cis,aes(x=t,ymin=lower,ymax=upper),width = 0,
                col = "#666666")+
  scale_y_continuous(limits=c(0,1))+
  labs(x='Time (days)',y='Prevalence')
prev_plot_true <- ggplot(long_prev_history)+
  geom_line(aes(x=time,y=value,group=variable),col = "#A6CEE3",alpha=0)+
  geom_line(data=data_sim_comptest3,aes(x=t,y=prev_true),size=1,col = "#1F78B4",alpha=1)+
  scale_y_continuous(limits=c(0,1))+
  labs(x='Time (days)',y='Prevalence')
prev_plot_true
prev_data_only <- ggplot(prev_cis)+
  geom_point(aes(x=t,y=mean),pch = 19,
             col = "#666666")+
  geom_errorbar(aes(x=t,ymin=lower,ymax=upper),width = 0,
                col = "#666666")+
  scale_y_continuous(limits=c(0,1))+
  labs(x='Time (days)',y='Prevalence')
prev_data_andtrue <- ggplot(prev_cis)+
  geom_point(aes(x=t,y=mean),pch = 19,
             col = "#666666")+
  geom_errorbar(aes(x=t,ymin=lower,ymax=upper),width = 0,
                col = "#666666")+
  geom_line(data=data_sim_comptest3,aes(x=t,y=prev_true),size=1,col = "#1F78B4",alpha=1)+
  scale_y_continuous(limits=c(0,1))+
  labs(x='Time (days)',y='Prevalence')
prev_plot_blank <- ggplot(long_prev_history)+
  geom_line(aes(x=time,y=value,group=variable),col = "#A6CEE3",alpha=0)+
  scale_y_continuous(limits=c(0,1))+
  labs(x='Time (days)',y='Prevalence')
prev_plot_blank
eir_history <- data.frame(t(result_32_200$history['EIR', 51:1000, -1]))
long_eir_history <- eir_history%>%
  mutate(time=c(1:nrow(eir_history))*30)%>%
  melt(id='time')
windows(7,5)
true_eir <- data.frame(EIR_true=rep(data_sim_comptest3$EIR_true,rep(30,61)))
true_eir$time <- 31:1860
eir_plot <- ggplot(long_eir_history)+
  geom_line(aes(x=time,y=value,group=variable),col = "#A6CEE3",alpha=0.2)+
  geom_line(data=true_eir[1:1800,],aes(x=time,y=EIR_true),size=1,
            col = "#1F78B4")+
  scale_y_continuous(limits=c(0,max(c(long_eir_history$value))))+
  scale_y_log10(breaks=c(.001,0.01,.1,1,10,100,1000),labels=c(.001,0.01,.1,1,10,100,1000))+
  labs(x='Time (days)',y='EIR')
eir_plot
windows(7,5)
true_eir_plot <- ggplot(true_eir[1:1800,])+
  geom_line(aes(x=time,y=EIR_true), size=1, col = "#1F78B4")+
  scale_y_log10(limits=c(min(long_eir_history$value),max(long_eir_history$value)),breaks=c(.001,0.01,.1,1,10,100,1000),labels=c(.001,0.01,.1,1,10,100,1000))+
  labs(x='Time (days)',y='EIR')
true_eir_plot
eir_plot_blank <- ggplot(long_eir_history)+
  geom_line(aes(x=time,y=value,group=variable),col = "#A6CEE3",alpha=0)+
  scale_y_continuous(limits=c(0,max(c(long_eir_history$value))))+
  scale_y_log10(breaks=c(.001,0.01,.1,1,10,100,1000),labels=c(.001,0.01,.1,1,10,100,1000))+
  labs(x='Time (days)',y='EIR')
eir_plot_blank

true_eir_plot <- ggplot(true_eir[1:1800,])+
  geom_line(aes(x=time,y=EIR_true),size=1,
            col = "#666666")+
  scale_y_log10(limits=c(min(long_eir_history$value),max(long_eir_history$value)),breaks=c(.001,0.01,.1,1,10,100,1000),labels=c(.001,0.01,.1,1,10,100,1000))+
  labs(x='Time (days)',y='EIR')
true_eir_plot

inc_history <- data.frame(t(result_32_200$history['inc', 51:1000, -1]))
long_inc_history <- inc_history%>%
  mutate(time=c(1:nrow(eir_history))*30)%>%
  melt(id='time')
inv_med <- long_inc_history%>%
  group_by(time)%>%
  summarise(med_inc=median(value))
windows(7,5)
inc_plot <- ggplot(long_inc_history)+
  geom_line(aes(x=time,y=value*10000,group=variable),col = "#A6CEE3",alpha=0.2)+
  geom_line(data=data_sim_comptest3,aes(x=t+30,y=inc_true*10000),size=1,col = "#1F78B4",alpha=1)+
  scale_y_continuous(limits=c(0,max(long_inc_history$value)*10000))+
  labs(x='Time (days)',y='Incidence (<5 year old)\nper 10,000 persons')
inc_plot
inc_plot_blank <- ggplot(long_inc_history)+
  geom_line(aes(x=time,y=value*10000,group=variable),col = "#A6CEE3",alpha=0)+
  geom_line(data=inv_med,aes(x=time,y=med_inc*10000),col = "#1F78B4",alpha=0)+
  scale_y_continuous(limits=c(0,max(long_inc_history$value)*10000))+
  labs(x='Time (days)',y='Incidence (<5 year old)\nper 10,000 persons')
inc_plot_blank

inc_plot_true <- ggplot(long_inc_history)+
  geom_line(aes(x=time,y=value*10000,group=variable),col = "#A6CEE3",alpha=0)+
  geom_line(data=data_sim_comptest3,aes(x=t+30,y=inc_true*10000),size=1,col = "#1F78B4",alpha=1)+
  scale_y_continuous(limits=c(0,max(long_inc_history$value)*10000))+
  labs(x='Time (days)',y='Incidence (<5 year old)\nper 10,000 persons')
inc_plot_true

windows(15,5)
eir_plot_blank + inc_plot_blank + prev_plot_blank
windows(15,5)
true_eir_plot + inc_plot_blank + prev_plot_blank
true_eir_plot + inc_plot_true + prev_plot_blank
true_eir_plot + inc_plot_true + prev_plot_true
true_eir_plot + inc_plot_true + prev_data_andtrue
true_eir_plot + inc_plot_true + prev_data_only
true_eir_plot + inc_plot_true + prev_plot
true_eir_plot + inc_plot + prev_plot
eir_plot + inc_plot_true+ prev_plot
eir_plot + inc_plot + prev_plot
