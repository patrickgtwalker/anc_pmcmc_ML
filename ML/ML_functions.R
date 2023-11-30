##Generate simulated data##

###### function to generate a bunch of simulations for a given volatility and EIR
### Random-walk on EIR but with a min =0.01 and max=1000
## TODO - check if any problem with only having a single volatility (i.e. draw from a distn rather than a set number)
generate_sim_compendium<-function(n_sims,volatility,init_EIR,duration,out_step){
  ## df with time-series of EIR,inc,prev and simmed volatility for each run
  sims_compendium<-data.frame(
    run=numeric(),t=numeric(),prev_true=numeric(),EIR_true=numeric(),vol_true=numeric(),inc_true=numeric()
  )
  ## run each sim and add it to df
  for(i in 1:n_sims){
    sims_compendium<-sims_compendium%>%add_row(
      cbind(
        data.frame(run=i),
        data_gen(EIR_volatility = volatility, init_EIR = init_EIR,time = duration,out_step = out_step)
      )
    )
  }
  ##return 'compendium' of sims
  return(sims_compendium)
}

## 
##
## I tried predicting the value at the end of a fitted window (hence 'lag' nomenclature that doesn't really work when predicting in the middle)
### this is the bit I think we can simplify/make less memory intensive using a rnn with Pete suggesting there might be a 'direction=both' argument that could help

#' function to generate torch data sets of a sliding window of observations through a set of runs with time index
#' I started by attempting to predict EIR at end of a window of lagged prevalence (hence 'lag' nomenclature)
#' but then added in option to place the predicted EIR in the middle of the window as prevalence trends before and after EIR provides information sooo
#'  sorry if lag nomenclature confusing
#' I think this can be simplifed/made less memory intensive using a rnn with Pete suggesting there might be a 'direction=both' argument that could help

#' @param sim_compendium set of simulations 
#' @param fit_var name of observed time series
#' @param pred_var name of time series wish to predict
#' @param t_var name of time variable 
#' @param lag_no number of 'lags' - see comment above about middle
#' @param middle whether we're predicting at the end of the window of observation or in the middle (middle)
#'
#' @return a tensor dataset with a  'fit_tensor' of all possible qualifying observed time-windows across all runs and 
#' a 'pred_tensor' being the value we wish to predict (i.e. either the value of the time series to predict at the end of the window or in the middle)
#' 
#' @export
#'
#' @examples
generate_torch_ds_lags<-function(sim_compendium,fit_var,pred_var,t_var,lag_no,middle=F){
  #create a 'fit_tensor' of all possible qualifying observed time-windows across all runs
  lags=1:lag_no
  map_lag <- lags %>% map(~partial(lag, n = .x))
  fit_tensor<-torch_tensor(
    as.matrix(
      ## arrange by run and time variable
      sim_compendium%>%group_by(run)%>%
        arrange(run,!!sym(t_var))%>%
        ### now add a window of lag_no observations going back in time from each timestep
        mutate(across(.cols = {{fit_var}}, .fns = map_lag, .names = "{.col}_lag{lags}"))%>%
        ### filter the earlier ones - would eventually fill these in based on equivalent assumption of previous dynamics (e.g. if in equilibrium) so can predict entire series
        filter(!!sym(t_var)>lag_no)%>%
        ## only select the windows to give the required tensor matrix
        ungroup()%>%
        select(contains("lag"))
    )
  )
  ### generate a tensor of values to be predicted 
  if(middle==T){
    ## pick the value in middle of observed window
    pred_tensor<-torch_tensor(
      as.matrix(sim_compendium%>%
                  ##arrange as before
                  group_by(run)%>%
                  arrange(run,!!sym(t_var))%>%
                  ### remove lag_no values but evenly from start and end (from memory think this means lag_no should be even) so that prediction is in the middle of the equivalent window
                  filter(!!sym(t_var)>floor(lag_no/2),!!sym(t_var)<(max(!!sym(t_var))-floor(lag_no/2)+1))%>%
                  ungroup()%>%
                  select(pred_var)))
  }
  else{
    ### same as above but just remove lag_no values from the start so that prediction is at end of relevant window
    pred_tensor<-torch_tensor(
      as.matrix(sim_compendium%>%
                  group_by(run)%>%
                  arrange(run,!!sym(t_var))%>%
                  filter(!!sym(t_var)>lag_no)%>%
                  ungroup()%>%
                  select(pred_var)))
  }
  ## return both the fit_tensor and pred_tensor as a tensor dataset
  return(tensor_dataset(fit_tensor, pred_tensor))
}


#' function that fits said model with sliding window structure with training, validation, testing dataset 

#'
#' @param model choice of model 
#' @param fit_var variable to fit
#' @param pred_var variable to predict
#' @param t_var time variable
#' @param middle whether or not 'predicting' in the middle of the window or at its end
#' @param sims_compendium_train set of sims for training
#' @param sims_compendium_test set of sims for testing
#' @param sims_compendium_valid set of sims for validation
#' @param d_hidden number of hidden states
#' @param epochs number of epochs
#' @param loss loss fn
#' @param optimizer optimiser fn
#' @param lag_no size of sliding window upon which to base prediction
#'
#' @return named list of the fitted model, predictions of the testing dataset and a default plot
#' @export
#'
#' @examples
generate_preds_valid_lag<-function(model,fit_var,pred_var,t_var,middle=F,sims_compendium_train,sims_compendium_test,sims_compendium_valid,
                                   d_hidden=100,epochs=50,loss=nn_mse_loss(),optimizer=optim_adam,lag_no=24){
  ### generate torch data sets for train,validate,test
  train_ds<-generate_torch_ds_lags(sim_compendium = sims_compendium_train,fit_var = fit_var,pred_var = pred_var,t_var=t_var,lag_no=lag_no,middle=middle)
  valid_ds<-generate_torch_ds_lags(sim_compendium = sims_compendium_valid,fit_var = fit_var,pred_var = pred_var,t_var=t_var,lag_no=lag_no,middle=middle)
  test_ds<-generate_torch_ds_lags(sim_compendium = sims_compendium_test,fit_var = fit_var,pred_var = pred_var,t_var=t_var,lag_no=lag_no,middle=middle)
  
  ## load datasets 
  train_dl<-dataloader(train_ds, batch_size = 100, shuffle = TRUE)
  valid_dl<-dataloader(valid_ds, batch_size = 100, shuffle = TRUE)
  test_dl<-dataloader(test_ds, batch_size = 100, shuffle = TRUE)
  
  # output dimensionality (number of predicted features - here just one)
  d_out <- ncol(train_ds$tensors[[2]])
  ## number of variables to use for prediction (here a window of observation of the data lag_no wide)
  d_in <-  ncol(train_ds$tensors[[1]])
  ## number of datapoints (number of windows*no of runs)
  n<-nrow(train_ds$tensors[[1]])

  ### fit the model
  fitted <- model %>%
    setup(loss = loss, 
          optimizer = optimizer,
          metrics = list(luz_metric_mae())) %>%
    set_hparams(
      d_in = d_in,
      d_hidden = d_hidden, d_out = d_out
    ) %>%
    fit(train_dl, epochs = epochs, valid_data = valid_dl)
  
  ### predict the test 
  pred<-fitted %>% predict(test_ds)
  if(middle==T){
    predict_data_df<-sims_compendium_test%>%group_by(run)%>%
      arrange(run,!!sym(t_var))%>%
      filter(!!sym(t_var)>floor(lag_no/2),!!sym(t_var)<(max(!!sym(t_var))-floor(lag_no/2)+1))%>%
      select(run,t,pred_var)
  }
  else{
    predict_data_df<-sims_compendium_test%>%group_by(run)%>%
      group_by(run)%>%
      arrange(run,!!sym(t_var))%>%
      filter(!!sym(t_var)>lag_no)%>%
      ungroup()%>%
      select(run,t,pred_var)
    
  }
  predict_data_df$predictions=as.vector(t(as.matrix(pred)))
  compare_plot<-ggplot(predict_data_df,aes(x=t,y=!!sym(pred_var)))+
    geom_line()+geom_line(aes(y=predictions),col="red")+
    facet_wrap(~run)+theme(strip.text.x = element_blank())
  return(list(fit_model=fitted,
              predictions=predict_data_df,
              compare_plot=compare_plot)
  )
}
