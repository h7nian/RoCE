# Derivative of corrected site scores with projection coefficients held fixed.
.source_score_derivative_rows <- function(state, block, site, coefficients) {
  stopifnot(site %in% c("target","source"), length(coefficients)==length(state$theta))
  index<-state$indices; theta<-state$theta
  W<-block$W; Z<-block$Z; A<-block$A; Y<-block$Y
  rows<-matrix(0,nrow(W),length(theta))
  alpha_final<-theta[index$alpha_final]; gamma_final<-theta[index$gamma_final]
  final_mean<-plogis(drop(W%*%alpha_final)); final_log_weight<-drop(Z%*%gamma_final)
  final_alpha_action<-drop(W%*%coefficients[index$alpha_final])
  final_gamma_action<-drop(Z%*%coefficients[index$gamma_final])
  source_weight_mean<-source_derivative_mean<-numeric(nrow(W))
  clip<-function(x,bound)pmax(-bound,pmin(bound,x))
  for(key in state$keys) {
    alpha_index<-index[[paste0("alpha_initial_",key)]]
    gamma_index<-index[[paste0("gamma_initial_",key)]]
    initial_predictor<-drop(W%*%theta[alpha_index])
    initial_mean<-plogis(initial_predictor)
    if(site=="target") {
      derivative<-initial_mean*(1-initial_mean)
      rows[,alpha_index]<-W*(A*derivative*drop(W%*%coefficients[alpha_index])-
        derivative*(1-2*initial_mean)*final_gamma_action/length(state$keys))
    } else {
      fraction<-nrow(state$pieces[[key]]$source_calibration$W)/nrow(state$source$W)
      initial_log_weight<-drop(Z%*%theta[gamma_index])
      initial_weight<-exp(-clip(initial_log_weight,state$M_fit))
      clipped_mean<-plogis(clip(initial_predictor,state$M_fit))
      derivative<-clipped_mean*(1-clipped_mean)
      rows[,alpha_index]<-W*(fraction*A*exp(-final_log_weight)*derivative*(1-2*clipped_mean)*
        (abs(initial_predictor)<state$M_fit)*final_gamma_action)
      rows[,gamma_index]<-Z*(-A*exp(-initial_log_weight)*drop(Z%*%coefficients[gamma_index])+
        fraction*A*initial_weight*(Y-final_mean)*(abs(initial_log_weight)<state$M_fit)*final_alpha_action)
      source_weight_mean<-source_weight_mean+fraction*initial_weight
      source_derivative_mean<-source_derivative_mean+fraction*derivative
    }
  }
  if(site=="target") rows[,index$alpha_final]<-W*(final_mean*(1-final_mean)) else {
    inference_weight<-exp(-clip(final_log_weight,state$M_inference))
    rows[,index$alpha_final]<-W*(A*final_mean*(1-final_mean)*
      (-inference_weight+source_weight_mean*final_alpha_action))
    rows[,index$gamma_final]<-Z*(-A*inference_weight*(Y-final_mean)*
      (abs(final_log_weight)<state$M_inference)-A*exp(-final_log_weight)*
      source_derivative_mean*final_gamma_action)
  }
  stopifnot(all(is.finite(rows)))
  rows
}
