#include <TMB.hpp>

template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}


template <class Type>
Type minus_one_to_one(Type x)
{
  return Type(2) * invlogit(x) - Type(1);
}

 // dlnorm
template<class Type>
Type dlnorm(Type x, Type meanlog, Type sdlog, int give_log=0){
  Type logres;
  logres = dnorm( log(x), meanlog, sdlog, true) - log(x);
  //return 1/(sqrt(2*M_PI)*sd)*exp(-.5*pow((x-mean)/sd,2));
  if(give_log)return logres; else return exp(logres);
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(obs_logR);   // observed log recruitment
  DATA_VECTOR(obs_S);    // observed  Spawner
  
  PARAMETER(alpha);
  PARAMETER(logbeta);
  PARAMETER(logSigObs);
  PARAMETER(rho);
  //PARAMETER_VECTOR(delta);
  
  int timeSteps=obs_logR.size();

  Type rhoo = minus_one_to_one(rho);

  //Type rho = exp(logrho);
  Type beta = exp(logbeta);
  Type SigObs = exp(logSigObs);
  Type Smax  = Type(1.0)/beta;
  
  Type tau  = Type(1.0)/(SigObs*SigObs);

  Type SigAR  = SigObs*sqrt(1-pow(rhoo,2));
  
  vector<Type> pred_logR(timeSteps), logRS(timeSteps), residuals(timeSteps), epsilon(timeSteps);

 

  //priors on precision and variance ratio
  //Type ans= -dbeta(rho,Type(3.0),Type(3.0),true);  
  //Type ans= -dnorm(logSigObs,Type(0.0),Type(5.0),true);   
  Type ans = Type(0);

 
  logRS(0) = alpha - beta * obs_S(0) ;
  pred_logR(0) = logRS(0) + log(obs_S(0)) ; 
  residuals(0) = obs_logR(0) - pred_logR(0);
  epsilon(0) =  residuals(0);//delta(0);
  
  //ans+= -dnorm(delta(0), Type(0),SigObs,true);  
  ans+= -dnorm(obs_logR(0),pred_logR(0),SigObs,true);

  for(int i=1;i<timeSteps;i++){
    if(!isNA(obs_logR(i))){
      logRS(i) = alpha - beta * obs_S(i) ;
      //pred_logR(i) = logRS(i) + log(obs_S(i)) + epsilon(i) ;
      pred_logR(i) = logRS(i) + log(obs_S(i)) + epsilon(i-1) * rhoo ;
      residuals(i) = obs_logR(i) - pred_logR(i);
      epsilon(i) = residuals(i);//epsilon(i-1) * rho; //+ delta(i)* sqrt(1-pow(rho,2));
      
      ans+=-dnorm(obs_logR(i),pred_logR(i),SigAR,true);
      //ans+=-dnorm(obs_logR(i),pred_logR(i),SigAR,true);
      //ans+=-dnorm(delta(i),Type(0),SigObs,true);
      //
    }
  
  }
  
  Type umsy     = Type(.5) * alpha - Type(0.07) * (alpha * alpha);
  Type Smsy     = alpha/beta * (Type(0.5) -Type(0.07) * alpha);  

  REPORT(pred_logR)
  REPORT(residuals)
  REPORT(alpha)
  //REPORT(delta)
  REPORT(tau)
  REPORT(rhoo)
  REPORT(beta)
  REPORT(SigObs)
  REPORT(SigAR)
  REPORT(Smax)
  REPORT(umsy)
  REPORT(epsilon)
  return ans;
}

