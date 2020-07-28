#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(S);
  DATA_VECTOR(logR);
  DATA_IVECTOR(stk);
  DATA_IVECTOR(flag);
  DATA_IVECTOR(cycle);
  DATA_INTEGER(N_Stks);
  
  PARAMETER_VECTOR(logA);
  PARAMETER_MATRIX(logB);
  PARAMETER_VECTOR(logSigma);
 // PARAMETER_MATRIX(logSgen);

  
  Type ans=0.0;
  vector <Type> sigma=exp(logSigma);
  int N_Obs = S.size(); 
  vector <Type> LogR_Pred(N_Obs);
  vector <Type> A = exp(logA);
  matrix <Type> B = exp(logB.array());

   matrix <Type> SMSY(N_Stks, 4);  //problem here maybe? since each year had different value?
   matrix <Type> LogSMSY(N_Stks, 4);
  // matrix <Type> Sgen = exp(logSgen.array());
  
  // for each stock need to start in year 4
  for(int i=0; i<N_Obs; i++){
    if(flag(i) == 1){
      
      LogR_Pred(i) = logA(stk(i)) + log(S(i)) - B(stk(i), 0)*S(i) - B(stk(i), 1)*S(i-1) - 
                      B(stk(i), 2)*S(i-2) - B(stk(i), 3)*S(i-3) ;
      ans += -dnorm(LogR_Pred(i), logR(i),  sigma(stk(i)), true);
    }
  }
  
  
  // for(int i=0; i<N_Obs; i++){
  //   if(flag(i) > 0  ){
      //SMSY(stk(i), cycle(i)) = (logA(stk(i)) - B(stk(i), 1)*S(i-1) - B(stk(i), 2)*S(i-2) - B(stk(i), 3)*S(i-3) )*
                              // (0.5-0.07*(logA(stk(i)) - B(stk(i), 1)*S(i-1) - B(stk(i), 2)*S(i-2) - 
                              // B(stk(i), 3)*S(i-3)))/ B(stk(i), 0);
      // LogSMSY(stk(i), cycle(i)) = (logA(stk(i)) - B(stk(i), 1)*S(i-1) - B(stk(i), 2)*S(i-2) - B(stk(i), 3)*S(i-3) ) +
      //                          logSgen(stk(i), cycle(i)) - B(stk(i), 0) * Sgen(stk(i), cycle(i));
      // Type Diff = LogSMSY(stk(i), cycle(i)) - log(SMSY(stk(i), cycle(i)));
      // ans += -dnorm(Diff, Type(0), Type(1));
  //   }
  // }

  
  //ADREPORT(SMSY);
  ADREPORT(A);
  //ADREPORT(Sgen);
  ADREPORT(sigma);
  ADREPORT(B);

  
  return ans;
  
}
