N = 10000
U = runif(N)
x = rep(0,N)
for(i in 1:N){
  if(U[i]<=0.5){
    x[i] = rnorm(1)
  }else{
    x[i] = 0
  }
}


loss_function = function(par,x,w){
  pi = par[1]
  tau = par[2]
  f_vec = c(mean(x^2)-pi*tau,mean(x^4)-3*pi*tau^2)
  loss = t(f_vec)%*%w%*%f_vec
  return(loss)
}



deriv_loss = function(par,x,w){
  pi = par[1]
  tau = par[2]
  G_mat = matrix(c(-tau,-3*tau^2,-pi,-6*pi*tau),2,2)
  f_vec = c(mean(x^2)-pi*tau,mean(x^4)-3*pi*tau^2)
  
  deriv = 2*t(G_mat)%*%w%*%f_vec
  return(deriv)
}



#step 1
w = diag(2)
lower = c(0,0)
upper = c(1,5)
result = optim(c(0.1,2),fn = loss_function,gr =deriv_loss ,method = 
        "L-BFGS-B",lower= lower,upper = upper,x=x,w=w)
coef_est = result$par

w = solve(var_f(coef_est,x))
result = optim(c(0.1,2),fn = loss_function,gr =deriv_loss ,method = 
                 "L-BFGS-B",lower= lower,upper = upper,x=x,w=w)
coef_est = result$par




var_f = function(par,x){
  pi = par[1]
  tau = par[2]
  vec1 = x^2-pi*tau
  vec2 = x^4-3*pi*tau^2
  result = cov(cbind(vec1,vec2))
  return(result)
}

