library(ggplot2)

simulate_once =function(Sym1){

n = nrow(Sym1)
K = 2
# Trimming
A_gam_indeg = apply(Sym1,2,sum)
trim_cond = A_gam_indeg<5*K*sum(Sym1)/nrow(Sym1)
A_gam = Sym1[trim_cond,trim_cond]



#Spectral decomp

s = svd(A_gam)
D = diag(s$d)
diag(D) = c(diag(D)[1:K],rep(0.0,nrow(D)-K))
A_hat = s$u%*% D%*%t(s$v)

Q_ivret = function(i,v,A){
  tmp = sapply(1:nrow(A),function(x)norm(as.matrix(A[x,])-as.matrix(A[v,]))^2)
  ret = which(tmp<=i*sum(A)/(100*nrow(A)^2))
  return(ret)
}
r_list = numeric(as.integer(log(n)))
T_list = list()
for(i in 1:as.integer(log(n))){
  Q_iv = sapply(1:nrow(A_hat),function(x)Q_ivret(i,x,A_hat))
  T_i = list()
  eps_ik = list()
  for (k in 1:K){
    if (k ==1){
    v_kstar = which.max(sapply(1:length(Q_iv),function(x) length(Q_iv[x][[1]])))
    T_i[[k]] = Q_iv[v_kstar][[1]]
    }else{
    v_kstar = which.max(sapply(1:length(Q_iv),
                     function(x) tryCatch(length(setdiff(Q_iv[x][[1]],
                                                Reduce(function(l,r) union(l,r),T_i))[[1]]), error=function(err) 0)))
    T_i[[k]] = setdiff(Q_iv[v_kstar][[1]],Reduce(function(x,y) union(x,y),T_i))[[1]]
    }
    eps_ik[[k]] = colSums(matrix(A_hat[T_i[[k]],]/length(T_i[[k]]),
                                    nrow = length(T_i[[k]]),ncol = ncol(A_hat)))
  }
  v_rem = setdiff(1:nrow(A_hat),Reduce(function(x,y) union(x,y),T_i[1:K]))
  for(v in v_rem){
    kstar = which.min(sapply(1:K,function(x) norm(as.matrix(A_hat[v,])-as.matrix(eps_ik[[x]]))))
    T_i[[kstar]] = union(T_i[[kstar]],v)
  }
  r_i = sum(sapply(1:K,function(x)sum(sapply(T_i[[x]],function(y)norm(as.matrix(A_hat[y,])-as.matrix(eps_ik[[x]]))^2))))
  r_list[i] = r_i
  T_list[[i]] = T_i
}
i_star = which.max(r_list)
S_k = list(T_list[[i_star]][[1]],T_list[[i_star]][[2]])

S_ik = list()
S_ik[[1]] = S_k
for (i in 1:as.integer(log(n))){
  first = TRUE
  for(v in 1:nrow(Sym1)){
    kstar = which.max(sapply(1:K,function(x) sum(Sym1[v,S_ik[[i]][[x]]])/length(S_ik[[i]][[x]])))
    if(first){
      S_ik[[i+1]] = list()
      S_ik[[i+1]][[kstar]] = v
      first = FALSE
    }else{
      S_ik[[i+1]][[kstar]] = tryCatch(union(S_ik[[i+1]][[kstar]],v), error=function(err) v)
    }
    
  }
}

N = length(S_ik)

return(S_ik[[N]])
}

num_nodes = 100
acc1_total = numeric(100)
accm_total = numeric(100)
accm2_total = numeric(100)
#one network
gen_network =function(num,across,within){
  Sym1 = sapply(1:num,
                function(x) sapply(1:num,
                                   function(y) ifelse(y<x,0,
                                                      ifelse(x<(num/2+1)&&y>(num/2),
                                                             rbinom(1,1,across),
                                                             rbinom(1,1,within)))))
  Sym1 = Sym1+t(Sym1)
  diag(Sym1) = 0
  return(Sym1)
}
true = c(rep(1,num/2),rep(2,num/2))
tics = seq(0.01,0.08,0.005)
first=TRUE
for(tic in tics){
  
  for(i in 1:100){
    Sk = simulate_once(gen_network(num_nodes,tic,0.1))
    acc1 = c((sum(1:(num_nodes/2) %in% Sk[[1]]) + sum((num_nodes/2+1):num_nodes %in% Sk[[2]]))/num_nodes,
             (sum(1:(num_nodes/2) %in% Sk[[2]]) + sum((num_nodes/2+1):num_nodes %in% Sk[[1]]))/num_nodes)
    acc1_total[i] = max(acc1)
  }
  
  for(i in 1:100){
    Sk = simulate_once(Reduce("+",lapply(1:3,function(x)gen_network(num_nodes,tic,0.1))))
    acc1 = c((sum(1:(num_nodes/2) %in% Sk[[1]]) + sum((num_nodes/2+1):num_nodes %in% Sk[[2]]))/num_nodes,
             (sum(1:(num_nodes/2) %in% Sk[[2]]) + sum((num_nodes/2+1):num_nodes %in% Sk[[1]]))/num_nodes)
    accm2_total[i] = max(acc1)
  }
  if(first){
    df = cbind(mean(acc1_total),mean(accm2_total))
    first = FALSE
  }else{
  df = rbind(df,cbind(mean(acc1_total),mean(accm2_total)))
  }
}  

sim_data = data.frame(var0 = df[,1],var1 = df[,2],across_prob = seq(0.01,0.08,0.005))
colnames(sim_data) = c("one_network","three_networks","across_prob")
ggplot(sim_data, aes(across_prob)) + 
  geom_line(aes(y = one_network, colour = "One Network")) + 
  geom_line(aes(y = three_networks, colour = "Three Networks"))+ 
  ylab('Recovery Accuracy')+xlab('Across Group Connection Prob (within is 0.1)')


for(i in 1:100){
  votes = numeric(num_nodes)
  for(j in 1:3){
    Sk = simulate_once(Sym1,num_nodes,0.05,0.1)
    acc1 = c((sum(1:(num_nodes/2) %in% Sk[[1]]) + sum((num_nodes/2+1):num_nodes %in% Sk[[2]]))/num_nodes,
             (sum(1:(num_nodes/2) %in% Sk[[2]]) + sum((num_nodes/2+1):num_nodes %in% Sk[[1]]))/num_nodes)
    rev = ifelse(which.max(acc1)==1,1,-1)
    votes = votes+sapply(1:num_nodes,function(x)ifelse(x %in%Sk[[1]],1*rev,-1*rev))
  }
  accm = c((sum(votes[1:(num_nodes/2)]>0)+sum(votes[(num_nodes/2+1):num_nodes]<0))/num_nodes,
           (sum(votes[1:(num_nodes/2)]<0)+sum(votes[(num_nodes/2+1):num_nodes]>0))/num_nodes)
  accm_total[i] = max(accm)
}