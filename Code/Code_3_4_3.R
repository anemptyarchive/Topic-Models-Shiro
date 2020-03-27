
# Stochastic Variational Bayes for LDA ------------------------------------

# 利用パッケージ
library(tidyverse)



# 確率的変分ベイズ推定 --------------------------------------------------------------


for(s in 1:S) {
  
  for(II in 1:InIter) {
    
    for(i in 1:n_d) {

      
    }
  }
  
  for(k in 1:K) { ## (各トピック)
    
    # :(3.159)
    xi_kv[k, ] <- xi_kv[k, ] + nu * (M * n_dkv[, k, ] + beta_v - xi_kv[k, ])
    
  } ## (/各トピック)
}
