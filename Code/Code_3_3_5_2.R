
# Variational Bayes for LDA (2) -------------------------------------------

# 利用パッケージ
library(tidyverse)


# パラメータの設定 ----------------------------------------------------------------

# イタレーション数
OutIter <- 5
InIter  <- 5
Abs_Err <- 0.1

# トピック数
K <- 5

# 事前分布のパラメータ
alpha_k <- rep(2, K)
beta_v  <- rep(2, V)

# 潜在トピック集合の分布
z_dvk <- array(0, dim = c(M, V, K))
for(d in 1:M) {
  for(v in 1:V) {
    tmp_q_z <- seq(0, 1, by = 0.01) %>% 
      sample(size = K, replace = TRUE)
    z_dvk[d, v, ] <- tmp_q_z / sum(tmp_q_z)
  }
}

## カウントの期待値
tmp_z <- array(0, dim = c(M, V, K))
for(k in 1:K) {
  tmp_z[, , k] <- z_dvk[, , k] * n_dv
}

# 文書ごとにおいて各トピックが割り当てられた単語数
n_dk <- apply(tmp_z, 1, sum) / K %>% 
        matrix(nrow = M, ncol = K)

# 全文書において各トピックが割り当てられた単語数
n_kv <- apply(tmp_z, c(3, 2), sum)

# 処理の検証用
sum(n_dk) == sum(n_dv)
sum(n_kv) == sum(n_dv)


# 事後分布パラメータの初期値(=事前分布のパラメータ)
eta_dk <- t(t(n_dk) + alpha_k)
eta_kv <- t(t(n_kv) + beta_v)



# 変分ベイズ -------------------------------------------------------------------

# 推移の確認用のデータフレームを作成
trace_eta_theta <- cbind(
  as.data.frame(eta_dk), 
  doc = as.factor(1:M),  # 文書番号
  Iter = 0 # 試行回数
)
trace_eta_phi <- cbind(
  as.data.frame(eta_kv), 
  topic = as.factor(1:K),  # トピック番号
  Iter = 0 # 試行回数
)

for(OI in 1:OutIter) { ## 
  
  for(d in 1:M) { ## (各文書)
    
    # 絶対誤差を計算
    if(OI > 1) { # 初回を飛ばす
      abs_err <- sum(abs(n_dk.new - n_dk) / n_d) / K
    } else {
      abs_err <- Abs_Err + 1
    }
    # インナーループ回数を初期化
    II <- 1
    while(abs_err >= Abs_Err && II <= InIter) { ## 
      
      # インナーループの回数を更新
      II <- II + 1
      
      for(v in 1:V) { ## (各語彙)
        if(n_dv[d, v] > 0) {
          for(n in 1:n_dv[d, v]) { ## (各単語)
            
            # 潜在トピック集合の事後分布を計算:式(3.99)
            term1 <- digamma(eta_kv[, v]) - digamma(apply(eta_kv, 1, sum))
            term2 <- digamma(eta_dk[d, ]) - digamma(apply(eta_dk, 2, sum))
            tmp_z_dvk <- exp(term1) * exp(term2)
            # 正規化
            z_dvk[d, v, ] <- tmp_z_dvk / sum(tmp_z_dvk)
            
          } ## (/各単語)
        }
      } ## (/各語彙)
      
      # カウントの期待値
      for(k in 1:K) {
        tmp_z[, , k] <- z_dvk[, , k] * n_dv
      }
      n_dk.new <- apply(tmp_z, 1, sum) / K %>% 
                  matrix(nrow = M, ncol = K)
      n_kv <- apply(tmp_z, c(3, 2), sum)

      # 事後分布のパラメータを計算:式(3.89)
      eta_dk[d, ] <- n_dk[d, ] + alpha_k
      
      print(paste0("OutIter=", OI, ", InIter=", II, ", Abs_Err=", abs_err))
    }
  } ## (/各文書)
  
  for(k in 1:K) { ## (各トピック)
    
    # 事後分布のパラメータ:式(3.95)
    eta_kv[k, ] <- n_kv[k, ] + beta_v
    
  } ## (/各トピック)
  
  # 推移の確認用のデータフレームを作成
  tmp_trace_eta_theta <- cbind(
    as.data.frame(eta_dk), 
    doc = as.factor(1:M),  # 文書番号
    Iter = I # 試行回数
  )
  tmp_trace_eta_phi <- cbind(
    as.data.frame(eta_kv), 
    topic = as.factor(1:K),  # トピック番号
    Iter = I # 試行回数
  )
  
  # データフレームを結合
  trace_eta_theta <- rbind(trace_eta_theta, tmp_trace_eta_theta)
  trace_eta_phi   <- rbind(trace_eta_phi, tmp_trace_eta_phi)
}


