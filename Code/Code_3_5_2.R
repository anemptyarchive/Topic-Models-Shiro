
# Particle filter for LDA -------------------------------------------------

# 利用パッケージ
library(tidyverse)


# パラメータの設定 ----------------------------------------------------------------

# サンプリング数
S <- 10

# リサンプリング数
R <- 30

# 閾値
threshold <- 20

# トピック数
K <- 5

# 事前分布のパラメータ
alpha_k <- rep(2, K)
beta_v <- rep(2, V)

# 潜在トピック集合の事後分布
z_di_k <- array(0, dim = c(M, V, max(n_dv), K))
for(d in 1:M) {
  for(v in 1:V) {
    if(n_dv[d, v] > 0) {
      for(n in 1:n_dv[d, v]) {
        # ランダムに値を生成
        tmp_q_z <- sample(seq(0, 1, by = 0.01), size = K, replace = TRUE)
        # 正規化
        z_di_k[d, v, n, ] <- tmp_q_z / sum(tmp_q_z)
      }
    }
  }
}

# 潜在トピック集合の受け皿
z_di_s <- array(0, dim = c(M, V, max(n_dv), S))

# 重みの受け皿
w_z_di_s <- array(1/S, dim = c(M, V, max(n_dv), S))

# 割り当てられたトピックに関する単語数の初期値
n_dk_s <- array(0, dim = c(M, K, S))
n_kv_s <- array(0, dim = c(K, V, S))


# パーティクルフィルタ --------------------------------------------------------------

for(d in 1:M) { ## (各文書)
  
  # 1期(語)前の値i-1を初期化
  old_v <- old_n <- 1 # (本当は0)

  for(v in 1:V) { ## (各語彙)
    if(n_dv[d, v] > 0) {
      for(n in 1:n_dv[d, v]) { ## (各単語)
        
        for(s in 1:S) { ## (サンプリング)
          
          # 潜在トピック集合の分布を計算:式(3.172)
          term1 <- (n_kv_s[, v, s] + beta_v[v]) / apply(t(n_kv_s[, , s]) + beta_v, 2, sum)
          term2 <- (n_dk_s[d, , s] + alpha_k) / sum(n_dk_s[d, , s] + alpha_k)
          tmp_q_z <- term1 * term2
          z_di_k[d, v, n, ] <- tmp_q_z / sum(tmp_q_z) # 正規化
          
          # 潜在トピックをサンプリング
          z_di_s[d, v, n, s] <- sample(1:K, size = 1, prob = z_di_k[d, v, n, ])
          
          ## 重みを計算:式(3.175)
          # 式(3.176)の計算
          term1 <- (n_kv_s[, v, s] + beta_v[v]) / apply(t(n_kv_s[, , s]) + beta_v, 2, sum)
          term2 <- (n_dk_s[d, , s] + alpha_k) / sum(n_dk_s[d, , s] + alpha_k)
          p_wz_s <- term1 * term2
          
          # 重みを計算:式(3.175)
          w_z_di_s[d, v, n, s] <- w_z_di_s[d, old_v, old_n, s] * sum(p_wz_s)
          
          # 割り当てられた潜在トピックに応じてカウントに加算
          tmp_n_k <- rep(0, K) # 初期化
          tmp_n_k[z_di_s[d, v, n, s]] <- 1 # k(=z_di)番目に1を代入
          n_dk_s[d, , s] <- n_dk_s[d, , s] + tmp_n_k
          n_kv_s[, v, s] <- n_kv_s[, v, s] + tmp_n_k
          
        } ## (/サンプリング)
        
        # 重みを正規化
        w_z_di_s[d, v, n, ] <- w_z_di_s[d, v, n, ] / sum(w_z_di_s[d, v, n, ])
        
        # ESSを計算
        ESS <- 1 / sum(w_z_di_s[d, v, n, ] ^ 2)
        
        if(ESS < threshold) {
          
          for(r in 1:R) { ## (活性化サンプル)
            
            # 活性化する重みをサンプリング
            l <- sample(1:M, size = 1)
            m_v <- sample(1:V, size = 1, prob = n_dv[l, ])
            m_n <- sample(1:n_dv[l, m_v], size = 1)
            
            for(s in 1:S) { ## (リサンプリング)
              
              # カウントを初期化
              n_dk.lm_s <- rep(0, K)
              n_kv.lm_s <- matrix(0, nrow = K, ncol = V)
              
              # l(d)=1のときに配列の構造が変わる対策(仮)
              if(l == 1) {
                mar <- 1
              } else if(l > 1) {
                mar <- 2
              }
              
              # 潜在トピックに関するカウントを計算
              for(k in 1:K) {
                n_dk.lm_s[k] <- sum(z_di_s[1:l, 1:m_v, 1:m_n, s] == k)
                n_kv.lm_s[k, ] <- apply(z_di_s[1:l, , , s] == k, mar, sum)
              }
              
              # l,m要素について取り除く(iをvとnに分けているためi-1の処理がめんどいため)
              tmp_n_k <- rep(0, K) # 初期化
              tmp_n_k[z_di_s[l, m_v, m_n, s]] <- 1 # k番目に1を代入
              n_dk.lm_s <- n_dk.lm_s - tmp_n_k
              n_kv.lm_s <- n_kv.lm_s - tmp_n_k
              
              # リサンプリング確率を計算:式(3.179)
              term1 <- (n_kv.lm_s[, m_v] + beta_v[m_v]) / apply(t(n_kv.lm_s) + beta_v, 2, sum)
              term2 <- (n_dk.lm_s + alpha_k) / sum(n_dk.lm_s + alpha_k)
              tmp_q_z <- term1 * term2
              
              # 潜在トピックをサンプリング
              z_di_s[l, m_v, m_n, s] <- sample(1:K, size = 1, prob = tmp_q_z)
            
            } ## (/リサンプリング)
          } ## (/活性化サンプル)
          
          # 重みを初期化
          w_z_di_s <- array(1 / S, dim = c(M, V, max(n_dv), S))
        }
        
        # 1期(語)前の添字を保存
        old_v <- v
        old_n <- n
        
      } ## (/各単語)
    }
  } ## (/各語彙)
} ## (/各文書)

warnings()

# 推定結果の確認 -----------------------------------------------------------------

w_z_di_k <- array(0, dim = c(M, V, max(n_dv), K))
for(d in 1:M) {
  for(v in 1:V) {
    if(n_dv[d, v] > 0) {
      for(n in 1:n_dv[d, v]) {
        for(k in 1:K) {
          tmp_w_s <- w_z_di_s[d, v, n, ]
          tmp_z_s <- z_di_s[d, v, n, ]
          w_z_di_k[d, v, n, k] <- sum(tmp_w_s[tmp_z_s == k])
        }
      }
    }
  }
}
w_z_dv_k <- apply(w_z_di_k, c(1, 2, 4), sum)

w_WideDF <- data.frame()
for(k in 1:K) {
  tmp_df <- cbind(
    as.data.frame(w_z_dv_k[, , k]), 
    doc = as.factor(1:M), 
    topic = as.factor(k)
  )
  w_WideDF <- rbind(w_WideDF, tmp_df)
}

w_LongDF <- pivot_longer(
  w_WideDF, 
  cols = -c(doc, topic), 
  names_to = "word", 
  names_prefix = "V", 
  names_ptypes = list(word = factor()), 
  values_to = "value"
)

ggplot(w_LongDF, aes(x = word, y = value, fill = word, color = word)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  facet_wrap(~ topic, labeller = label_both) + 
  theme(legend.position = "none") + 
  scale_x_discrete(breaks = seq(1, V, by = 10)) + 
  labs(title = "Particle filter for LDA", 
       subtitle = "w")




n_dvk <- array(0, dim = c(M, V, K))
for(k in 1:K) {
  n_dvk[, , k] <- apply(z_di_s == k, c(1, 2), sum)
}

n_kv <- matrix(0, nrow = K, ncol = V)
for(k in 1:K) {
  n_kv[k, ] <- apply(z_di_s == k, 2, sum)
}

p_n_kv <- t(t(n_kv) / apply(n_kv, 2, sum))


n_kv_WideDF <- cbind(
  as.data.frame(t(n_kv)), 
  word = v_index[["TERM"]]
)

n_kv_LongDF <- pivot_longer(
  n_kv_WideDF, 
  cols = -word, 
  names_to = "topic", 
  names_prefix = "V", 
  names_ptypes = list(topic = factor()), 
  values_to = "value"
)

ggplot(n_kv_LongDF, aes(x = topic, y = value, fill = topic)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  facet_wrap(~ word) + 
  labs(title = "Particle filter for LDA")

ggplot(n_kv_LongDF, aes(x = word, y = value, fill = word, color = word)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  facet_wrap(~ topic, labeller = label_both) + 
  theme(legend.position = "none", axis.text.x = element_text(angle = 90)) + 
  labs(title = "Particle filter for LDA")

# try ---------------------------------------------------------------------

z_di_k <- seq(0, 1, by = 0.01) %>% 
  sample(size = M * V * max(n_dv) * S, replace = TRUE) %>% 
  array(c(M, V, max(n_dv), K))
for(k in 1:K) {
  z_di_k[, , , k] <- z_di_k[, , , k] / apply(z_di_k, c(1, 2, 3), sum)
  #z_di_k[, , , k][n_dv == 0] <- 0
}
array(rep(2, 8), dim = c(2, 2, 2)) - array(rep(1, 8), dim = c(2, 2, 2))

