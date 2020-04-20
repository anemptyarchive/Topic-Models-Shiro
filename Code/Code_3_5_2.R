
# Particle filter for LDA -------------------------------------------------

# 利用パッケージ
library(tidyverse)


# パラメータの設定 ----------------------------------------------------------------

# サンプリング数(粒子の数)
S <- 50

# リサンプリング数
R <- 10

# 閾値
threshold <- 10

# トピック数
K <- 5

# 事前分布のパラメータ
alpha_k <- rep(2, K)
beta_v <- rep(2, V)

# 潜在トピック集合の事後分布
z_dv_k <- array(0, dim = c(M, V, K))
for(d in 1:M) {
  for(v in 1:V) {
    if(n_dv[d, v] > 0) {
      # ランダムに値を生成
      tmp_q_z <- sample(seq(0, 1, by = 0.01), size = K, replace = TRUE)
      # 正規化
      z_dv_k[d, v, ] <- tmp_q_z / sum(tmp_q_z)
    }
  }
}

# 潜在トピック集合の受け皿
z_di_s <- array(0, dim = c(M, V, max(n_dv), S))

# 重みの受け皿
w_z_di_s <- array(1/S, dim = c(M, V, max(n_dv), S))

# 割り当てられたトピックに関する単語数の受け皿
n_dk <- rep(0, K) # (ベクトルで扱うことに注意)
n_kv <- matrix(0, nrow = K, ncol = V)


# パーティクルフィルタ --------------------------------------------------------------

for(d in 1:M) { ## (各文書)
  
  # 動作確認用
  start_time <- Sys.time()
  
  # 文書dにおいて各トピックが割り当て単語数
  n_dk <- rep(0, K) # (ここではベクトルで扱うことに注意)
  
  # 1期(語)前の値i-1を初期化
  old_v <- old_n <- 1 # (本当は0)

  for(v in 1:V) { ## (各語彙)
    if(n_dv[d, v] > 0) {
      for(n in 1:n_dv[d, v]) { ## (各単語)
        
        for(s in 1:S) { ## (サンプリング)
          
          # 潜在トピック集合の分布を計算:式(3.172)
          term1 <- (n_kv[, v] + beta_v[v]) / apply(t(n_kv) + beta_v, 2, sum)
          term2 <- (n_dk + alpha_k) / sum(n_dk + alpha_k)
          tmp_q_z <- term1 * term2
          z_dv_k[d, v, ] <- tmp_q_z / sum(tmp_q_z) # 正規化
          
          # 潜在トピックをサンプリング
          k <- sample(1:K, size = 1, prob = z_dv_k[d, v, ])
          z_di_s[d, v, n, s] <- k
          
          # 割り当てられたトピックに関する単語数に加算
          n_dk[k] <- n_dk[k] + 1
          n_kv[k, v] <- n_kv[k, v] + 1
          
          ## 重みを計算:式(3.175),(3.176)
          term1 <- (n_kv[, v] + beta_v[v]) / apply(t(n_kv) + beta_v, 2, sum)
          term2 <- (n_dk + alpha_k) / sum(n_dk + alpha_k)
          w_z_di_s[d, v, n, s] <- w_z_di_s[d, old_v, old_n, s] * sum(term1 * term2)
          
        } ## (/サンプリング)
        
        # 重みを正規化
        w_z_di_s[d, v, n, ] <- w_z_di_s[d, v, n, ] / sum(w_z_di_s[d, v, n, ])
        
        # ESSを計算
        ESS <- 1 / sum(w_z_di_s[d, v, n, ] ^ 2)
        
        if(ESS < threshold) {
          
          for(r in 1:R) { ## (活性化サンプル)
            
            # 活性化する重みをサンプリング(仮)
            vec_n_dv <- as.vector(t(n_dv))
            re_w <- sample(
              1:((d - 1) * V + v), 
              size = 1, 
              prob = vec_n_dv[1:((d - 1) * V + v)] # 出現回数を確率として使用
            )
            l <- re_w %/% V + 1
            m_v <- re_w %% V
            m_n <- sample(1:n_dv[l, m_v], size = 1)
            
            for(s in 1:S) { ## (リサンプリング)
              
              # カウントを初期化
              n_dk.lm <- rep(0, K)
              n_kv.lm <- matrix(0, nrow = K, ncol = V)
              
              # 潜在トピックに関するカウントを計算
              for(k in 1:K) {
                n_dk.lm[k] <- sum(z_di_s[1:l, 1:m_v, 1:m_n, s] == k)
                # 各添字が1のときに配列の構造が変わる対策(仮)
                for(cv in 1:old_v) {
                  if(cv < old_v){
                    tmp_count <- sum(z_di_s[1:d, cv, , s] == k)
                  } else if(cv == old_v) {
                    tmp_count <- sum(z_di_s[1:d, cv, 1:old_n, s] == k)
                  }
                }
                n_kv.lm[k, cv] <- tmp_count
              }
              
              # l,m要素について取り除く(iをvとnに分けているためi-1の処理がめんどいため)
              tmp_n_k <- rep(0, K) # 初期化
              tmp_n_k[z_di_s[l, m_v, m_n, s]] <- 1 # k番目に1を代入
              n_dk.lm <- n_dk.lm - tmp_n_k
              n_kv.lm[, m_v] <- n_kv.lm[, m_v] - tmp_n_k
              
              # リサンプリング確率を計算:式(3.179)
              term1 <- (n_kv.lm[, m_v] + beta_v[m_v]) / apply(t(n_kv.lm) + beta_v, 2, sum)
              term2 <- (n_dk.lm + alpha_k) / sum(n_dk.lm + alpha_k)
              tmp_q_z <- term1 * term2
              
              # 潜在トピックをサンプリング
              z_di_s[l, m_v, m_n, s] <- sample(1:K, size = 1, prob = tmp_q_z)
                
              # 重みを初期化
              w_z_di_s[l, m_v, m_n, ] <- 1 / S
              
            } ## (/リサンプリング)
          } ## (/活性化サンプル)
          
          # 動作確認
          print(paste0("d=", d, ", v=", v, "...Resampling"))
        }
        
        # 1期(語)前の添字を保存
        old_v <- v
        old_n <- n
        
      } ## (/各単語)
    }
  } ## (/各語彙)
  
  # 動作確認
  print(paste0("d=", d, "(", round(d / M * 100, 1), "%)...", round(Sys.time() - start_time, 2)))

} ## (/各文書)


# 推定結果の確認 -----------------------------------------------------------------

df_n_dv_wide <- cbind(
  as_tibble(n_dv), 
  doc = 1:M
)
df_n_dv_long <- pivot_longer(
  df_n_dv_wide, 
  cols = -doc, 
  names_to = "word", 
  names_prefix = "V", 
  names_ptypes = list(word = numeric()), 
  values_to = "freq"
)
df_n_dv <- tibble()
for(i in 1:nrow(df_n_dv_long)) {
  tmp_df <- tibble(
    doc = df_n_dv_long[["doc"]][i], 
    word = df_n_dv_long[["word"]][i], 
    freq = rep(1, df_n_dv_long[["freq"]][i]), 
    topic = as.factor(z_di_s[df_n_dv_long[["doc"]][i], df_n_dv_long[["word"]][i], 1:df_n_dv_long[["freq"]][i], S])
  )
  df_n_dv <- rbind(df_n_dv, tmp_df)
}

ggplot(df_n_dv, aes(x = word, y = doc, color = topic)) + 
  geom_point(position = "jitter") + 
  scale_x_continuous(breaks = seq(0, V, by = 5)) + 
  scale_y_continuous(breaks = seq(1, M))


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
w_z_di_k <- sum(w_z_di_s[d, v, n, ][z_di_s[d, v, n, ] == k])
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
  values_to = "prob"
)

ggplot(w_LongDF, aes(x = word, y = prob, fill = word, color = word)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  facet_wrap(~ topic, labeller = label_both) + 
  theme(legend.position = "none") + 
  scale_x_discrete(breaks = seq(1, V, by = 10)) + 
  labs(title = "Particle filter for LDA", 
       subtitle = "w")

w_LongDF %>% 
  filter(doc == 1:5) %>% 
  ggplot(aes(x = word, y = prob, fill = word, color = word)) + 
    geom_bar(stat = "identity", position = "dodge") + 
    facet_wrap(topic ~ doc, labeller = label_both, scales = "free_y") + 
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
  facet_wrap(~ word, scales = "free_y") + 
  labs(title = "Particle filter for LDA")

ggplot(n_kv_LongDF, aes(x = word, y = value, fill = word, color = word)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  facet_wrap(~ topic, labeller = label_both) + 
  theme(legend.position = "none", axis.text.x = element_text(angle = 90)) + 
  labs(title = "Particle filter for LDA")

# try ---------------------------------------------------------------------


# カウントを計算
for(k in 1:K) {
  n_dk[k] <- sum(z_di_s[d, 1:old_v, 1:old_n, s] == k)
  # 各添字が1のときに配列の構造が変わる対策(仮)
  for(cv in 1:old_v) {
    if(cv < old_v){
      tmp_count <- sum(z_di_s[1:d, cv, , s] == k)
    } else if(cv == old_v) {
      tmp_count <- sum(z_di_s[1:d, cv, 1:old_n, s] == k)
    }
  }
  n_kv[k, cv] <- tmp_count
}


z_di_k <- seq(0, 1, by = 0.01) %>% 
  sample(size = M * V * max(n_dv) * S, replace = TRUE) %>% 
  array(c(M, V, max(n_dv), K))
for(k in 1:K) {
  z_di_k[, , , k] <- z_di_k[, , , k] / apply(z_di_k, c(1, 2, 3), sum)
  #z_di_k[, , , k][n_dv == 0] <- 0
}
array(rep(2, 8), dim = c(2, 2, 2)) - array(rep(1, 8), dim = c(2, 2, 2))

