
# Particle filter for LDA -------------------------------------------------

# 利用パッケージ
library(tidyverse)

## 簡易文書データ
# 文書数
M <- 10
# 語彙数
V <- 20
# 文書ごとの各語彙数
n_dv <- matrix(sample(0:10, M * V, replace = TRUE), M, V)


# パラメータの設定 ----------------------------------------------------------------

# サンプリング数(粒子の数)を指定
S <- 50

# リサンプリング数を指定
R <- 10

# 閾値を指定
threshold <- 10

# トピック数を指定
K <- 5

# 事前分布のパラメータを指定
alpha_k <- rep(2, K)
beta_v  <- rep(2, V)

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
w_di_s <- array(1/S, dim = c(M, V, max(n_dv), S))

# 割り当てられたトピックに関する単語数の受け皿
n_kv <- matrix(0, nrow = K, ncol = V)


# パーティクルフィルタ --------------------------------------------------------------

for(d in 1:M) { ## (各文書)
  
  # 動作確認用
  start_time <- Sys.time()
  
  # 文書dにおいて各トピックが割り当て単語数の初期化
  n_dk <- rep(0, K) # (ここではベクトルで扱うことに注意)
  
  # 1期(語)前の値i-1を初期化
  old_v <- old_n <- 1 # (本当は0でなきゃ)

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
          
          # 重みを計算:式(3.175),(3.176)
          term1 <- (n_kv[, v] + beta_v[v]) / apply(t(n_kv) + beta_v, 2, sum)
          term2 <- (n_dk + alpha_k) / sum(n_dk + alpha_k)
          w_di_s[d, v, n, s] <- w_di_s[d, old_v, old_n, s] * sum(term1 * term2)
          
        } ## (/サンプリング)
        
        # 重みを正規化
        w_di_s[d, v, n, ] <- w_di_s[d, v, n, ] / sum(w_di_s[d, v, n, ])
        
        # ESSを計算
        ESS <- 1 / sum(w_di_s[d, v, n, ] ^ 2)
        
        if(ESS < threshold) {
          
          for(r in 1:R) { ## (活性化サンプル)
            
            # 活性化する重みをサンプリング(仮):(n_dv[1:(d-1), 1:V(大文字)]とn_dv[d, 1:v(小文字)]だけ取り出したい)
            vec_n_dv <- as.vector(t(n_dv))
            re_w <- sample(
              1:((d - 1) * V + v), 
              size = 1, 
              prob = vec_n_dv[1:((d - 1) * V + v)] # 各語彙の出現回数を確率として使用(単語レベルで一様にするため)
            )
            l <- re_w %/% V + 1
            mv <- re_w %% V
            mn <- sample(1:n_dv[l, mv], size = 1)
            
            for(s in 1:S) { ## (リサンプリング)
              
              # 割り当てられたトピックに関する単語数を初期化
              n_dk.lm <- rep(0, K)
              n_kv.lm <- matrix(0, nrow = K, ncol = V)
              
              # 潜在トピックに関するカウントを計算
              for(k in 1:K) {
                n_dk.lm[k] <- sum(z_di_s[1:d, 1:v, 1:n, s] == k)
                # 各添字が1のときに配列の構造が変わる対策(仮)
                for(cv in 1:v) {
                  if(cv < v){
                    tmp_count <- sum(z_di_s[1:d, cv, , s] == k)
                  } else if(cv == v) {
                    tmp_count <- sum(z_di_s[1:d, cv, 1:n, s] == k)
                  }
                }
                n_kv.lm[k, cv] <- tmp_count
              }
              
              # l,m要素について取り除く(iをvとnに分けているためi-1の処理がめんどいため)
              tmp_n_k <- rep(0, K) # 初期化
              tmp_n_k[z_di_s[l, mv, mn, s]] <- 1 # k番目に1を代入
              n_dk.lm <- n_dk.lm - tmp_n_k
              n_kv.lm[, mv] <- n_kv.lm[, mv] - tmp_n_k
              
              # リサンプリング確率を計算:式(3.179)
              term1 <- (n_kv.lm[, v] + beta_v[v]) / apply(t(n_kv.lm) + beta_v, 2, sum)
              term2 <- (n_dk.lm + alpha_k) / sum(n_dk.lm + alpha_k)
              tmp_q_z <- term1 * term2
              
              # 潜在トピックをサンプリング
              z_di_s[l, mv, mn, s] <- sample(1:K, size = 1, prob = tmp_q_z)
                
              # 重みを初期化
              w_di_s <- array(1 / S, dim = c(M, V, max(n_dv), S))
              
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


# try:パラメータの推定 ------------------------------------------------------------

# トピック集合の分布を近似:式(3.177)
tmp_p_z <- array(0, dim = c(M, V, max(n_dv), K))
for(k in 1:K) {
  # トピックがkの要素についてs方向に和をとる
  tmp_p_z[, , , k] <- apply(w_di_s * (z_di_s == k), c(1, 2, 3), sum)
}

# k方向の和が1となるように正規化
p_z_di_k <- array(0, dim = c(M, V, max(n_dv), K))
for(k in 1:K) {
  p_z_di_k[, , , k] <- tmp_p_z[, , , k] / apply(tmp_p_z, c(1, 2, 3), sum)
}

# n_dvが0(出現回数が0)について0除算によるNaNとなるのでそれを0に置換
p_z_di_k[is.nan(p_z_di_k)] <- 0

# 処理の検証
near(apply(p_z_di_k, c(1, 2), sum), n_dv)


# 統計量を計算
E_n_dk <- apply(p_z_di_k, c(1, 4), sum)
E_n_kv <- apply(p_z_di_k, c(4, 2), sum)

# ハイパーパラメータを計算:式(3.89),(3.95)
xi_theta_dk <- t(t(E_n_dk) + alpha_k)
xi_phi_kv   <- t(t(E_n_kv) + beta_v)


# try:推定結果の確認 -----------------------------------------------------------------

## トピック分布(期待値)
# thetaの期待値を計算:式(2.10)
theta_dk <- xi_theta_dk / apply(xi_theta_dk, 1, sum)

# 作図用のデータフレームを作成
theta_df_wide <- cbind(
  as.data.frame(theta_dk), 
  doc = as.factor(1:M)
)

# データフレームをlong型に変換
theta_df_long <- pivot_longer(
  theta_df_wide, 
  cols = -doc, # 変換せずにそのまま残す現列名
  names_to = "topic", # 現列名を格納する新しい列の名前
  names_prefix = "V", # 現列名から取り除く文字列
  names_ptypes = list(topic = factor()),  # 現列名を要素とする際の型
  values_to = "prob" # 現要素を格納する新しい列の名前
)

# 作図
ggplot(theta_df_long, aes(x = topic, y = prob, fill = topic)) + 
  geom_bar(stat = "identity", position = "dodge") + # 棒グラフ
  facet_wrap( ~ doc, labeller = label_both) + # グラフの分割
  labs(title = "Particle Filter for LDA", 
       subtitle = expression(Theta)) # ラベル


## 単語分布(期待値)
# phiの期待値を計算:式(2.10)
phi_kv <- xi_phi_kv / apply(xi_phi_kv, 1, sum)

# 作図用のデータフレームを作成
phi_df_wide <- cbind(
  as.data.frame(phi_kv), 
  topic = as.factor(1:K) # トピック番号
)

# データフレームをlong型に変換
phi_df_long <- pivot_longer(
  phi_df_wide, 
  cols = -topic, # 変換せずにそのまま残す現列名
  names_to = "word", # 現列名を格納する新しい列の名前
  names_prefix = "V", # 現列名から取り除く文字列
  names_ptypes = list(word = factor()),  # 現列名を要素とする際の型
  values_to = "prob" # 現要素を格納する新しい列の名前
)

# 作図
ggplot(phi_df_long, aes(x = word, y = prob, fill = word, color = word)) + 
  geom_bar(stat = "identity", position = "dodge") + # 棒グラフ
  facet_wrap( ~ topic, labeller = label_both) + # グラフの分割
  scale_x_discrete(breaks = seq(0, V, by = 10)) + # x軸目盛
  theme(legend.position = "none") + # 凡例
  labs(title = "Particle Filter for LDA", 
       subtitle = expression(Phi)) # ラベル


# try ---------------------------------------------------------------------

n_df_wide <- cbind(
  as_tibble(n_dv), 
  doc = 1:M
)
n_df_long <- pivot_longer(
  n_df_wide, 
  cols = -doc, 
  names_to = "word", 
  names_prefix = "V", 
  names_ptypes = list(word = numeric()), 
  values_to = "freq"
)

particle_df <- tibble()
for(i in 1:nrow(n_df_long)) {
  if(n_df_long[["freq"]][i] > 0) {
    for(s in 1:S) {
      tmp_df <- tibble(
        doc = n_df_long[["doc"]][i], 
        word = n_df_long[["word"]][i], 
        topic = as.factor(z_di_s[doc, word, 1:n_df_long[["freq"]][i], s]), 
        particle = s
      )
      particle_df <- rbind(particle_df, tmp_df)
    }
  }
}
particle_df %>% 
  filter(particle == S) %>% 
  ggplot(aes(x = word, y = doc, color = topic)) + 
    geom_point(position = "jitter") + 
    scale_x_continuous(breaks = seq(0, V, by = 5)) + 
    scale_y_continuous(breaks = seq(1, M))


library(gganimate)

particle_graph <- ggplot(particle_df, aes(x = word, y = doc, color = topic)) + 
  geom_point(position = "jitter") + 
  scale_x_continuous(breaks = seq(0, V, by = 5)) + 
  scale_y_continuous(breaks = seq(1, M)) + 
  transition_manual(particle) + 
  labs(title = "Particle Filter for LDA", 
       subtitle = "particle={current_frame}")

animate(particle_graph, nframes = S, fps = 10)