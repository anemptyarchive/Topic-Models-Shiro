
# Collapsed Variational Bayes Method for LDA -------------------------------------------

# 利用パッケージ
library(tidyverse)

## 簡易文書データ
# 文書数
M <- 10
# 語彙数
V <- 20
# 文書ごとの各語彙数
n_dv <- matrix(sample(1:5, M * V, replace = TRUE), M, V)


# パラメータの設定 ----------------------------------------------------------------

# イタレーション数
Iter <- 10

# トピック数
K <- 5

# 事前分布のパラメータ
alpha_k <- rep(2, K)
beta_v  <- rep(2, V)

# 潜在トピック集合の分布
z_dv_k <- new_z_dv_k <- array(0, dim = c(M, V, K))
for(d in 1:M) { ## (各文書)
  for(v in 1:V) { ## (各語彙)
    if(n_dv[d, v] > 0) {
      
      # ランダムに値を生成
      tmp_q_z <- sample(seq(0, 1, by = 0.01), size = K, replace = TRUE)
      
      # 値を正規化
      z_dv_k[d, v, ] <- tmp_q_z / sum(tmp_q_z)
    }
  }
}

## カウントの期待値
tmp_p_n <- array(0, dim = c(M, V, K))
for(k in 1:K) {
  tmp_p_n[, , k] <- z_dv_k[, , k] * n_dv
}

# 文書ごとにおいて各トピックが割り当てられた単語数
E_n_dk <- apply(tmp_p_n, c(1, 3), sum)

# 全文書において各トピックが割り当てられた単語数
E_n_kv <- apply(tmp_p_n, c(3, 2), sum)


# 周辺化変分ベイズ ----------------------------------------------------------------

# 推移の確認用
trace_alpha <- matrix(0, nrow = K, ncol = Iter + 1)
trace_beta <- matrix(0, nrow = V, ncol = Iter + 1)
# 初期値を代入
trace_alpha[, 1] <- alpha_k
trace_beta[, 1] <- beta_v

for(I in 1:Iter) { ## (イタレーション)
  
  # 動作確認用
  start_time <- Sys.time()
  
  for(d in 1:M) { ## (各文書)
    
    for(v in 1:V) { ## (各語彙)
      if(n_dv[d, v] > 0) {
        
        ## カウントの期待値を更新
        # n_{d,v}からd,i要素を除いたn_{d,v}^{\d,i}を作成
        n_dv.di <- n_dv
        n_dv.di[d, v] <- n_dv.di[d, v] - 1
        
        # 各単語数と掛け合わせる
        tmp_p_n <- array(0, dim = c(M, V, K))
        tmp_q_n <- array(0, dim = c(M, V, K))
        for(k in 1:K) {
          tmp_p_n[, , k] <- z_dv_k[, , k] * n_dv.di
          tmp_q_n[, , k] <- z_dv_k[, , k] * (1 - z_dv_k[, , k]) * n_dv.di
        }
        
        # d,i要素を除いたカウントの期待値と分散を計算
        E_n_dk.di <- apply(tmp_p_n, c(1, 3), sum)
        V_n_dk.di <- apply(tmp_q_n, c(1, 3), sum)
        E_n_kv.di <- apply(tmp_p_n, c(3, 2), sum)
        V_n_kv.di <- apply(tmp_q_n, c(3, 2), sum)
        
        # 潜在トピック集合の分布を計算:式(3.130)
        term1  <- (E_n_kv.di[, v] + beta_v[v]) / apply(t(E_n_kv.di) + beta_v, 2, sum) * (E_n_dk.di[d, ] + alpha_k)
        term21 <- V_n_kv.di[, v] / (2 * (E_n_kv.di[, v] + beta_v[v])^2)
        term22 <- V_n_dk.di[d, ] / (2 * (E_n_dk.di[d, ] + alpha_k)^2)
        term3  <- apply(V_n_kv.di, 1, sum) / (2 * apply(t(E_n_kv.di) + beta_v, 2, sum)^2)
        tmp_q_z <- term1 * exp(- term21 - term22) * exp(term3)
        
        # 値を正規化
        new_z_dv_k[d, v, ] <- tmp_q_z / sum(tmp_q_z)
        
      }
    } ## (/各語彙)
  } ## (/各文書)
  
  # 潜在トピック集合の分布を更新
  z_dv_k <- new_z_dv_k
  
  ## カウントの期待値を更新
  tmp_p_n <- array(0, dim = c(M, V, K))
  for(k in 1:K) {
    tmp_p_n[, , k] <- z_dv_k[, , k] * n_dv
  }
  
  # 文書ごとにおいて各トピックが割り当てられた単語数
  E_n_dk <- apply(tmp_p_n, c(1, 3), sum)
  
  # 全文書において各トピックが割り当てられた単語数
  E_n_kv <- apply(tmp_p_n, c(3, 2), sum)
  
  # 事後分布のパラメータを更新:式(3.191)
  alpha_numer <- apply(digamma(t(E_n_dk) + alpha_k) - digamma(alpha_k), 1, sum)
  alpha_denom <- sum(digamma(apply(t(E_n_dk) + alpha_k, 2, sum)) - digamma(sum(alpha_k)))
  alpha_k <- alpha_k * alpha_numer / alpha_denom
  
  # 単語分布のパラメータを更新:式(3.192)
  beta_numer <- apply(digamma(t(E_n_kv) + beta_v) - digamma(beta_v), 1, sum)
  beta_denom <- sum(digamma(apply(t(E_n_kv) + beta_v, 2, sum)) - digamma(sum(beta_v)))
  beta_v <- beta_v * beta_numer / beta_denom
  
  # 推移の確認用
  trace_alpha[, I + 1] <- alpha_k
  trace_beta[, I + 1] <- beta_v
  
  # 動作確認
  print(paste0(
    I, "th try...", 
    round(Sys.time() - start_time)
  ))
}

# 処理の検証用
sum(apply(E_n_dk, 1, sum) == apply(n_dv, 1, sum)) == M
sum(as.integer(apply(E_n_kv, 2, sum)) == as.integer(apply(n_dv, 2, sum))) == V
apply(E_n_dk.di, 1, sum) == apply(n_dv, 1, sum)
as.integer(apply(E_n_kv.di, 2, sum)) == as.integer(apply(n_dv, 2, sum))


# 推定結果の確認 -----------------------------------------------------------------

## トピック分布(期待値)
# thetaの期待値を計算:式(2.10)
theta_k <- alpha_k / sum(alpha_k)

# 作図用のデータフレームを作成
theta_df <- tibble(
  prob = theta_k, 
  topic = as.factor(1:K) # トピック番号
)

# 作図
ggplot(theta_df, aes(x = topic, y = prob, fill = topic)) + 
  geom_bar(stat = "identity", position = "dodge") + # 棒グラフ
  labs(title = "Collapsed Variational Bayes Method for LDA", 
       subtitle = expression(theta)) # ラベル


## 単語分布(期待値)
# phiの期待値を計算:式(2.10)
phi_v <- beta_v / sum(beta_v)

# 作図用のデータフレームを作成
phi_df <- tibble(
  prob = phi_v, 
  word = as.factor(1:V) # 語彙番号
)

# 作図
ggplot(phi_df, aes(x = word, y = prob, fill = word, color = word)) + 
  geom_bar(stat = "identity", position = "dodge") + # 棒グラフ
  scale_x_discrete(breaks = seq(1, V, by = 10)) + # x軸目盛
  theme(legend.position = "none") + # 凡例
  labs(title = "Collapsed Variational Bayes Method for LDA", 
       subtitle = expression(phi)) # ラベル


# 推移の確認用gif ---------------------------------------------------------------------

# 利用パッケージ
library(gganimate)


# データフレームに変換
trace_alpha_df <- data.frame()
trace_beta_df <- data.frame()
for(I in 1:(Iter + 1)) {
  
  # データフレームを作成
  tmp_trace_alpha <- data.frame(
    value = trace_alpha[, I], 
    topic = as.factor(1:K),  # トピック番号
    Iter = I - 1 # 試行回数
  )
  tmp_trace_beta <- data.frame(
    value = trace_beta[, I], 
    word = as.factor(1:V),  # 語彙番号
    Iter = I - 1 # 試行回数
  )
  
  # データフレームを結合
  trace_alpha_df <- rbind(trace_alpha_df, tmp_trace_alpha)
  trace_beta_df  <- rbind(trace_beta_df, tmp_trace_beta)
}


## トピック分布
# 作図
graph_alpha <- ggplot(trace_alpha_df, aes(topic, value, fill = topic)) + 
  geom_bar(stat = "identity", position = "dodge") + # 棒グラフ
  transition_manual(Iter) + # フレーム
  labs(title = "Collapsed Variational Bayes Method for LDA", 
       subtitle = "Iter={current_frame}") # ラベル

# 描画
animate(graph_alpha, nframes = (Iter + 1), fps = 10)


## 単語分布
# 作図
graph_beta <- ggplot(trace_beta_df, aes(word, value, fill = word, color = word)) + 
  geom_bar(stat = "identity", position = "dodge") + # 棒グラフ
  scale_x_discrete(breaks = seq(1, V, by = 10)) + # x軸目盛
  theme(legend.position = "none") + # 凡例
  transition_manual(Iter) + # フレーム
  labs(title = "Collapsed Variational Bayes Method for LDA", 
       subtitle = "Iter={current_frame}") # ラベル

# 描画
animate(graph_beta, nframes = (Iter + 1), fps = 10)


### 折れ線グラフ
## トピック分布のパラメータ
ggplot(trace_alpha_df, aes(x = Iter, y = value, color = topic)) + 
  geom_line() + # 折れ線グラフ
  labs(title = "Collapsed Variational Bayes Method for LDA", 
       subtitle = expression(alpha)) # ラベル


## 単語分布のパラメータ
ggplot(trace_beta_df, aes(x = Iter, y = value, color = word)) + 
  geom_line(alpha = 0.5) + # 
  theme(legend.position = "none") + 
  labs(title = "Collapsed Variational Bayes Method for LDA", 
       subtitle = expression(beta)) # ラベル



# try ---------------------------------------------------------------------

z_df <- data.frame()
for(k in 1:K) {
  tmp_z_df <- cbind(
    as.data.frame(z_dv_k[, , k]), 
    doc = as.factor(1:M), 
    topic = as.factor(k)
  )
  z_df <- rbind(z_df, tmp_z_df)
}
z_df_long <- pivot_longer(
  z_df, 
  cols = -c(doc, topic), 
  names_to = "word", 
  names_prefix = "V", 
  names_ptypes = list(word = factor()), 
  values_to = "prob"
)
z_df_long %>% 
  filter(doc %in% as.factor(1:5)) %>% 
  filter(word %in% as.factor(1:5)) %>% 
  ggplot(aes(x = topic, y = prob, fill = topic)) + 
    geom_bar(stat = "identity", position = "dodge") + 
    facet_wrap(doc ~ word, labeller = label_both)
