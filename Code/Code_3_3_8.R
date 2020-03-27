
# Collapsed Variational Bayes Method for LDA -------------------------------------------

# 利用パッケージ
library(tidyverse)


# パラメータの設定 ----------------------------------------------------------------

# イタレーション数
Iter <- 10

# トピック数
K <- 5

# 事前分布のパラメータ
alpha_k <- rep(2, K)
beta_v  <- rep(2, V)

# 潜在トピック集合の分布
z_dvk <- array(1 / K, dim = c(M, V, K))

## カウントの期待値
tmp_z <- array(0, dim = c(M, V, K))
for(k in 1:K) {
  tmp_z[, , k] <- z_dvk[, , k] * n_dv
}

# 文書ごとにおいて各トピックが割り当てられた単語数
n_dk <- apply(tmp_z, c(1, 3), sum)

# 全文書において各トピックが割り当てられた単語数
n_kv <- apply(tmp_z, c(3, 2), sum)

# 処理の検証用
sum(n_dk) == sum(n_dv)
sum(n_kv) == sum(n_dv)


# 周辺化変分ベイズ ----------------------------------------------------------------

# 受け皿
new_z_dvk <- array(0, c(M, V, K))
tmp_z_dvk.di <- array(0, c(M, V, K))

# 推移の確認用データフレームを作成
trace_alpha <- tibble(
  value = alpha_k, 
  topic = as.factor(1:K),  # トピック番号
  Iter = 0 # 試行回数
)
trace_beta <- tibble(
  value = beta_v, 
  word = as.factor(1:V),  # 語彙番号
  Iter = 0 # 試行回数
)

for(I in 1:Iter) { ## (イタレーション)
  
  for(d in 1:M) { ## (各文書)
    
    for(v in 1:V) { ## (各語彙)
      if(n_dv[d, v] > 0) {
        for(n in 1:n_dv[d, v]) { ## (各単語)
          
          # カウントを更新
          z_dvk.di <- z_dvk
          z_dvk.di[d, v, ] <- 0
          
          E_n_dk.di <- apply(z_dvk.di, c(1, 3), sum)
          V_n_dk.di <- apply(z_dvk.di * (1 - z_dvk.di), c(1, 3), sum)
          
          for(k in 1:K) {
            tmp_z_dvk.di[, , k] <- z_dvk.di[, , k] * n_dv
          }
          E_n_kv.di <- apply(tmp_z_dvk.di, c(3, 2), sum)
          V_n_kv.di <- apply(tmp_z_dvk.di * (1 - tmp_z_dvk.di), c(3, 2), sum)
          
          
          # 潜在トピック集合の分布を計算:式(3.130)
          term1 <- (E_n_kv.di[, v] + beta_v[v]) / apply(t(E_n_kv.di) + beta_v, 2, sum) * (E_n_dk.di[d, ] + alpha_k)
          term21 <- V_n_kv.di[, v] / (2 * (E_n_kv.di[, v] + beta_v[v])^2)
          term22 <- V_n_dk.di[d, ] / (2 * (E_n_dk.di[d, ] + alpha_k)^2)
          term3 <- apply(V_n_kv.di, 1, sum) / (2 * apply(t(E_n_kv.di) + beta_v, 2, sum)^2)
          new_z_dvk[d, v, ] <- term1 * exp(- term21 - term22) * exp(term3)
          
        } ## (/各単語)
      }
    } ## (/各語彙)
  } ## (/各文書)
  
  # 潜在トピック集合の分布を更新
  z_dvk <- new_z_dvk
  
  # カウントを更新
  n_dk <- apply(z_dvk, c(1, 3), sum)
  n_kv <- apply(z_dvk, c(3, 2), sum)
  
  # 事後分布のパラメータを更新:式(3.191)
  alpha_numer <- apply(digamma(t(n_dk) + alpha_k) - digamma(alpha_k), 1, sum)
  alpha_denom <- sum(digamma(apply(t(n_dk) + alpha_k, 2, sum)) - digamma(sum(alpha_k)))
  alpha_k <- alpha_k * alpha_numer / alpha_denom
  
  # 単語分布のパラメータを更新:式(3.192)
  beta_numer <- apply(digamma(t(n_kv) + beta_v) - digamma(beta_v), 1, sum)
  beta_denom <- sum(digamma(apply(t(n_kv) + beta_v, 2, sum)) - digamma(sum(beta_v)))
  beta_v <- beta_v * beta_numer / beta_denom
  
  # 推移の確認用データフレームを作成
  tmp_trace_alpha <- tibble(
    value = alpha_k, 
    topic = as.factor(1:K),  # トピック番号
    Iter = I # 試行回数
  )
  tmp_trace_beta <- tibble(
    value = beta_v, 
    word = as.factor(1:V),  # 語彙番号
    Iter = I # 試行回数
  )
  
  # データフレームを結合
  trace_alpha <- rbind(trace_alpha, tmp_trace_alpha)
  trace_beta  <- rbind(trace_beta, tmp_trace_beta)
}



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
ggplot(phi_df, aes(x = word, y = prob, fill = word)) + 
  geom_bar(stat = "identity", position = "dodge") + # 棒グラフ
  scale_x_discrete(breaks = seq(1, V, by = 10)) + # x軸目盛
  theme(legend.position = "none") + # 凡例
  labs(title = "Collapsed Variational Bayes Method for LDA", 
       subtitle = expression(phi)) # ラベル


# 推移の確認用gif ---------------------------------------------------------------------

# 利用パッケージ
library(gganimate)

## トピック分布
# 作図
graph_alpha <- ggplot(trace_alpha, aes(topic, value, fill = topic)) + 
  geom_bar(stat = "identity", position = "dodge") +  # 棒グラフ
  transition_manual(Iter) + 
  labs(title = "Collapsed Variational Bayes Method for LDA", 
       subtitle = "Iter={current_frame}") # ラベル

# 描画
animate(graph_alpha)


## 単語分布
# 作図
graph_beta <- ggplot(trace_beta, aes(word, value, fill = word)) + 
  geom_bar(stat = "identity", position = "dodge") +  # 棒グラフ
  scale_x_discrete(breaks = seq(1, V, by = 10)) +    # x軸目盛
  theme(legend.position = "none") +                  # 凡例
  transition_manual(Iter) + 
  labs(title = "Collapsed Variational Bayes Method for LDA", 
       subtitle = "Iter={current_frame}") # ラベル

# 描画
animate(graph_beta)




# 推定結果の確認 -----------------------------------------------------------------


# 
z_WideDF <- NULL
for(k in 1:k) {
  tmp_z_df <- cbind(
    as.data.frame(z_dvk), 
    doc = as.factor(1:M),  # 文書番号
    topic = as.factor(k) # トピック番号
  )
  z_WideDF <- rbind(z_WideDF, tmp_z_df)
}

# データフレームをlong型に変換
z_LongDF <- pivot_longer(
  z_WideDF, 
  cols = -c(doc, topic),  # 変換せずにそのまま残す現列名
  names_to = "word",      # 現列名を格納する新しい列の名前
  names_prefix = "V",     # 現列名から取り除く文字列
  names_ptypes = list(word = factor()),  # 現列名を要素とする際の型
  values_to = "prob"      # 現要素を格納する新しい列の名前
)

# 作図
z_LongDF %>% 
  filter(doc %in% 1:5) %>% 
  filter(word %in% 1:5) %>% 
  ggplot(aes(x = topic, y = prob, fill = topic)) + 
    geom_bar(stat = "identity", position = "dodge") + # 棒グラフ
    facet_wrap(doc ~ word, labeller = label_both) + # グラフの分割
#  scale_x_discrete(breaks = seq(1, V, by = 10)) + # x軸目盛
    theme(legend.position = "none") + # 凡例
    labs(title = "Collapsed Variational Bayes Method for LDA", 
         subtitle = expression(z[dn])) # ラベル



