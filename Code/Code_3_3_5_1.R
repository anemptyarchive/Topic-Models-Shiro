
# Variational Bayes for LDA (1) -------------------------------------------

# 利用パッケージ
library(tidyverse)

## 動作検証用
# 文書数
M <- 10
# 語彙数
V <- 20
# 文書ごとの各語彙数
n_dv <- matrix(sample(1:3, M * V, replace = TRUE), M, V)


# パラメータの設定 ----------------------------------------------------------------

# イタレーション数
Iter <- 50

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
n_dk <- apply(tmp_z, c(1, 3), sum)

# 全文書において各トピックが割り当てられた単語数
n_kv <- apply(tmp_z, c(3, 2), sum)

# 処理の検証用
sum(n_dk) == sum(n_dv)
sum(n_kv) == sum(n_dv)

# 事後分布パラメータの初期値(=事前分布のパラメータ)
xi_dk <- t(t(n_dk) + alpha_k)
xi_kv <- t(t(n_kv) + beta_v)


# 変分ベイズ -------------------------------------------------------------------

# 推移の確認用のデータフレームを作成
trace_xi_theta <- cbind(
  as.data.frame(xi_dk), 
  doc = as.factor(1:M),  # 文書番号
  Iter = 0 # 試行回数
)
trace_xi_phi <- cbind(
  as.data.frame(xi_kv), 
  topic = as.factor(1:K),  # トピック番号
  Iter = 0 # 試行回数
)

for(I in 1:Iter) { ## (試行回数)
  
  for(d in 1:M) { ## (各文書)
    
    for(v in 1:V) { ## (各語彙)
      if(n_dv[d, v] > 0) {
        for(n in 1:n_dv[d, v]) { ## (各単語)
          
          # 潜在トピック集合の事後分布を計算:式(3.99)
          term1 <- digamma(xi_kv[, v]) - digamma(apply(xi_kv, 1, sum))
          term2 <- digamma(xi_dk[d, ]) - digamma(apply(xi_dk, 2, sum))
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
    n_dk <- apply(tmp_z, c(1, 3), sum)
    n_kv <- apply(tmp_z, c(3, 2), sum)
    
    # 事後分布のパラメータを計算:式(3.89)
    xi_dk[d, ] <- n_dk[d, ] + alpha_k
    
  } ## (/各文書)
  
  for(k in 1:K) { ## (各トピック)
    
    # 事後分布のパラメータ:式(3.95)
    xi_kv[k, ] <- n_kv[k, ] + beta_v
    
  } ## (/各トピック)
  
  # 推移の確認用のデータフレームを作成
  tmp_trace_xi_theta <- cbind(
    as.data.frame(xi_dk), 
    doc = as.factor(1:M),  # 文書番号
    Iter = I # 試行回数
  )
  tmp_trace_xi_phi <- cbind(
    as.data.frame(xi_kv), 
    topic = as.factor(1:K),  # トピック番号
    Iter = I # 試行回数
  )
  
  # データフレームを結合
  trace_xi_theta <- rbind(trace_xi_theta, tmp_trace_xi_theta)
  trace_xi_phi   <- rbind(trace_xi_phi, tmp_trace_xi_phi)
}


# 推定結果の確認 ----------------------------------------------------------------------

## トピック分布(期待値)
# thetaの期待値を計算:式(2.10)
theta_dk <- xi_dk / apply(xi_dk, 1, sum)

# 作図用のデータフレームを作成
theta_WideDF <- cbind(
  as.data.frame(theta_dk), 
  doc = as.factor(1:M) # 文書番号
)

# データフレームをlong型に変換
theta_LongDF <- pivot_longer(
  theta_WideDF, 
  cols = -doc,         # 変換せずにそのまま残す現列名
  names_to = "topic",  # 現列名を格納する新しい列の名前
  names_prefix = "V",  # 現列名から取り除く文字列
  names_ptypes = list(topic = factor()),  # 現列名を要素とする際の型
  values_to = "prob"   # 現要素を格納する新しい列の名前
)

# 作図
ggplot(theta_LongDF, aes(x = topic, y = prob, fill = topic)) + 
  geom_bar(stat = "identity", position = "dodge") + # 棒グラフ
  facet_wrap( ~ doc, labeller = label_both) + # グラフの分割
  labs(title = "Variational Bayes for LDA (1)", 
       subtitle = expression(Theta)) # ラベル


## 単語分布(期待値)
# phiの期待値を計算:式(2.10)
phi_kv <- xi_kv / apply(xi_kv, 1, sum)

# 作図用のデータフレームを作成
phi_WideDF <- cbind(
  as.data.frame(phi_kv), 
  topic = as.factor(1:K) # トピック番号
)

# データフレームをlong型に変換
phi_LongDF <- pivot_longer(
  phi_WideDF, 
  cols = -topic,       # 変換せずにそのまま残す現列名
  names_to = "word",   # 現列名を格納する新しい列の名前
  names_prefix = "V",  # 現列名から取り除く文字列
  names_ptypes = list(word = factor()),  # 現列名を要素とする際の型
  values_to = "prob"   # 現要素を格納する新しい列の名前
)

# 作図
ggplot(phi_LongDF, aes(x = word, y = prob, fill = word, color = word)) + 
  geom_bar(stat = "identity", position = "dodge") + # 棒グラフ
  facet_wrap( ~ topic, labeller = label_both) + # グラフの分割
  scale_x_discrete(breaks = seq(1, V, by = 10)) + # x軸目盛
  theme(legend.position = "none") + # 凡例
  labs(title = "Variational Bayes for LDA (1)", 
       subtitle = expression(Phi)) # ラベル


# 推移の確認用gif ---------------------------------------------------------------------

# 利用パッケージ
library(gganimate)


## トピック分布
# データフレームをlong型に変換
trace_xi_theta_LongDF <- pivot_longer(
  trace_xi_theta, 
  cols = -c(doc, Iter),          # 変換せずにそのまま残す現列名
  names_to = "topic",   # 現列名を格納する新しい列の名前
  names_prefix = "V",  # 現列名から取り除く文字列
  names_ptypes = list(topic = factor()),  # 現列名を要素とする際の型
  values_to = "value"    # 現要素を格納する新しい列の名前
)

# 作図
graph_xi_theta <- ggplot(trace_xi_theta_LongDF, aes(x = topic, y = value, fill = topic)) + 
  geom_bar(stat = "identity", position = "dodge") +  # 棒グラフ
  facet_wrap( ~ doc, labeller = label_both) +        # グラフの分割
  transition_manual(Iter) + 
  labs(title = "Variational Bayes for LDA (1):xi^theta_dk", 
       subtitle = "Iter={current_frame}") # ラベル

# 描画
animate(graph_xi_theta, nframes = (Iter + 1), fps = 10)


## 単語分布
# 作図用のデータフレームを作成
trace_xi_phi_LongDF <- pivot_longer(
  trace_xi_phi, 
  cols = -c(topic, Iter),        # 変換せずにそのまま残す現列名
  names_to = "word",    # 現列名を格納する新しい列の名前
  names_prefix = "V",  # 現列名から取り除く文字列
  names_ptypes = list(word = factor()),  # 現列名を要素とする際の型
  values_to = "value"    # 現要素を格納する新しい列の名前
)

# 作図
graph_xi_phi <- ggplot(trace_xi_phi_LongDF, aes(x = word, y = value, fill = word, color = word)) + 
  geom_bar(stat = "identity", position = "dodge") +  # 棒グラフ
  facet_wrap( ~ topic, labeller = label_both) +      # グラフの分割
  scale_x_discrete(breaks = seq(1, V, by = 10)) +    # x軸目盛
  theme(legend.position = "none") +                  # 凡例
  transition_manual(Iter) + 
  labs(title = "Variational Bayes for LDA (1):xi^phi_kv", 
       subtitle = "Iter={current_frame}") # ラベル

# 描画
animate(graph_xi_phi, nframes = (Iter + 1), fps = 10)


