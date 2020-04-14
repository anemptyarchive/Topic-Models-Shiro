
# Variational Bayes for LDA (1) -------------------------------------------

# 利用パッケージ
library(tidyverse)

## 簡易文書データ
# 文書数
M <- 10
# 語彙数
V <- 20
# 文書ごとの各語彙数
n_dv <- matrix(sample(1:10, M * V, replace = TRUE), M, V)


# パラメータの設定 ----------------------------------------------------------------

# イタレーション数を指定
Iter <- 50

# トピック数を指定
K <- 5

# 事前分布のパラメータを指定
alpha_k <- rep(2, K)
beta_v  <- rep(2, V)

# 潜在トピック集合の分布
z_dv_k <- array(0, dim = c(M, V, K))
for(d in 1:M) {
  for(v in 1:V) {
    tmp_q_z <- seq(0, 1, by = 0.01) %>% 
      sample(size = K, replace = TRUE)
    z_dv_k[d, v, ] <- tmp_q_z / sum(tmp_q_z)
  }
}

## カウントの期待値
tmp_n <- array(0, dim = c(M, V, K))
for(k in 1:K) {
  tmp_n[, , k] <- z_dv_k[, , k] * n_dv
}

# 文書ごとにおいて各トピックが割り当てられた単語数の期待値の初期値
E_n_dk <- apply(tmp_n, c(1, 3), sum)

# 全文書において各トピックが割り当てられた単語数の期待値の初期値
E_n_kv <- apply(tmp_n, c(3, 2), sum)

# 事後分布パラメータの初期値
xi_dk <- t(t(E_n_dk) + alpha_k)
xi_kv <- t(t(E_n_kv) + beta_v)


# 変分ベイズ -------------------------------------------------------------------

# 推移の確認用
trace_xi_theta <- array(0, dim = c(M, K, Iter + 1))
trace_xi_phi   <- array(0, dim = c(K, V, Iter + 1))
# 初期値を代入
trace_xi_theta[, , 1] <- xi_dk
trace_xi_phi[, , 1]   <- xi_kv

for(I in 1:Iter) { ## (試行回数)
  
  for(d in 1:M) { ## (各文書)
    
    for(v in 1:V) { ## (各語彙)
      if(n_dv[d, v] > 0) {
        
        # 潜在トピック集合の事後分布を計算:式(3.99)
        term1 <- digamma(xi_kv[, v]) - digamma(apply(xi_kv, 1, sum))
        term2 <- digamma(xi_dk[d, ]) - digamma(apply(xi_dk, 2, sum))
        tmp_q_z <- exp(term1) * exp(term2)
        # 正規化
        z_dv_k[d, v, ] <- tmp_q_z / sum(tmp_q_z)
        
      }
    } ## (/各語彙)
    
    # カウントの期待値
    for(k in 1:K) {
      tmp_n[, , k] <- z_dv_k[, , k] * n_dv
    }
    E_n_dk <- apply(tmp_n, c(1, 3), sum)
    E_n_kv <- apply(tmp_n, c(3, 2), sum)
    
    # 事後分布のパラメータを計算:式(3.89)
    xi_dk[d, ] <- E_n_dk[d, ] + alpha_k
    
  } ## (/各文書)
  
  for(k in 1:K) { ## (各トピック)
    
    # 事後分布のパラメータ:式(3.95)
    xi_kv[k, ] <- E_n_kv[k, ] + beta_v
    
  } ## (/各トピック)
  
  # 推移の確認用
  trace_xi_theta[, , I + 1] <- xi_dk
  trace_xi_phi[, , I + 1]   <- xi_kv
  
}

# 処理の検証用
sum(E_n_dk) == sum(n_dv)
sum(E_n_kv) == sum(n_dv)


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


# 推移の確認 ---------------------------------------------------------------------

# 利用パッケージ
library(gganimate)


# データフレームに変換
trace_xi_theta_WideDF <- data.frame()
trace_xi_phi_WideDF <- data.frame()
for(I in 1:(Iter + 1)) {
  # データフレームに変換
  tmp_trace_xi_theta <- cbind(
    as.data.frame(trace_xi_theta[, , I]), 
    doc = as.factor(1:M), # 文書番号
    Iter = I - 1 # 試行回数
  )
  tmp_trace_xi_phi <- cbind(
    as.data.frame(trace_xi_phi[, , I]), 
    topic = as.factor(1:K), # トピック番号
    Iter = I - 1 # 試行回数
  )
  # 結合
  trace_xi_theta_WideDF <- rbind(trace_xi_theta_WideDF, tmp_trace_xi_theta)
  trace_xi_phi_WideDF <- rbind(trace_xi_phi_WideDF, tmp_trace_xi_phi)
}

# データフレームをlong型に変換
trace_xi_theta_LongDF <- pivot_longer(
  trace_xi_theta_WideDF, 
  cols = -c(doc, Iter),          # 変換せずにそのまま残す現列名
  names_to = "topic",   # 現列名を格納する新しい列の名前
  names_prefix = "V",  # 現列名から取り除く文字列
  names_ptypes = list(topic = factor()),  # 現列名を要素とする際の型
  values_to = "value"    # 現要素を格納する新しい列の名前
)
trace_xi_phi_LongDF <- pivot_longer(
  trace_xi_phi_WideDF, 
  cols = -c(topic, Iter),        # 変換せずにそのまま残す現列名
  names_to = "word",    # 現列名を格納する新しい列の名前
  names_prefix = "V",  # 現列名から取り除く文字列
  names_ptypes = list(word = factor()),  # 現列名を要素とする際の型
  values_to = "value"    # 現要素を格納する新しい列の名前
)


## トピック分布
# 作図
graph_xi_theta <- ggplot(trace_xi_theta_LongDF, aes(x = topic, y = value, fill = topic)) + 
  geom_bar(stat = "identity", position = "dodge") +  # 棒グラフ
  facet_wrap( ~ doc, labeller = label_both) +        # グラフの分割
  transition_manual(Iter) + 
  labs(title = "Variational Bayes for LDA (1)", 
       subtitle = "Iter={current_frame}") # ラベル

# 描画
animate(graph_xi_theta, nframes = (Iter + 1), fps = 10)


## 単語分布
# 作図
graph_xi_phi <- ggplot(trace_xi_phi_LongDF, aes(x = word, y = value, fill = word, color = word)) + 
  geom_bar(stat = "identity", position = "dodge") +  # 棒グラフ
  facet_wrap( ~ topic, labeller = label_both) +      # グラフの分割
  scale_x_discrete(breaks = seq(1, V, by = 10)) +    # x軸目盛
  theme(legend.position = "none") +                  # 凡例
  transition_manual(Iter) + 
  labs(title = "Variational Bayes for LDA (1)", 
       subtitle = "Iter={current_frame}") # ラベル

# 描画
animate(graph_xi_phi, nframes = (Iter + 1), fps = 10)


### 折れ線グラフ
## トピック分布のパラメータ
# 文書番号を指定
doc_num <- 10

# 作図
trace_xi_theta_LongDF %>% 
  filter(doc == doc_num) %>% 
  ggplot(aes(x = Iter, y = value, color = topic)) + 
    geom_line() + 
    labs(title = "Variational Bayes for LDA (1)", 
         subtitle = expression(xi^theta)) # ラベル


## トピック分布のパラメータ
# トピック番号を指定
topic_num <- 1

# 作図
trace_xi_phi_LongDF %>% 
  filter(topic == topic_num) %>% 
  ggplot(aes(x = Iter, y = value, color = word)) + 
    geom_line(alpha = 0.5) + 
    theme(legend.position = "none") + # 凡例
    labs(title = "Variational Bayes for LDA (1)", 
         subtitle = expression(xi^phi)) # ラベル


