
# Text Handling -----------------------------------------------------------


# パッケージの読み込み
library(RMeCab)
library(tidyverse)


## 抽出する単語の指定
# 品詞(大分類)を指定
PoS_1 <- "名詞|^動詞|形容詞"

# 品詞(細分類)を指定
PoS_2 <- "一般|^自立"

# 最低出現頻度を指定
Freq <- 5

# 抽出しない単語を指定
stop_words <- "[a-z]"

# 形態素解析
mecab_df <- docDF("text_data", type = 1) # テキストファイルの保存先を指定する

# 文書dの語彙vの出現回数(N_dv)の集合
N_dv <- mecab_df %>% 
  filter(grepl(PoS_1, POS1)) %>%        # 指定した品詞(大分類)を取り出す
  filter(grepl(PoS_2, POS2)) %>%        # 指定した品詞(細分類)を取り出す
  filter(!grepl(stop_words, TERM)) %>%  # ストップワードを除く
  select(-c(TERM, POS1, POS2)) %>%      # 数値列のみを残す
  filter(apply(., 1, sum) >= Freq) %>%  # 指定した頻度以上の語彙を取り出す
  t()                                   # 転置

# 確認用の行列名
dimnames(N_dv) <- list(paste0("d=", 1:nrow(N_dv)), # 行名
                       paste0("v=", 1:ncol(N_dv))) # 列名

# 文書dの単語数(N_d)のベクトル
N_d <- apply(N_dv, 1, sum) # 行方向に和をとる

# 文書数
M <- nrow(N_dv)

# 総語彙数
V <- ncol(N_dv)


