
# Text Handling -----------------------------------------------------------


# 利用パッケージ
library(RMeCab)
library(tidyverse)


# 形態素解析
mecab_df <- docDF("text_data", type = 1) # テキストファイルの保存先を指定する


## 抽出する単語の指定
# 品詞(大分類)を指定
PoS_1 <- "名詞|^動詞|形容詞"

# 品詞(細分類)を指定
PoS_2 <- "一般|^自立"

# 最低出現頻度を指定
Freq <- 5

# 抽出しない単語を指定
stop_words <- "[a-z]"


# 文書ごとの各語彙の出現回数(単語数)
n_dv <- mecab_df %>% 
  filter(grepl(PoS_1, POS1)) %>%        # 指定した品詞(大分類)を取り出す
  filter(grepl(PoS_2, POS2)) %>%        # 指定した品詞(細分類)を取り出す
  filter(!grepl(stop_words, TERM)) %>%  # ストップワードを除く
  select(-c(TERM, POS1, POS2)) %>%      # 数値列のみを残す
  filter(apply(., 1, sum) >= Freq) %>%  # 指定した頻度以上の語彙を取り出す
  t()                                   # 転置


# 各文書の単語数
n_d <- apply(n_dv, 1, sum) # 行方向に和をとる

# 各語彙の単語数
n_v <- apply(n_dv, 2, sum) # 列方向に和をとる


# 文書数
M <- nrow(n_dv)

# 総語彙数
V <- ncol(n_dv)


