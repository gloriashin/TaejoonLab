library(tidyverse)
library(ggplot2)

# 1. TSV 파일 목록 가져오기
file_list <- list.files(path = "/media/seungyun/8TBHardDisk1/XL_Plou2017/hisat2-tsv/st12", pattern = "\\.tsv$", full.names = TRUE)

# 2. 각 파일에서 TPM 데이터 추출 및 컬럼명 생성
tpm_data_list <- lapply(file_list, function(file) {
  # 파일 읽기
  df <- read_tsv(file, show_col_types = FALSE)
  
  # 파일명에서 컬럼명 추출 (st 뒤의 숫자와 그 다음 _ 사이 문자)
  filename <- basename(file)
  # st와 _ 사이의 패턴 추출 (예: st12_NPp에서 12_NPp)
  col_name <- str_extract(filename, "st\\d+_[^_]+")
  col_name <- str_remove(col_name, "st")  # st 제거하여 12_NPp 형태로
  
  # Gene ID와 TPM 컬럼만 선택하고 컬럼명 변경
  df_selected <- df %>%
    select('Gene ID', TPM) %>%
    rename(!!col_name := TPM)
  
  return(df_selected)
})

# 3. Gene ID 기준으로 모든 데이터 합치기
tpm_combined <- tpm_data_list %>%
  reduce(full_join, by = "Gene ID")

# 4. Gene ID에서 "GeneID:" 문자열 제거
tpm_combined <- tpm_combined %>%
  mutate(`Gene ID` = str_replace(`Gene ID`, "GeneID:", ""))

# 5. 최소 2개 replicate에서 발현하는 유전자만 필터링
# (각 stage의 replicate가 3개씩 있다고 가정)
# TPM > 0을 발현으로 간주
tpm_filtered <- tpm_combined %>%
  rowwise() %>%
  mutate(n_expressed = sum(c_across(-`Gene ID`) > 0, na.rm = TRUE)) %>%
  ungroup() %>%
  filter(n_expressed >= 2) %>%
  select(-n_expressed)

# 6. PCA 수행을 위한 데이터 준비
# Gene ID를 행 이름으로 설정하고 숫자 데이터만 추출
tpm_matrix <- tpm_filtered %>%
  column_to_rownames("Gene ID") %>%
  as.matrix()

# NA를 0으로 대체
tpm_matrix[is.na(tpm_matrix)] <- 0

# 7. PCA 수행 (샘플이 컬럼에 있으므로 transpose 필요)
tpm_matrix_t <- t(tpm_matrix)
pca_result <- prcomp(tpm_matrix_t, scale. = TRUE, center = TRUE)

# 8. PCA 결과 시각화
# PCA 스코어 데이터프레임 생성
pca_scores <- as.data.frame(pca_result$x)
pca_scores$Sample <- rownames(pca_scores)

# Stage 정보 추출 (숫자 부분)
pca_scores$Stage <- str_extract(pca_scores$Sample, "^\\d+")

# 분산 설명 비율 계산
var_explained <- summary(pca_result)$importance[2,] * 100

# PCA plot (PC1 vs PC2)
p1 <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = Stage, label = Sample)) +
  geom_point(size = 3) +
  geom_text(vjust = -0.5, hjust = 0.5, size = 3) +
  labs(
    title = "PCA of TPM Expression Data",
    x = paste0("PC1 (", round(var_explained[1], 2), "%)"),
    y = paste0("PC2 (", round(var_explained[2], 2), "%)")
  ) +
  theme_bw() +
  theme(legend.position = "right")

print(p1)

# 9. Scree plot (분산 설명 비율)
scree_data <- data.frame(
  PC = paste0("PC", 1:length(var_explained)),
  Variance = var_explained[1:min(10, length(var_explained))]
)
scree_data$PC <- factor(scree_data$PC, levels = scree_data$PC)

p2 <- ggplot(scree_data, aes(x = PC, y = Variance)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_line(aes(group = 1), color = "red") +
  geom_point(color = "red", size = 2) +
  labs(
    title = "Scree Plot",
    x = "Principal Component",
    y = "Variance Explained (%)"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p2)

# 10. 결과 저장
# 필터링된 TPM 데이터 저장
write_csv(tpm_filtered, "tpm_filtered_combined.csv")

# PCA 스코어 저장
write_csv(pca_scores, "pca_scores.csv")

# 요약 정보 출력.
cat("\n=== 분석 요약 ===\n")
cat("총 파일 수:", length(file_list), "\n")
cat("총 유전자 수 (필터링 전):", nrow(tpm_combined), "\n")
cat("필터링된 유전자 수:", nrow(tpm_filtered), "\n")
cat("샘플 수:", ncol(tpm_filtered) - 1, "\n")
cat("\n상위 5개 PC의 분산 설명 비율:\n")
print(var_explained[1:5])

