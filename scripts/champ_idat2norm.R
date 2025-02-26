# R腳本 (champ_idat_processor.R)


# ✅ **檢查是否透過命令行傳遞參數**
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 5) {
    cat("✅ 檢測到 4 個命令行參數:\n")
    cat("  1. pd_file_path:", args[1], "\n")
    cat("  2. idat_file_dir:", args[2], "\n")
    cat("  3. arraytype_file:", args[3], "\n")
    cat("  4. output_path:", args[4], "\n")
    cat("  5. cores:", args[5], "\n")
    # ✅ 如果有傳遞 4 個參數，則使用命令行參數
    pd_file_path <- args[1]
    idat_file_dir <- args[2]
    arraytype_file <- args[3]
    output_path <- args[4]
    cores <- as.integer(args[5])  # 設定 cores 數量
    # ✅ gse從 `pd_file_path` 推測
    gse <- strsplit(basename(pd_file_path), "_")[[1]][1]  
} else {
    stop("usage: Rscript champ_idat2norm.R pd_file_path idat_file_dir arraytype_file output_path cores")
}

# ✅ **從文件中獲取 arraytype**
if (file.exists(arraytype_file)) {
    arraytype <- readLines(arraytype_file, warn = FALSE)
    arraytype <- trimws(arraytype)  # 去除首尾空白字符
    cat("📄 從 arraytype.txt 獲取 arraytype:", arraytype, "\n")
} else {
    stop("❌ 找不到 arraytype.txt 文件")
}
# ✅ **檢查arraytype是否正確**
if (arraytype == "" || !(arraytype %in% c("EPICv1", "EPICv2", "450K"))) {
  stop(paste("❌ 無效的 arraytype:", arraytype, "請確認 arraytype 是否正確，應為 EPICv1, EPICv2, 或 450K"))
}

# ✅ **檢查 cores 是否有效**
if (is.na(cores) || cores <= 0) {
    cores <- 1  # 預設為1 核心
    cat("⚠️ 核心數無效，預設為:", cores, "\n")
}

cat("🚀 開始處理 GSE:", gse, "\n")
cat("📄 Sample Sheet:", pd_file_path, "\n")
cat("📂 IDAT 目錄:", idat_file_dir, "\n")
cat("🔬 Array Type:", arraytype, "\n")
cat("📂 輸出檔案:", output_path, "\n")
cat("🔥 使用核心數:", cores, "\n")


# ✅ **載入必要的 library**
suppressMessages(library(ChAMP))
suppressMessages(library(data.table))

# ✅ **定義處理函數**
process_idat <- function(pd_file_path, idat_file_dir, arraytype, output_path, cores) {
    tryCatch({
        # ✅ 檢查文件是否存在
        if (!file.exists(pd_file_path)) {
            stop(paste("❌ 找不到 Sample Sheet:", pd_file_path))
        }
        if (!dir.exists(idat_file_dir)) {
            stop(paste("❌ 找不到 IDAT 目錄:", idat_file_dir))
        }

        # ✅ 讀取樣本表
        pd_file <- fread(pd_file_path, header = TRUE, sep = ",")

        # ✅ 加載 IDAT 數據
        myLoad <- champ.load(directory = idat_file_dir, arraytype = arraytype)

        # ✅ 正規化數據
        myNorm <- champ.norm(beta = myLoad$beta, arraytype = arraytype, cores = cores)

        # ✅ 釋放記憶體
        rm(myLoad)
        gc()

        # ✅ 轉換為 data.table 格式，並添加 probeID 作為第一列
        myNorm_dt <- data.table(probeID = rownames(myNorm), myNorm)
        setorder(myNorm_dt, probeID)

        # ✅ 高效寫入 CSV
        fwrite(myNorm_dt, output_path, row.names = FALSE)

        cat("✅ 處理完成，輸出儲存於:", output_path, "\n")

    }, error = function(e) {
        cat("❌ 發生錯誤:", e$message, "\n")
    })
}

# ✅ 執行 IDAT 處理
process_idat(pd_file_path, idat_file_dir, arraytype, output_path, cores)
