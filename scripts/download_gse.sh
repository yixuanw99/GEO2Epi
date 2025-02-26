#!/bin/bash
# download_gse.sh
# 用途：根據傳入的 gse_id 下載檔案，並將檔案儲存到指定的下載目錄

if [ "$#" -lt 2 ]; then
    echo "Usage: $0 <gse_id> <download_directory>"
    exit 1
fi

gse_id="$1"
output_dir="$2"

# 若目錄不存在則建立
mkdir -p "$output_dir"

# 移除 "GSE" 前綴以取得純數字
num=${gse_id#GSE}

# 計算需要的前綴長度
prefix_len=$((${#num} - 3))
prefix=${gse_id:0:$((prefix_len + 3))}

# 構建 FTP 路徑（改為目錄）
miniml_ftp_path="ftp://ftp.ncbi.nlm.nih.gov/geo/series/${prefix}nnn/$gse_id/miniml/${gse_id}_family.xml.tgz"
suppl_ftp_path="ftp://ftp.ncbi.nlm.nih.gov/geo/series/${prefix}nnn/$gse_id/suppl/${gse_id}_RAW.tar"

echo "======================================"
echo "開始下載: $gse_id"
echo "下載miniml: $miniml_ftp_path"
echo "下載suppl: $suppl_ftp_path"
# 進入指定的輸出目錄進行下載
cd "$output_dir" || exit 1

# 設定 User-Agent 避免被封鎖
USER_AGENT="Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36"

# 下載 miniml，啟用續傳、避免重複下載、限制速率
wget --wait=5 --timeout=30 --user-agent="$USER_AGENT" \
     --timestamping --continue -nd -np -v --show-progress "$miniml_ftp_path" || {
    echo "下載 miniml 失敗，錯誤碼: $?"; exit 1;
}

# 下載 suppl，啟用續傳、避免重複下載、限制速率
wget --wait=5 --timeout=30 --user-agent="$USER_AGENT" \
     --timestamping --continue -nd -np -v --show-progress "$suppl_ftp_path" || {
    echo "下載 suppl 失敗，錯誤碼: $?"; exit 1;
}

echo "下載時間: $(date '+%Y-%m-%d %H:%M:%S')"
echo "下載檔案儲存至: $output_dir"
echo "======================================"