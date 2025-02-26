#!/bin/bash
# extract_raw.sh
# 用途：解壓縮 RAW tar 檔到原始目錄，然後將所有 .idat.gz 檔案搬移到指定目錄並解壓

set -euo pipefail

if [ "$#" -lt 2 ]; then
    echo "Usage: $0 <tar_file> <outdir>"
    exit 1
fi

tar_file="$1"
outdir="$2"
tar_dir=$(dirname "$tar_file")

if [ -f "$tar_file" ]; then
    echo "找到 $tar_file, 開始解壓縮到 $tar_dir..."
    tar -xf "$tar_file" -C "$tar_dir"
else
    echo "❌ 找不到 $tar_file"
    exit 1
fi

echo "搜尋 $tar_dir 中的 .idat.gz 檔案..."

if compgen -G "${tar_dir}/*.idat.gz" > /dev/null; then
    echo "已找到 .idat.gz 檔案，開始移動至 $outdir 並執行 gunzip..."
    mkdir -p "$outdir"
    # 將所有找到的 .idat.gz 檔案搬移到 outdir
    mv "$tar_dir"/*.idat.gz "$outdir"
    gunzip "$outdir"/*.idat.gz
    echo "完成，檔案已解壓並保存於: $outdir"
else
    echo "未找到 .idat.gz 檔案"
fi
echo "--------------------------------------"
