#!/bin/bash
# extract_miniml.sh
# 用途：解壓縮下載的 XML 檔案，並檢查是否成功

if [ "$#" -lt 2 ]; then
    echo "Usage: $0 <xml_tgz> <xml>"
    exit 1
fi

xml_tgz="$1"
xml="$2"
xml_dir=$(dirname "$xml")

echo -n "正在解壓縮... "
if [ -f "$xml_tgz" ]; then
    echo "找到 $xml_tgz, 開始解壓縮..."
    mkdir -p "$xml_dir"
    tar -xzf "$xml_tgz" -C "$xml_dir"
else
    echo "❌ 找不到 $xml_tgz"
    exit 1
fi

if [ -f "$xml" ]; then
    echo "檔案已保存到: $xml"
else
    echo "❌ 解壓縮失敗"
fi
echo "--------------------------------------"
