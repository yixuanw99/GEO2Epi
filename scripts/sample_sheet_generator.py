# -*- coding: utf-8 -*-

import argparse
import xml.etree.ElementTree as ET
import pandas as pd
import re

def process_gse_xml(xml_file, arraytype_outpath, sample_sheet_outpath):
    """解析 GSE XML 並產生 Sample Sheet"""
    
        # 解析 XML
    tree = ET.parse(xml_file)
    root = tree.getroot()
    
    # 初始化資料存儲
    data = []
    characteristics_labels = set()
    ns = {'miniml': 'http://www.ncbi.nlm.nih.gov/geo/info/MINiML'}

    # 設定有效的平台
    valid_platforms = {
        "GPL13534": "450K",
        "GPL16304": "450K",
        "GPL21145": "EPICv1",
        "GPL23976": "EPICv1",
        "GPL33022": "EPICv2"
    }

    platform = ""

    # 解析 Sample 元素
    for sample in root.findall("miniml:Sample", namespaces=ns):
        sample_iid = sample.attrib['iid']
        if sample.find("miniml:Type", namespaces=ns).text != "genomic":
            print(f"⚠️ 跳過非基因樣本: {sample_iid}")
            continue
        
        sample_accession = sample.find("miniml:Accession[@database='GEO']", namespaces=ns).text
        if sample_accession != sample_iid:
            raise ValueError(f"❌ 樣本 ID 不匹配: {sample_accession} vs {sample_iid}")
        
        sample_platform = sample.find("miniml:Platform-Ref", namespaces=ns).attrib['ref']
        if sample_platform not in valid_platforms:
            print(f"⚠️ 不支援的平台: {sample_platform}，跳過樣本 {sample_accession}")
            continue
        
        if not platform:
            platform = sample_platform
        elif sample_platform != platform:
            raise ValueError(f"❌ 發現多個不同的平台: {platform} vs {sample_platform}")

        # 解析 IDAT 文件名稱
        idat_info = sample.find("miniml:Supplementary-Data[@type='IDAT']", namespaces=ns)
        if idat_info is not None:
            file_name = idat_info.text.strip().split("/")[-1]
            parts = file_name.split("_")
            if len(parts) >= 3:
                extracted_sample_name = parts[0]
                sentrix_id = parts[1]
                sentrix_position = parts[2].split(".")[0]
                if extracted_sample_name != sample_accession:
                    raise ValueError(f"❌ 樣本名稱不匹配: {sample_accession} vs {extracted_sample_name} in {file_name}")
            else:
                raise ValueError(f"❌ IDAT 檔案名稱格式錯誤: {file_name}")
        else:
            sentrix_id, sentrix_position = "N/A", "N/A"

        # 提取 Characteristics
        sample_channel = sample.find("miniml:Channel", namespaces=ns)
        characteristics = {}
        for char in sample_channel.findall("miniml:Characteristics", namespaces=ns):
            if char is None:
                print(f"⚠️ 樣本 {sample_accession} 無 Characteristics")
            tag = re.sub(r"[ -]", "_", char.attrib['tag'])
            column_name_mapping = {"age": "Age", "Sex": "Gender", "sex": "Gender", "gender": "Gender"}
            tag = column_name_mapping.get(tag, tag)
            value = re.sub(r"[ -]", "_", char.text.strip()) if char.text else "N/A"
            characteristics[tag] = value
            characteristics_labels.add(tag)

        # 儲存資料
        row = {"Sample_Name": sample_accession, "Sentrix_ID": sentrix_id, "Sentrix_Position": sentrix_position, **characteristics}
        data.append(row)

    # 轉換為 DataFrame 並存儲
    df = pd.DataFrame(data)
    ordered_columns = ["Sample_Name", "Sentrix_ID", "Sentrix_Position"] + sorted(characteristics_labels)
    df = df.reindex(columns=ordered_columns, fill_value="N/A")
    df.to_csv(sample_sheet_outpath, index=False)

    print(f"✅ Sample Sheet 生成成功: {sample_sheet_outpath}")

    # 取得對應的 `arraytype`
    arraytype = valid_platforms.get(platform, "UNKNOWN")

    # ✅ 存入 `arraytype.txt`
    with open(arraytype_outpath, "w") as f:
        f.write(arraytype)

    print(f"✅ Array type: {arraytype} (已寫入 {arraytype_outpath})")

# **主程式**
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process GSE XML to generate Sample Sheet")
    parser.add_argument("--xml_file", type=str, required=False, help="Path to GSE XML file")
    parser.add_argument("--arraytype_outpath", type=str, required=False, help="Path to save arraytype")
    parser.add_argument("--sample_sheet_outpath", type=str, required=False, help="Path to save sample sheet")

    args = parser.parse_args()

    xml_file = args.xml_file
    arraytype_outpath = args.arraytype_outpath
    sample_sheet_outpath = args.sample_sheet_outpath

    process_gse_xml(xml_file, arraytype_outpath, sample_sheet_outpath)
