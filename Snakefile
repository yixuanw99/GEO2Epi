# Snakefile
import os

configfile: "config/config.yaml"

PROJECT_ROOT = config["PROJECT_ROOT"]
GSE_IDS      = config["GSE_ID_LIST"]
CORES        = config.get("CORES", 1)

# 定義目錄的函式（根據 gse_id 來構造路徑）
def raw_dir(gse_id):
    return os.path.join(PROJECT_ROOT, "data", "raw", gse_id)

def idat_dir(gse_id):
    return os.path.join(raw_dir(gse_id), "idat")

def processed_dir():
    return os.path.join(PROJECT_ROOT, "data", "processed")

def biolearn_dir():
    return os.path.join(PROJECT_ROOT, "data", "results", "biolearn")

def methyl_dir():
    return os.path.join(PROJECT_ROOT, "data", "results", "MethylDetectR")

rule all:
    input:
        expand(os.path.join(biolearn_dir(), "{gse_id}_biolearn_output.csv"), gse_id=GSE_IDS),
        expand(os.path.join(methyl_dir(), "{gse_id}_MethylDetectR_output.csv"), gse_id=GSE_IDS),
        expand(os.path.join(raw_dir("{gse_id}"), "raw_removed.txt"), gse_id=GSE_IDS)

# 針對全域資料夾 (processed, biolearn, methyl)
rule create_global_dirs:
    output:
        processed_dir = directory(processed_dir()),
        biolearn_dir = directory(biolearn_dir()),
        methyl_dir   = directory(methyl_dir())
    shell:
        """
        mkdir -p {output.processed_dir}
        mkdir -p {output.biolearn_dir}
        mkdir -p {output.methyl_dir}
        """

rule download_gse:
    output:
        xml_tgz = os.path.join(raw_dir("{gse_id}"), "{gse_id}_family.xml.tgz"),
        raw_tar = os.path.join(raw_dir("{gse_id}"), "{gse_id}_RAW.tar")
    params:
        gse_id = "{gse_id}",
        outdir = directory(raw_dir("{gse_id}"))
    log:
        "logs/{gse_id}/{gse_id}_download_gse.log"
    resources:
        download_jobs=1  # 限制最多執行 1 個
    shell:
        """
        bash {PROJECT_ROOT}/scripts/download_gse.sh {params.gse_id} {params.outdir} 2>&1 | tee -a {log}
        """

rule extract_miniml:
    input:
        xml_tgz = os.path.join(raw_dir("{gse_id}"), "{gse_id}_family.xml.tgz")
    output:
        xml = os.path.join(raw_dir("{gse_id}"), "{gse_id}_family.xml")
    log:
        "logs/{gse_id}/{gse_id}_extract_miniml.log"
    shell:
        """
        bash {PROJECT_ROOT}/scripts/extract_miniml.sh {input.xml_tgz} {output.xml} 2>&1 | tee -a {log}
        """

rule extract_raw:
    input:
        raw_tar = os.path.join(raw_dir("{gse_id}"), "{gse_id}_RAW.tar")
    output:
        extraction_done = os.path.join(raw_dir("{gse_id}"), "extraction_done.txt")
    params:
        idat_dir = directory(idat_dir("{gse_id}"))
    log:
        "logs/{gse_id}/{gse_id}_extract_raw.log"
    shell:
        """
        bash {PROJECT_ROOT}/scripts/extract_raw.sh {input.raw_tar} {params.idat_dir} 2>&1 | tee -a {log}
        touch {output.extraction_done}
        """

rule generate_sample_sheet:
    input:
        xml = os.path.join(raw_dir("{gse_id}"), "{gse_id}_family.xml")
    output:
        sample_sheet = os.path.join(idat_dir("{gse_id}"), "{gse_id}_Sample_Sheet.csv"),
        arraytype    = os.path.join(raw_dir("{gse_id}"), "arraytype.txt")
    log:
        "logs/{gse_id}/{gse_id}_generate_sample_sheet.log"
    shell:
        """
        python3 {PROJECT_ROOT}/scripts/sample_sheet_generator.py --xml_file {input.xml} --arraytype_outpath {output.arraytype} --sample_sheet_outpath {output.sample_sheet} 2>&1 | tee -a {log}
        """

rule champ_idat2norm:
    input:
        extraction_done = os.path.join(raw_dir("{gse_id}"), "extraction_done.txt"),
        sample_sheet = os.path.join(idat_dir("{gse_id}"), "{gse_id}_Sample_Sheet.csv"),
        arraytype    = os.path.join(raw_dir("{gse_id}"), "arraytype.txt")
    output:
        output_path = os.path.join(processed_dir(), "{gse_id}_champ_norm.csv")
    params:
        idat_dir = directory(idat_dir("{gse_id}")),
        cores = CORES
    log:
        "logs/{gse_id}/{gse_id}_champ_idat2norm.log"
    shell:
        """
        Rscript {PROJECT_ROOT}/scripts/champ_idat2norm.R {input.sample_sheet} {params.idat_dir} {input.arraytype} {output.output_path} {params.cores} 2>&1 | tee -a {log}
        """

rule clean_raw:
    input:
        extraction_done = os.path.join(raw_dir("{gse_id}"), "extraction_done.txt"),
        champ_done = os.path.join(processed_dir(), "{gse_id}_champ_norm.csv")
    output:
        raw_removed = os.path.join(raw_dir("{gse_id}"), "raw_removed.txt")
    params:
        raw_dir_path = lambda wildcards: raw_dir(wildcards.gse_id),
        idat_dir_path = lambda wildcards: idat_dir(wildcards.gse_id),
        gse_id = "{gse_id}"
    log:
        "logs/{gse_id}/{gse_id}_clean_raw.log"
    shell:
        """
        find {params.raw_dir_path} -maxdepth 1 ! -name "arraytype.txt" ! -name "*_family.xml" -type f -exec rm -f {{}} +
        find {params.idat_dir_path} -maxdepth 1 ! -name "{params.gse_id}_Sample_Sheet.csv" -type f -exec rm -f {{}} +
        echo "Cleaned raw data for {params.gse_id} - $(date)" | tee -a {log}
        touch {output.raw_removed}
        """

rule run_biolearn:
    input:
        beta_matrix  = os.path.join(processed_dir(), "{gse_id}_champ_norm.csv"),
        sample_sheet = os.path.join(idat_dir("{gse_id}"), "{gse_id}_Sample_Sheet.csv")
    output:
        output_path = os.path.join(biolearn_dir(), "{gse_id}_biolearn_output.csv")
    log:
        "logs/{gse_id}/{gse_id}_run_biolearn.log"
    shell:
        """
        python3 {PROJECT_ROOT}/scripts/biolearn_run.py --beta_matrix {input.beta_matrix} --sample_sheet {input.sample_sheet} --output_path {output.output_path} 2>&1 | tee -a {log}
        """

rule run_MethylDetectR:
    input:
        beta_matrix  = os.path.join(processed_dir(), "{gse_id}_champ_norm.csv"),
        sample_sheet = os.path.join(idat_dir("{gse_id}"), "{gse_id}_Sample_Sheet.csv")
    output:
        output_path = os.path.join(methyl_dir(), "{gse_id}_MethylDetectR_output.csv")
    params:
        project_root = PROJECT_ROOT
    log:
        "logs/{gse_id}/{gse_id}_run_MethylDetectR.log"
    shell:
        """
        Rscript {PROJECT_ROOT}/scripts/MethylDetectR_run.R {input.beta_matrix} {input.sample_sheet} {output.output_path} {params.project_root} 2>&1 | tee -a {log}
        """

