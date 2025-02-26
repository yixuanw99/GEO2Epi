import argparse
import pandas as pd
from biolearn.model_gallery import ModelGallery
from biolearn.data_library import GeoData

def run_biolearn(beta_matrix_path, sample_sheet_path, output_path):
    
    gallery = ModelGallery()
    
    # Load data
    data = GeoData.from_methylation_matrix(beta_matrix_path)
    sample_sheet = pd.read_csv(sample_sheet_path)
    data.metadata = sample_sheet.set_index('Sample_Name')
    data.metadata.index.name = None
    data.metadata.rename(columns={'Gender': 'sex', 'Age': 'age'}, inplace=True)
    if 'sex' in data.metadata.columns:
        data.metadata['sex'] = data.metadata['sex'].map({'Female': 2, 'Male': 1})
    
    models = [
        "Horvathv1", "Hannum", "Lin", "PhenoAge", "YingCausAge", "YingDamAge", "YingAdaptAge", "Horvathv2", "PEDBE", 
        "Zhang_10", "DunedinPoAm38", "DunedinPACE", "GrimAgeV1", "GrimAgeV2", "AlcoholMcCartney", "BMI_McCartney", 
        "DNAmTL", "HRSInCHPhenoAge", "Knight", "LeeControl", "SexEstimation", "LeeRefinedRobust", "LeeRobust", 
        "SmokingMcCartney", "DownSyndrome", "EducationMcCartney", "TotalCholesterolMcCartney", 
        "HDLCholesterolMcCartney", "LDLCholesterolMcCartney", "BodyFatMcCartney", "BMI_Reed", "ProstateCancerKirby", 
        "HepatoXu", "CVD_Westerman", "AD_Bahado-Singh", "DepressionBarbu"
    ]
    
    # Validate metadata
    if 'age' not in data.metadata.columns or 'sex' not in data.metadata.columns:
        models.remove("GrimAgeV1")
        models.remove("GrimAgeV2")
    
    results = []
    
    # Run models
    for i, model in enumerate(models, start=1):
        print(f"Processing model {i}/{len(models)}: {model}")
        try:
            df = gallery.get(model).predict(data)
            df = df.add_prefix(f"{model}_")
            results.append(df)
        except Exception as e:
            print(f"⚠️ Model {model} failed: {e}")
    
    # Save results
    if results:
        combined_df = pd.concat(results, axis=1)
        combined_df = pd.merge(combined_df, sample_sheet, left_index=True, right_on='Sample_Name')
        combined_df.to_csv(output_path, index=True)
    
    print(f"✅ BioLearn Analysis Complete. Results saved at: {output_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run BioLearn Analysis")
    parser.add_argument("--beta_matrix", type=str, help="Path to Beta Matrix CSV", required=False)
    parser.add_argument("--sample_sheet", type=str, help="Path to Sample Sheet CSV", required=False)
    parser.add_argument("--output_path", type=str, help="Path to save BioLearn output", required=False)
    args = parser.parse_args()
    
    run_biolearn(args.beta_matrix, args.sample_sheet, args.output_path)