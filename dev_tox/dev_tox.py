#
import base64

from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Draw
from rdkit.Chem.Draw import SimilarityMaps
from scipy.spatial.distance import cdist
import numpy as np
from scipy.spatial.distance import pdist, cdist

"""import glob
import gzip
import bz2"""
import os
#import _pickle as cPickle

import io
import matplotlib.pyplot as plt


# for setting confidence-based AD:
AD_THRESH = 0.6


import joblib  # Ensure this import is at the beginning of your script

MODEL_DIR = os.path.join(os.path.dirname(__file__), "models")  # Directory where models are stored

MODEL_DICT = {
    'Overall Toxicity': [joblib.load(os.path.join(MODEL_DIR, 'ECFP4_overall_tox_svm.joblib'))],
    'First Trimester Toxicity': [joblib.load(os.path.join(MODEL_DIR, 'ECFP4_first_trimester_svm.joblib'))],
    'Second Trimester Toxicity': [joblib.load(os.path.join(MODEL_DIR, 'ECFP4_second_trimester_svm.joblib'))],
    'Third Trimester Toxicity': [joblib.load(os.path.join(MODEL_DIR, 'ECFP4_third_trimester_svm.joblib'))],
}


# CHECK THIS
AD_DICT = {
    'Overall Toxicity': [joblib.load(os.path.join(MODEL_DIR, 'overall_toxicity_AD.pkl'))],
    'First Trimester Toxicity': [joblib.load(os.path.join(MODEL_DIR, 'first_tri_toxicity_AD.pkl'))],
    'Second Trimester Toxicity': [joblib.load(os.path.join(MODEL_DIR, 'second_tri_toxicity_AD.pkl'))],
    'Third Trimester Toxicity': [joblib.load(os.path.join(MODEL_DIR, 'third_tri_toxicity_AD.pkl'))],
}


# # lol I'm just like screw code readability sorry
# MODEL_DICT_INVERT = {v["model"]: key for key, val in MODEL_DICT.items() for v in val}

CLASSIFICATION_DICT = {
    'Overall Toxicity': {
        0: "Non-toxic",
        1: "Toxic"
    },
    'First Trimester Toxicity': {
        0: "Non-toxic",
        1: "Toxic"
    },
    'Second Trimester Toxicity': {
        0: "Non-toxic",
        1: "Toxic"
    },
    'Third Trimester Toxicity': {
        0: "Non-toxic",
        1: "Toxic"
    },
    }


AD_DICT_BOOL = {
    True: "Inside",
    False: "Outside"
}


def _get_AD_thresh(training_smiles, file_name):
    fps = np.array([list(AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(smi), radius=2, nBits=2048, useFeatures=False))
                    for smi in training_smiles])

    dists = pdist(fps)
    mean_1 = dists.mean()
    dists_2 = dists[np.where(dists < mean_1)]
    mean_2 = dists_2.mean()
    std_2 = dists_2.std()

    threshold = mean_2 + (0.5 * std_2)

    import pickle
    pickle.dump((threshold, fps),  open(file_name, "wb"))


def calc_ad(query_fp, ad_tuple, selected_features):
    dist = cdist(query_fp.reshape(1, -1), ad_tuple[1][:, selected_features], "euclidean")
    return "TRUE" if bool((dist < ad_tuple[0]).any()) else "FALSE"


def run_prediction(model, smi, selected_features, calculate_ad=True, ad_tup=None, threshold=0.5):
    """_summary_

    Args:
        model (_type_): _description_
        model_data (_type_): _description_
        smi (_type_): _description_
        calculate_ad (bool, optional): _description_. Defaults to True.

    Returns:
        _type_: _description_
    """
    fp = np.zeros((2048, 1))
    # sub in your FP function
    _fp = AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(smi), radius=2, nBits=2048, useFeatures=False)
    DataStructs.ConvertToNumpyArray(_fp, fp)
    fp = fp[selected_features]

    pred_proba = model.predict_proba(fp.reshape(1, -1))[:, 1]
    pred = 1 if pred_proba > threshold else 0

    # if pred == 0:
    #     pred_proba = 1 - float(pred_proba)

    # used to get proba of the inactive class if deemed inactive
    # if pred == 0:
    #     pred_proba = 1-pred_proba

    if calculate_ad:
        ad = calc_ad(fp, ad_tup, selected_features)
        return pred, pred_proba, ad

    return pred, float(pred_proba), ""


def get_prob_map(model, smi, selected_features_file):
    """Generates a probability map for a given molecule.

    Args:
        model: The trained model.
        smi: The SMILES string for the molecule.
        selected_features_file (str): Path to the selected features file.

    Returns:
        SVG string representing the probability map.
    """
    # Load the selected features indices from the correct file
    selected_features = np.load(selected_features_file)

    def get_fp(mol, idx):
        """Generates the fingerprint for the molecule at a specific atom index."""
        fps = np.zeros((2048, 1))
        _fps = SimilarityMaps.GetMorganFingerprint(mol, idx, radius=3, nBits=2048)
        DataStructs.ConvertToNumpyArray(_fps, fps)

        # Ensure fps is reshaped to one-dimensional
        fps = fps.flatten()

        # Select only the features that were used during model training
        fps_selected = fps[selected_features]

        return fps_selected

    def get_proba(fps):
        """Gets the probability prediction for a specific fingerprint."""
        return float(model.predict_proba(fps.reshape(1, -1))[:, 1])

    mol = Chem.MolFromSmiles(smi)
    d2d = Draw.MolDraw2DCairo(500, 500)
    SimilarityMaps.GetSimilarityMapForModel(mol, get_fp, get_proba, draw2d=d2d)
    fig = base64.b64encode(d2d.GetDrawingText()).decode("ascii")
    return fig


def main(smi, calculate_ad=True, make_prop_img=False, **kwargs):
    values = {}
    trimester_toxic = -1  # Flag to check if any trimester toxicity is predicted as toxic

    for key, val in kwargs.items():
        if key in MODEL_DICT.keys():  # check if this kwarg is for a model
            if val:  # check if model is turned on
                model_data = MODEL_DICT[key][0]
                # print(f"Loading model from: {model_file}")
                
                # Load model and threshold from the joblib file
                # model_data = joblib.load(model_file)
                model = model_data['model']  # Extract model
                threshold = model_data.get('threshold', 0.5)  # Extract threshold or default to 0.5
                
                # Determine the correct selected features file based on the model type
                if key == 'Overall Toxicity':
                    selected_features_file = os.path.join(MODEL_DIR, 'selected_features_overall_tox_ecfp4.npy')
                elif key == 'First Trimester Toxicity':
                    selected_features_file = os.path.join(MODEL_DIR, 'selected_features_first_tri_ecfp4.npy')
                elif key == 'Second Trimester Toxicity':
                    selected_features_file = os.path.join(MODEL_DIR, 'selected_features_second_tri_ecfp4.npy')
                elif key == 'Third Trimester Toxicity':
                    selected_features_file = os.path.join(MODEL_DIR, 'selected_features_third_tri_ecfp4.npy')
                else:
                    raise ValueError(f"Unknown model type: {key}")
                selected_features = np.load(selected_features_file)

                print(threshold)
                # Run the prediction with the correct selected features file
                pred, pred_proba, ad = run_prediction(model, smi, selected_features, calculate_ad=calculate_ad, ad_tup=AD_DICT[key][0], threshold=threshold)

                # Generate probability map if needed
                contrib_base64_bytes = ""
                if make_prop_img:
                    contrib_base64_bytes = get_prob_map(model, smi, selected_features_file)

                pred_proba = round(float(pred_proba), 3)

                values[key] = [pred, float(pred_proba), round(threshold, 3), ad, contrib_base64_bytes]

                # Check if any trimester prediction is "Toxic"
                if key in ['First Trimester Toxicity', 'Second Trimester Toxicity', 'Third Trimester Toxicity'] and pred == 1:
                    if float(pred_proba) > trimester_toxic:
                        trimester_toxic = float(pred_proba)
                        trimester_thresh = threshold

    # Adjust "Overall Toxicity" prediction if any trimester is toxic
    if trimester_toxic != -1:
        # If overall toxicity has been predicted, update it to "Toxic"
        if 'Overall Toxicity' in values:
            values['Overall Toxicity'][0] = 1  # Set prediction to "Toxic"
            values['Overall Toxicity'][1] = trimester_toxic  # Adjust probability if needed
            values['Overall Toxicity'][2] = round(trimester_thresh, 3)  # Adjust probability if needed
        else:
            # If "Overall Toxicity" wasn't predicted earlier, create an entry for it as "Toxic"
            values['Overall Toxicity'] = [1, trimester_toxic, trimester_thresh, None, ""]  # Assuming full certainty if trimester is toxic

    # Prepare results for output
    processed_results = []
    for key, val in values.items():
        processed_results.append([key, CLASSIFICATION_DICT[key][val[0]], val[1], val[2], val[3], val[4]])

    return processed_results




# def write_csv_file(smiles_list, calculate_ad=False):
#     headers = list(MODEL_DICT.keys())
#
#     if calculate_ad:
#         headers = headers + [_ + "_AD" for _ in headers]
#
#     string_file = StringIO()
#     writer = csv.DictWriter(string_file, fieldnames=['SMILES', *headers])
#     writer.writeheader()
#
#     for smiles in tqdm(smiles_list):
#         molecule = MolFromSmiles(smiles)
#
#         row = {'SMILES': smiles}
#
#         if molecule is None:
#             row['SMILES'] = f"(invalid){smiles}"
#             writer.writerow(row)
#             continue
#
#         data = main(smiles, calculate_ad=calculate_ad, **MODEL_DICT)
#
#         for model_name, pred, pred_proba, ad, _ in data:
#             try:
#                 pred_proba = float(pred_proba[:-1]) / 100  # covert back to 0-1 float
#                 row[
#                     model_name] = pred_proba if pred == 1 else 1 - pred_proba  # this is to make sure its proba for class 1
#             except ValueError:
#                 row[model_name] = "No prediction"  # if pred_proba is string skip
#             if calculate_ad:
#                 row[model_name + "_AD"] = ad
#
#         writer.writerow(row)
#
#     return string_file.getvalue()


# if __name__ == "__main__":
#     import argparse
#     import csv
#     from io import StringIO
#     from rdkit.Chem import MolFromSmiles
#     from tqdm import tqdm
#
#     parser = argparse.ArgumentParser()
#     parser.add_argument("--infile", type=str, required=True,
#                         help="location to csv of SMILES")
#     parser.add_argument("--outfile", type=str, default=os.path.join(os.getcwd(), "phakin_output.csv"),
#                         help="location and file name for output")
#     parser.add_argument("--smiles_col", type=str, default="SMILES",
#                         help="column name containing SMILES of interest"),
#     parser.add_argument("--ad", action="store_true",
#                         help="calculate the AD")
