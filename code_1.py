from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import SimilarityMaps
from scipy.spatial.distance import cdist
import numpy as np

import glob
import gzip
import bz2
import os
import _pickle as cPickle

import io
import matplotlib.pyplot as plt

# god hates me so in my version of python I cannot supress these damn user warning so I do this nuclear option instead
"""import warnings
def warn(*args, **kwargs):
    pass
warnings.warn = warn"""

MODEL_DICT = {
    'Overall Toxicity': ['DT_overall_model.joblib'],
    'First Trimester Toxicity': ['DT_first_trimester_model.joblib'],
    'Second Trimester Toxicity': ['DT_second_trimester_model.joblib'],
    'Third Trimester Toxicity': ['DT_third_trimester_model.joblib'],
}

# lol I'm just like screw code readability sorry
MODEL_DICT_INVERT = {v: key for key, val in MODEL_DICT.items() for v in val}

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
    },


AD_DICT = {
    True: "Inside",
    False: "Outside"
}


def run_prediction(model, model_data, smiles, calculate_ad=True):
    """_summary_

    Args:
        model (_type_): _description_
        model_data (_type_): _description_
        smiles (_type_): _description_
        calculate_ad (bool, optional): _description_. Defaults to True.

    Returns:
        _type_: _description_
    """
    fp = np.zeros((2048, 1))
    _fp = AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(smiles), radius=3, nBits=2048)
    DataStructs.ConvertToNumpyArray(_fp, fp)

    pred_proba = model.predict_proba(fp.reshape(1, -1))[:, 1]
    pred = 1 if pred_proba > 0.5 else 0

    if pred == 0:
        pred_proba = 1-pred_proba

    if calculate_ad:
        ad = model_data["D_cutoff"] > np.min(cdist(model_data['Descriptors'].to_numpy(), fp.reshape(1, -1)))
        return pred, pred_proba, ad
    return pred, pred_proba, None


def get_prob_map(model, smiles):
    def get_fp(mol, idx):
        fps = np.zeros((2048, 1))
        _fps = SimilarityMaps.GetMorganFingerprint(mol, idx, radius=3, nBits=2048)
        DataStructs.ConvertToNumpyArray(_fps, fps)
        return fps

    def get_proba(fps):
        return float(model.predict_proba(fps.reshape(1, -1))[:, 1])

    mol = Chem.MolFromSmiles(smiles)
    fig, _ = SimilarityMaps.GetSimilarityMapForModel(mol, get_fp, get_proba)
    imgdata = io.StringIO()
    fig.savefig(imgdata, format='svg')
    imgdata.seek(0)  # rewind the data
    plt.savefig(imgdata, format="svg", bbox_inches="tight")

    return imgdata.getvalue()


def multiclass_ranking(ordered_preds):
    idx = 0
    one_detected = False
    for i, o in enumerate(ordered_preds):
        if int(o) == 1:
            if not one_detected:
                idx = i+1
                one_detected = True
        if int(o) == 0:
            if one_detected:
                idx = 0
                return idx
    return idx if idx != 0 else len(ordered_preds)+1


def main(smiles, calculate_ad=True, make_prop_img=False, **kwargs):

    print(smiles)

    def default(key, d):
        if key in d.keys():
            return d[key]
        else:
            return False

    models = sorted([f for f in glob.glob("./PhaKinPro/models/*.joblib")], key=lambda x: x.split("_")[1])
    models_data = sorted([f for f in glob.glob("./PhaKinPro/models/*.joblib")], key=lambda x: x.split("_")[1])

    values = {}

    for model_endpoint, model_data in zip(models, models_data):
        if not default(MODEL_DICT_INVERT[os.path.basename(model_endpoint)], kwargs):
            continue
        with gzip.open(model_endpoint, 'rb') as f:
            model = cPickle.load(f)

        with bz2.BZ2File(model_data, 'rb') as f:
            model_data = cPickle.load(f)

        pred, pred_proba, ad = run_prediction(model, model_data, smiles, calculate_ad=calculate_ad)

        svg_str = ""
        if make_prop_img:
            svg_str = get_prob_map(model, smiles)

        values.setdefault(MODEL_DICT_INVERT[os.path.basename(model_endpoint)], []).append([int(pred), str(round(float(pred_proba) * 100, 2)) + "%", AD_DICT[ad], svg_str])

    processed_results = []
    for key, val in values.items():
        if key in ['Overall Toxicity', 'First Trimester Toxicity', 'Second Trimester Toxicity', 'Third Trimester Toxicity']:
            new_pred = multiclass_ranking([_[0] for _ in val]) ### CHECK THIS IS CORRECT
            if new_pred == 0:
                processed_results.append([key, "Inconsistent result: no prediction", "Very unconfident", "NA", ""])
            else:
                # this is because of how the hierarchical model works
                if new_pred in [1, 2]:
                    p = 0
                else:
                    p = new_pred - 2
                processed_results.append([key, CLASSIFICATION_DICT[key][new_pred], val[p][1], val[p][2], val[p][3]])
        else:
            processed_results.append([key, CLASSIFICATION_DICT[key][val[0][0]], val[0][1], val[0][2], val[0][3]])

    return processed_results


def write_csv_file(smiles_list, calculate_ad=False):
    headers = list(MODEL_DICT.keys())

    if calculate_ad:
        headers = headers + [_ + "_AD" for _ in headers]

    string_file = StringIO()
    writer = csv.DictWriter(string_file, fieldnames=['SMILES', *headers])
    writer.writeheader()

    for smiles in tqdm(smiles_list):
        molecule = MolFromSmiles(smiles)

        row = {'SMILES': smiles}

        if molecule is None:
            row['SMILES'] = f"(invalid){smiles}"
            writer.writerow(row)
            continue

        data = main(smiles, calculate_ad=calculate_ad, **MODEL_DICT)

        for model_name, pred, pred_proba, ad, _ in data:
            try:
                pred_proba = float(pred_proba[:-1]) / 100  # covert back to 0-1 float
                row[
                    model_name] = pred_proba if pred == 1 else 1 - pred_proba  # this is to make sure its proba for class 1
            except ValueError:
                row[model_name] = "No prediction"  # if pred_proba is string skip
            if calculate_ad:
                row[model_name + "_AD"] = ad

        writer.writerow(row)

    return string_file.getvalue()


if __name__ == "__main__":
    import argparse
    import csv
    from io import StringIO
    from rdkit.Chem import MolFromSmiles
    from tqdm import tqdm

    parser = argparse.ArgumentParser()
    parser.add_argument("--infile", type=str, required=True,
                        help="location to csv of SMILES")
    parser.add_argument("--outfile", type=str, default=os.path.join(os.getcwd(), "phakin_output.csv"),
                        help="location and file name for output")
    parser.add_argument("--smiles_col", type=str, default="SMILES",
                        help="column name containing SMILES of interest"),
    parser.add_argument("--ad", action="store_true",
                        help="calculate the AD")
