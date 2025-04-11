import csv
from io import StringIO
from rdkit.Chem import MolFromSmiles
from dev_tox.dev_tox import main


def get_csv_from_smiles(smiles_list, options):
    # CSV writer expects a file object, not a string. 
    # StringIO can be used to store a string as a file-like object.

    options["make_prop_img"] = False  # do not need to create images for csv

    headers = []
    for key, val in options.items():
        if (key not in ["calculate_ad", "make_prop_img"]) and val:
            headers.extend([key, key + " Probability", key + " Threshold"])
            if options["calculate_ad"]:
                headers.append(key + " in AD")

    string_file = StringIO()
    writer = csv.DictWriter(string_file, fieldnames=['SMILES', *headers])
    writer.writeheader()

    for smiles in smiles_list:
        molecule = MolFromSmiles(smiles)

        row = {'SMILES': smiles}

        if molecule is None:
            row['SMILES'] = f"(invalid){smiles}"
            writer.writerow(row)
            continue

        data = main(smiles, **options)

        for model_name, pred, pred_proba, thresh, ad, _ in data:
            try:
                row[model_name] = pred
                row[model_name+" Probability"] = pred_proba
                row[model_name+" Threshold"] = thresh
                if options["calculate_ad"]:
                    row[model_name + " in AD"] = ad
            except ValueError:
                row[model_name] = "No prediction"  # if pred_proba is string skip
                row[model_name + " Probability"] = "NA"
                row[model_name + " Threshold"] = "NA"
                if options["calculate_ad"]:
                    row[model_name + " in AD"] = "NA"

        writer.writerow(row)

    return string_file.getvalue()
