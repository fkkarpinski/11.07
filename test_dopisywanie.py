import pandas as pd
import pubchempy as pcp
import requests
from tqdm import tqdm
import os

def getIDs(smiles):
    """get CIDs from SMILES"""
    ids = []
    for smile in tqdm(smiles, desc="getIDs()"):
        try:
            compound = pcp.get_compounds(smile, 'smiles')[0]
            ids.append(str(compound.cid))
        except Exception:
            ids.append("")
    return ids

def getNames(ids_lista):
    """get IUPAC name from CID, create a list with names"""
    names = []
    for cid in tqdm(ids_lista, desc="getNames()"):
        try:
            compound = pcp.Compound.from_cid(cid)
            names.append(compound.iupac_name)
        except Exception:
            names.append("")
    return names

def getAssayInfo(CID):
    """
    get data regarding AIDs, activities and BioAssay titles through an API connection,
    iterate through JSON rows to find desired values
    """
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{CID}/assaysummary/json"
    try:
        r = requests.get(url)
        if r.status_code == 200:
            data = r.json()
            aids, activities, titles = [], [], []
            rows = data.get("Table", {}).get("Row", [])
            for entry in rows:
                cells = entry.get("Cell", [])
                if cells and len(cells) > 0:
                    aids.append(cells[0])
                    activities.append(cells[4])
                    titles.append(cells[9])
            return aids, activities, titles
        else:
            return [], [], []
    except Exception:
        return [], [], []

def getAssayDescription(AID):
    """get BioAssay descriptions basing on AID, directly access JSON container with the description"""
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/aid/{AID}/description/JSON"
    try:
        r = requests.get(url)
        if r.status_code == 200:
            data = r.json()
            desc = data["PC_AssayContainer"][0]["assay"]["descr"]["description"]
            if isinstance(desc, list):
                return " ".join(str(d) for d in desc)
            return str(desc)
        else:
            return ""
    except Exception:
        return ""
    
def populate(smiles, ids_lista, names, output_file="partial_output.csv"):
    """populate the dataframe, consider a case when there are no details to include"""
    rows = []
    # Sprawdź, czy istnieje plik z częściowym outputem
    if os.path.exists(output_file):
        existing_data = pd.read_csv(output_file).to_dict(orient="records")
        rows.extend(existing_data)

    for i, (smile, cid, name) in enumerate(tqdm(zip(smiles, ids_lista, names), total=len(smiles), desc="populate()")):
        try:
            aids, activities, titles = getAssayInfo(cid)
            descriptions = []
            for aid in tqdm(aids, desc="getAssayDescription()", leave=False):
                descriptions.append(getAssayDescription(aid))
            if aids:
                for aid, activity, title, desc in zip(aids, activities, titles, descriptions):
                    rows.append({
                        "SMILES": smile,
                        "CID": cid,
                        "Name": name,
                        "AID": aid,
                        "ACTIVITY": activity,
                        "TITLE": title,
                        "DESCRIPTION": desc
                    })
            else:
                rows.append({
                    "SMILES": smile,
                    "CID": cid,
                    "Name": name,
                    "AID": "",
                    "ACTIVITY": "",
                    "TITLE": "",
                    "DESCRIPTION": ""
                })
        except Exception:
            rows.append({
                "SMILES": smile,
                "CID": cid,
                "Name": name,
                "AID": "",
                "ACTIVITY": "",
                "TITLE": "",
                "DESCRIPTION": ""
            })

        # Zapisywanie danych co 10 iteracji
        if i % 10 == 0 or i == len(smiles) - 1:
            pd.DataFrame(rows).to_csv(output_file, index=False, header=True)

    return rows

def toCSV(rows):
    tempToDataFrame = pd.DataFrame(rows)
    tempToDataFrame_unique = tempToDataFrame[
        (tempToDataFrame["AID"] == "") | (~tempToDataFrame["AID"].duplicated(keep="first"))
    ].reset_index(drop=True)
    
    tempToDataFrame_unique.to_csv("out.csv", index=False, header=True)

#################
#################

#define input file, column header "smiles"
inputFile = "cluster_x_0515_y_-30-39.csv" 
df_input = pd.read_csv(inputFile)

#choose rows to fill based on empty "target_seq_label_1" column
mask = df_input["target_seq_label_1"].isnull() | (df_input["target_seq_label_1"] == "")
df_missing = df_input[mask]

inputData = df_missing["SMILES"].tolist()
partial_output_file = "partial_output.csv"

#check for partial file
if os.path.exists(partial_output_file):
    processed_data = pd.read_csv(partial_output_file)
    processed_smiles = set(processed_data["SMILES"])
    inputData = [smile for smile in inputData if smile not in processed_smiles]
else:
    processed_data = pd.DataFrame()

#main funcs
IDlist = getIDs(inputData)
names = getNames(IDlist)
rows = populate(inputData, IDlist, names, output_file=partial_output_file)

final_rows = pd.concat([processed_data, pd.DataFrame(rows)], ignore_index=True)
toCSV(final_rows.to_dict(orient="records"))

df_filled = pd.read_csv("out.csv")

#merge input output files
result = pd.merge(df_input, df_filled, on="SMILES", how="left", suffixes=("", "_filled"))

# Zachowaj wszystkie kolumny z df_input, nawet jeśli nie ma danych w df_filled
result = result[[col for col in result.columns if not col.endswith("_filled")]]

#output final file
result.to_csv("pubchem_scraper_out.csv", index=False)