import sys, os
import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, rdFMCS, Descriptors, PandasTools,SaltRemover
remover = SaltRemover.SaltRemover()

def GetLarge(string):
    '''Get largest string'''
    List = string.split('.')
    List = sorted(List, key=lambda x:len(x), reverse=True)
    return List[0]

def GetSimilarityScore(query, fps):
    '''Calcluate maximum similarity between query and known binders'''
    query_fp = AllChem.GetMorganFingerprint(Chem.MolFromSmiles(query),3)
    Similarity = [DataStructs.TanimotoSimilarity(x,query_fp) for x in fps]
    Index = np.argmax(Similarity)
    return (Similarity[Index], Index)

def apply_and_concat(dataframe, field, func, column_names, arg):
    '''Incorporate similarity and cluster'''
    return pd.concat((
        dataframe,
        dataframe[field].apply(lambda x: pd.Series(func(x, arg), index=column_names))), axis=1)

def Run(df, temp_fps):
    Check = []
    new_smile = []

    for i in df.smiles.tolist():
        try:
            mol = Chem.MolFromSmiles(GetLarge(i))
            if mol is None:
                Check.append(1)
                new_smile.append(0)
            else:
                Check.append(2)
                mol = remover.StripMol(mol, dontRemoveEverything=True)
                new_smile.append(GetLarge(Chem.MolToSmiles(mol)))
        except:
            Check.append(0)
            new_smile.append(0)

    df['Check'] = pd.Series(Check, index=df.index)
    df['new_smile'] = pd.Series(new_smile, index=df.index)
    df = df.loc[df['Check']==2]
    df = df[['zinc_id','new_smile']]
    df = df.drop_duplicates('new_smile')
    df = apply_and_concat(df, 'new_smile', GetSimilarityScore, ['tanimoto_sim','cluster'], temp_fps)
    df['tanimoto_sim'] = round(df['tanimoto_sim'],3)
    df['cluster'] = df['cluster'].astype(int)
    
    return df 

def main():
    args = sys.argv[1:]
    if not args:
        print ('usage: python cal_similar.py input output')

        sys.exit(1)

    elif sys.argv[1] == '--help':
        print ('usage: python cal_similar.py input output')

        sys.exit(1)

    elif sys.argv[1].endswith('.csv'):
        if len(args) == 2:
            df = pd.read_csv(sys.argv[1])
            temp = pd.read_csv('/home/cyang/random/VantAl/dataset/curated_binders/binders_after_cluster.csv')
            fps = [AllChem.GetMorganFingerprint(Chem.MolFromSmiles(x), 3) for x in temp.smiles.tolist()]
            df = Run(df, fps)
            df.to_csv(sys.argv[2], index=False)
            
        else:
            sys.exit(1)
    else:
        sys.exit(1)

if __name__ == '__main__':
    main()
