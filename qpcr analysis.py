import pandas as pd
import openpyxl
import numpy as np
   
def trunc(plate, name):
    """
    remove unnecessary headers of each raw data table until the first row starts with your chosen name
    re-define the header using the line that your chosen name first appears 

    Args:
    - `plate` (pd.dataframe): raw data
    - `name` (str): a word which the header row starts

    Returns:
    - `new_table`(pd.dataframe): dataframe with redefined header
    """
    mask = plate.iloc[:, 0].str.lower() == name
    line = mask.idxmax() if mask.any() else None
    plate = plate.iloc[line:]
    old_column_names = plate.columns.tolist()
    column_names = [plate[name].iloc[0] for name in old_column_names]
    new_table = {}
    for (oname, name) in zip(old_column_names, column_names):
        new_table[name] = list(plate[oname].iloc[1:])
    return pd.DataFrame.from_dict(new_table)


def cleanup(plate, colmean):
    """
    for the input plate data, clean up the input column by converting all values in the column to float
    the value is set to 0 if it is NA in the original data

    Args:
    - `plate` (pd.dataframe): qPCR raw data
    - `colmean` (str): column name of the column needing to be clean up

    Returns:
    - `new_plate`(pd.dataframe): clean table
    """

    new_plate = plate.loc[:, ['Well Position', colmean]]
    new_plate.loc[:, colmean] = new_plate.loc[:, colmean].apply(pd.to_numeric, errors='coerce')
    new_plate.loc[:, colmean] = new_plate.loc[:, colmean].fillna(0)
    return new_plate


# add qPCR Quantity Mean of housekeeping gene and gene of interest to output table
def matchmaker(df, plates):
    '''
    split the key file to match with the plate setup
    match quantity means to well position in the output file

    Args:
    - `split` (pd.dataframe): the dataframe contaning plate number keys to be edited
    - `plates`(pd.dataframe): the qPCR results plate

    Returns:
    - `split_dfs`(pd.dataframe): dataframe containing HK gene Quantity Mean and GI Quantity Mean matched by 
    PCR well setup from PCR raw data
'''

    row_numbers = [index for index, row in key.iterrows() if 'A3' in row['qPCR Well']]
    split_dfs = []
    for i in range(len(row_numbers) - 1):
        split_dfs.append(df[row_numbers[i]: row_numbers[i + 1]].reset_index(drop = True))

    split_dfs.append(df[row_numbers[-1]:].reset_index(drop = True))

    for i in range(len(split_dfs)):
        hk_key = f'hk_plate{i}'
        gi_key = f'gi_plate{i}'
        df_hk = plates[hk_key]
        df_gi = plates[gi_key]
        split_dfs[i] = merge(split_dfs[i], df_hk, "Quantity Mean", "Quantity_Mean_HK")
        split_dfs[i] = merge(split_dfs[i], df_gi, "Quantity Mean", "Quantity_Mean_GI")
        #print(split_dfs[i])
        
    return split_dfs

def merge(split, plate, src_col, tgt_col):
    '''
    match quantity means to well position in the output file

    Args:
    - `split` (pd.dataframe): the dataframe contaning plate number keys to be edited
    - `plate`(pd.dataframe): the qPCR results plate
    - `scr_col`(str): the column name of the qPCR results plate contaning Quantity Mean
    - `tgt_col`(str): name of the target column to put the Quantity Mean in the output dataframe

    Returns:
    - `split`(pd.dataframe): dataframe containing HK gene Quantity Mean and GI Quantity Mean matched by 
    PCR well setup from PCR raw data
    '''
    split[tgt_col] = split[tgt_col].astype(float)

    for i in range(split.shape[0]):
        for j in range(plate.shape[0]):
            if plate.loc[j, 'Well Position'] in split.loc[i, 'qPCR Well']:
                split.loc[i, tgt_col] = max(float(plate.loc[j, src_col]), float(split.loc[i, tgt_col]))

    return split 

# 

def calculator(df):
    '''
    calculate quantity mean of house keeping gene by tissue type, calculate z = GI/(HK/HK by Tissue)
    calculate z*10
    calculate log10(z)
    if z ==0, z is set to 9 for enabling taking the log
    if HK==0, z is set to NaN

    Args:
    - `df` (pd.dataframe): the dataframe contaning tissue type, HK, GI quantity mean

    Returns:
    - `df`(pd.dataframe): dataframe containing calculated results in new columns
    '''
    
    mean_hk_by_tissue = df.groupby('TISSUE')['Quantity_Mean_HK'].mean().reset_index()
    df = df.merge(mean_hk_by_tissue, on='TISSUE', suffixes=('', '_TissueMean'))
    df['z'] = df.apply(lambda row: 9 if row['Quantity_Mean_GI'] == 0 else 
                                  (np.nan if row['Quantity_Mean_HK'] == 0 else 
                                   row['Quantity_Mean_GI'] / (row['Quantity_Mean_HK'] / row['Quantity_Mean_HK_TissueMean'])), axis=1)

    df['z * 10'] = df['z'] * 10

    df['log(z)'] = df['z * 10'].apply(lambda x: np.log10(x) if x > 0 else np.nan)
    
    return df        


if __name__ == "__main__":
   
    ####################
    # input files here #
    ####################

    # input file path to house keeping gene plates in order
    hk_plates = [
        "C:/Users/kingh/OneDrive/Desktop/QPCR/2021-02-21_LNT020_RPS29_Plate1.xls",
        "C:/Users/kingh/OneDrive/Desktop/QPCR/2021-02-18_LNT020_RPS29_plate2.xls",
        "C:/Users/kingh/OneDrive/Desktop/QPCR/2021-02-21_LNT020_RPS29_plate3.xls",
    ]

    # input file path to gene of interest plates in order
    gi_plates = [
        "C:/Users/kingh/OneDrive/Desktop/QPCR/2021-02-18_LNT020_cr6_plate1.xls",
        "C:/Users/kingh/OneDrive/Desktop/QPCR/2021-02-21_LNT020_cr6_plate2.xls",
        "C:/Users/kingh/OneDrive/Desktop/QPCR/2021-02-18_LNT020_cr6_plate3.xls",
    ]

    # Import the qPCR plate setup key
    # paste the path to your qPCR plate setup key after pd.read_excel
    key = pd.read_excel("C:/Users/kingh/OneDrive/Desktop/QPCR/LN-T020 - LN051 tissue.xlsx")

    # input a sheet that contains genotype, virus and other information
    info = pd.read_excel("C:/Users/kingh/OneDrive/Desktop/QPCR/LN051_information.xlsx", header=0)

    assert len(hk_plates) == len(gi_plates), "hk and gi should have equal number files"
    n = len(hk_plates)

    plates = {
        f"hk_plate{i}": pd.read_excel(path, sheet_name='Results') for i, path in enumerate(hk_plates)
    }

    plates.update({
        f"gi_plate{i}": pd.read_excel(path, sheet_name='Results') for i, path in enumerate(gi_plates)
    })

    for k, v in plates.items():
        plates[k] = trunc(v, 'well')
        plates[k] = cleanup(plates[k],'Quantity Mean')

    

# modify "key" to separate the list of qPCR Well and drop unnessary columns
    key = trunc(key, 'well')
    key['qPCR Well'] = key['qPCR Well'].str.split(' & ')
    key.drop(columns = ['Well', '#'], axis = 1, inplace = True)
    key.rename(columns={'SAMPLE #': 'Sample'}, inplace= True)

    # add genotype, virus... information to the plate set-up key
    output = key.copy()
    output['Quantity_Mean_HK'] = 0
    output['Quantity_Mean_GI'] = 0
    output['z'] = 0
    output['z * 10'] = 0
    output['log(z)'] = 0
    new_info = info.iloc[:, [0, 2, 8]]
    new_info = new_info.rename(columns={'No': 'Sample'})
    output = pd.merge(output, new_info, on='Sample', how='left')

    # add quantity mean
    output = matchmaker(output, plates)
    # do calculations
    output = pd.concat(output, ignore_index= True)
    output = calculator(output)

    #############
    # output
    #############

    # name your desired output file with name.xlsx
    output_name = "C:/Users/kingh/OneDrive/Desktop/LN051_qPCR_analysis.xlsx"
    output.to_excel(output_name, sheet_name="Sheet1", index=False)



