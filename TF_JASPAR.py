#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 21 21:48:31 2022

@author: alyion
"""
import coreapi
import numpy as np

##############################
#     Function功能設定區域     #
##############################

#輸入轉錄因子名稱，輸出JASPAR ID
def get_matrix_id(tf):
    # Initialize a client & load the schema document
    client = coreapi.Client()
    schema = client.get("https://jaspar.genereg.net/api/v1/docs/")
    # Interact with the API endpoint
    action = ["matrix", "list"]
    params = {
        "search": tf,
        "tax_id": "9606",
        "version": 'latest',
        "data_type": "ChIP-seq",
        "release": "2022"
    }
    result = client.action(schema, action, params=params)
    if tf.lower() in result['results'][0]['name'].lower() or result['results'][0]['name'].lower() in tf.lower():
        MOTIF = result['results'][0]['matrix_id']
        return MOTIF

#輸入JASPAR ID，輸出MOTIF
def get_matrix_ATCG(m_ID, Gfreq = 0.2, Afreq = 0.3, Cfreq = 0.2, Tfreq = 0.3): #預設MOTIF中ATCG的百分比切點為AT30%、CG20%
    # Initialize a client & load the schema document
    client = coreapi.Client()
    schema = client.get("https://jaspar.genereg.net/api/v1/docs/")
    # Interact with the API endpoint
    action = ["matrix", "read"]
    params = {
        "matrix_id": m_ID,
    }
    result = client.action(schema, action, params=params)
    ATCG_Matrix = dict(result['pfm'])
    GACT_array = np.append(np.array([ATCG_Matrix["G"],np.array(ATCG_Matrix["A"])]), np.array([ATCG_Matrix["C"],np.array(ATCG_Matrix["T"])]), axis=0)
    perc = GACT_array / GACT_array.sum(axis=0)
    #print(np.array([[f"{i:.2%}" for i in val] for val in perc]))
    seq1 = []
    for i in range(len(perc[0])):
        seq2 = []
        if perc[0][i] > Gfreq:
            seq2.append('G')
        if perc[1][i] > Afreq:
            seq2.append('A')
        if perc[2][i] > Cfreq:
            seq2.append('C')
        if perc[3][i] > Tfreq:
            seq2.append('T')
        seq2 = str("[" + "".join(seq2) + "]")
        seq1.append(seq2)
    seq1 = "".join(seq1).replace("[A]","A").replace("[T]","T").replace("[C]","C").replace("[G]","G")
    return seq1

def replace_values(ATCG):
    motif_sym_dict = {"[GA]":"R", "[CT]":"Y", "[AC]":"M", "[GT]":"K", "[GC]":"S", "[AT]":"W", "[ACT]":"H", "[GCT]":"B", "[GAC]":"V", "[GAT]":"D", "[GACT]":"N"}
    for motifs,sym in motif_sym_dict.items():
        ATCG = ATCG.replace(motifs,sym)
    return ATCG

#以轉錄因子nfkb確認api是否能撈取成功
get_matrix_id("nfkb")

##############################
#  Step.1 輸入TF 找JASPAR ID  #
##############################
inpath = "/Users/alyion/Google Drive/研究所/❖Su LAB❖/9.開會報告/20220913Meeting（DEGs）"
outpath = "/Users/alyion/Google Drive/研究所/❖Su LAB❖/9.開會報告/20220913Meeting（DEGs）"

####Step.1-1輸入檔案，製作TF串列

tf_list = []
path = inpath + '/total_TFs.txt'
with open(path) as f:
    for line in f.readlines():
        tf_list.append(line.rsplit())

print(tf_list)

tf2_list =[x for item in tf_list for x in item]

print(tf2_list)


####Step.1-2迴圈TF串列用get_matrix_id功能取得JASPAR ID，有對應的TF和ID存入字典，找不到的存入串列
#參數設定
tfYesmotif_dict = {}
tfNomotif_list = []
i = len(tf2_list)

for tfs in tf2_list:
    try:
        MOTIF = get_matrix_id(tfs)
        if MOTIF == None:
            print(f'Can not find JASPAR MOTIF ID for {tfs}!')
        else:
            tfYesmotif_dict.setdefault(tfs,MOTIF)
            print(f'The JASPAR MOTIF ID for {tfs} is {MOTIF}!')
    except:
        tfNomotif_list.append(tfs)
        print(f'Can not find JASPAR MOTIF ID for {tfs}!')
    i -= 1
    if i != 0:
        print(f'{i} TFs left!')
    else:
        print('Scraping done!!')
            
print(f'We found {len(tfYesmotif_dict)} MOTIFs from our uploaded TFs by scrapying JASPAR!')
print(tfYesmotif_dict)
print(f'{len(tfNomotif_list)} other TFs failed, listed below:')
print(tfNomotif_list)

####Step.1-3儲存為CSV檔案
with open(outpath + '/JASPAR_motif_ID.csv', 'w') as f:
    for key in tfYesmotif_dict.keys():
            f.write("%s,%s\n" % (key, tfYesmotif_dict[key]))
            
#########################
#       黏貼bed.url      #
#########################

tfmotifURL_dict = {}
for j in tfYesmotif_dict.keys():
    motif_id = tfYesmotif_dict[j]
    motif_url = 'https://jaspar.genereg.net/download/data/2022/bed/' + motif_id + '.bed'
    tfmotifURL_dict.setdefault(j, motif_url)
print(tfmotifURL_dict)
    
with open(outpath + '/JASPAR_motif_URL.csv', 'w') as f:
    for key in tfmotifURL_dict.keys():
        f.write("%s,%s\n" % (key, tfmotifURL_dict[key]))

#########################
#  Step.2 輸入ID 找MOTIF #
#########################
####Step.2-1輸入JASPAR ID，以功能get_matrix_ATCG迴圈輸出MOTIF序列，並將對應的ID與MOTIF存入字典
#參數設定
w = len(tfYesmotif_dict)
MOTIFSeq_dict = {}

for J_ID in tfYesmotif_dict.values():
    motifseq = get_matrix_ATCG(J_ID)
    print(motifseq)
    MOTIFSeq_dict.setdefault(J_ID, motifseq)
    w -= 1
    if w != 0:
        print(f'{w} MOTIFs left!')
    else:
        print('All done!')
print('After setting AGCT threshold to 30% and 20%, the result of new MOTIFs are shown below:')
print(MOTIFSeq_dict)

####Step.2-2儲存為CSV檔案
with open(outpath +'/JASPAR_motif2030.csv', 'w') as f:
    for key in MOTIFSeq_dict.keys():
        MOTIFSeq_dict[key] = "".join(MOTIFSeq_dict[key]).replace("[A]","A").replace("[T]","T").replace("[C]","C").replace("[G]","G")
        f.write("%s,%s\n" % (key, MOTIFSeq_dict[key]))
 
#########################
#     Step.3 合併檔案    #
#     以利下一步上傳到R   #
#########################           
 
import pandas as pd

ID_df = pd.read_csv(outpath + "/JASPAR_motif_ID.csv", header=None)
URL_df = pd.read_csv(outpath + "/JASPAR_motif_URL.csv", header=None)
ATCG_df = pd.read_csv(outpath + "/JASPAR_motif2030.csv", header=None)
ID_df.rename(columns={0: "TF", 1: "JASPAR ID"}, inplace=True)
print(ID_df)

URL_df.rename(columns={0: "TF", 1: "BEDURL"}, inplace=True)
print(URL_df)

ATCG_df.rename(columns={0: "JASPAR ID", 1: "MOTIF (20%30%)"}, inplace=True)
print(ATCG_df)

df1 = pd.DataFrame.merge(ID_df,ATCG_df,on='JASPAR ID')
df2 = pd.DataFrame.merge(df1,URL_df,on='TF')
df2['No'] = df2.index + 1
MOTIF_symbol = []
for item in df2["MOTIF (20%30%)"]:
    MOTIF_symbol.append(replace_values(item))
df2['MOTIF Symbol'] = MOTIF_symbol
df2['mers'] = df2['MOTIF Symbol'].apply(len)
df2 = df2[["No","TF","JASPAR ID","BEDURL","mers","MOTIF Symbol","MOTIF (20%30%)"]]
print(df2)

df2.to_csv(outpath + "/For_R_MOTIFs.csv", index = False)
