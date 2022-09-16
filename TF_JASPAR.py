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
        "version": 'latest',
        "data_type": "ChIP-seq",
        "release": 2022
    }
    result = client.action(schema, action, params=params)
    MOTIF = result['results'][0]['matrix_id']
    return MOTIF

#輸入JASPAR ID，輸出MOTIF
def get_matrix_ATCG(m_ID):
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
    GACT_array = np.append(np.array([list(ATCG_Matrix.values())[0],np.array(list(ATCG_Matrix.values())[1])]), np.array([list(ATCG_Matrix.values())[2],np.array(list(ATCG_Matrix.values())[3])]), axis=0)
    perc = GACT_array / GACT_array.sum(axis=0)
    #print(np.array([[f"{i:.2%}" for i in val] for val in perc]))
    seq1 = []
    for i in range(len(perc[0])):
        seq2 = []
        Gfreq, Afreq, Cfreq, Tfreq = 0, 0, 0, 0 #設定MOTIF中ATCG的百分比切點
        if perc[0][i] > Gfreq:
            seq2.append('G')
        if perc[1][i] > Afreq:
            seq2.append('A')
        if perc[2][i] > Cfreq:
            seq2.append('C')
        if perc[3][i] > Tfreq:
            seq2.append('T')
        seq1.append(seq2)
    return seq1

#以轉錄因子nfkb設定
get_matrix_id("nfkb")

##############################
#  Step.1 輸入TF 找JASPAR ID  #
##############################

####Step.1-1輸入檔案，製作TF串列

tf_list = []
path = '/Users/alyion/Google Drive/研究所/❖Su LAB❖/9.開會報告/20220913Meeting（DEGs）/254_TFs.txt'
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
with open('/Users/alyion/Google Drive/研究所/❖Su LAB❖/9.開會報告/20220913Meeting（DEGs）/JASPAR_motif_ID.csv', 'w') as f:
    for key in tfYesmotif_dict.keys():
            f.write("%s, %s\n" % (key, tfYesmotif_dict[key]))
            
#########################
#       黏貼bed.url      #
#########################

tfmotifURL_dict = {}
for j in tfYesmotif_dict.keys():
    motif_id = tfYesmotif_dict[j]
    motif_url = 'https://jaspar.genereg.net/download/data/2022/bed/' + motif_id + '.bed'
    tfmotifURL_dict.setdefault(j, motif_url)
print(tfmotifURL_dict)
    
with open('/Users/alyion/Google Drive/研究所/❖Su LAB❖/9.開會報告/20220913Meeting（DEGs）/JASPAR_motif_URL.csv', 'w') as f:
    for key in tfmotifURL_dict.keys():
            f.write("%s, %s\n" % (key, tfmotifURL_dict[key]))

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
with open('/Users/alyion/Google Drive/研究所/❖Su LAB❖/9.開會報告/20220913Meeting（DEGs）/JASPAR_motif2030.csv', 'w') as f:
    for key in MOTIFSeq_dict.keys():
            f.write("%s, %s\n" % (key, MOTIFSeq_dict[key]))
            
print(get_matrix_ATCG('MA1106.1'))
print(get_matrix_ATCG('MA0488.1'))
print(get_matrix_ATCG('MA0491.1'))
print(get_matrix_ATCG('MA0079.3'))
print(get_matrix_ATCG('MA0083.2'))



