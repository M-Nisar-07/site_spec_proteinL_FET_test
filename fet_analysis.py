import pandas as pd
import sys
from utils import (get_query, get_d , generate_matrix, get_query_exp, get_query_pr)
from utils_profile import generate_matrix_pr
from itertools import chain

KINASE = 'PAK1'
SITE = 'S223'

            # ======================differential=====================

# df = pd.read_csv("input/pak1 data.csv")

# df["Exp_cond_id"] = df["Exp_cond_id"].apply(lambda x:str(x)+'_ph')

# df1,df2 = generate_matrix(df,KINASE,SITE)

# df1.to_excel(f"output/{KINASE}_{SITE}_UUDD_cross_talk.xlsx", index = False)

# df2.to_excel(f"output/{KINASE}_{SITE}_UDDU_cross_talk.xlsx", index = False)



            # ==========================PROFILING====================

# Please request the data management team for the raw data of each kianase
df = pd.read_csv("input/count_pak1_s223.csv")
df_exp = pd.read_csv("input/distinct_exp_condition_id_pak1.csv")

df["gene_site"] =  df.apply(lambda x: str(x["mapped_genesymbol"]) + '_' +str(x['mapped_phosphosite']) , axis=1 )

df["Exp_cond_id"] = df["all_exp_conditions"].apply(lambda x:x.split(';'))

sample = set(list(chain.from_iterable(df["Exp_cond_id"].tolist())))

print(len(sample))

df = generate_matrix_pr(df,df_exp,KINASE,SITE)

df.to_excel(f"output/{KINASE}_{SITE}_profile_cross_talk.xlsx", index = False)

# =================================================================








