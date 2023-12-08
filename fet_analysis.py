import pandas as pd
import sys
from utils import (get_query, get_d , generate_matrix, get_query_exp)

KINASE = 'PAK1'
SITE = 'S204'

# ==========================PROFILING=============================
# q = get_query(KINASE , SITE)
# df = get_d(q)

# q = get_query_exp(KINASE , SITE)
# dfe = get_d(q)
# =================================================================


df = pd.read_csv("pak1 data.csv")

df["Exp_cond_id"] = df["Exp_cond_id"].apply(lambda x:str(x)+'_ph')

df1,df2 = generate_matrix(df,KINASE,SITE)

df1.to_excel(f"{KINASE}_{SITE}_UUDD_cross_talk.xlsx", index = False)

df2.to_excel(f"{KINASE}_{SITE}_UDDU_cross_talk.xlsx", index = False)
