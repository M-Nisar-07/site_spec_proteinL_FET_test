import pandas as pd
import sys
from utils import (get_query, get_d , generate_matrix, get_query_exp)

KINASE = 'PAK1'
SITE = 'S144'

# ==========================PROFILING=============================
# q = get_query(KINASE , SITE)
# df = get_d(q)

q = get_query_exp(KINASE , SITE)
dfe = get_d(q)

df = pd.read_excel("pak1_diff_data.xlsx", engine="openpyxl")

df = generate_matrix(df,KINASE,SITE, dfe)

df.to_excel(f"{KINASE}_{SITE}_cross_talk.xlsx")

