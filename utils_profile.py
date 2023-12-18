import pandas as pd 
from scipy.stats import fisher_exact
import numpy as np 
import os 
from itertools import chain


def get_combinations(k,s):
    li = [[k,i] for i in s]
    df = pd.DataFrame(li , columns=['site1', 'site2'])
    df = df.loc[df["site2"] != k]

    return df


def get_n00(s1,s2,t):
    return len(t - (s1.union(s2)))
    
def get_n01(s1,s2,t):
    return len(s2 - (s1.intersection(s2)))

def get_n10(s1,s2,t):
    return len(s1 - (s1.intersection(s2)))

def get_n11(s1,s2,t):
    return len(s1.intersection(s2))

def get_n(s1,s2,df, total):
    s1_con = df.at[s1,'Exp_cond_id']
    s2_con = df.at[s2,'Exp_cond_id']
    s1_con = set(s1_con)
    s2_con = set(s2_con)
    total = set(total)

    n_00 = get_n00(s1_con,s2_con,total) 
    n_01 = get_n01(s1_con,s2_con,total)
    n_10 = get_n10(s1_con,s2_con,total)
    n_11 = get_n11(s1_con,s2_con,total)
    data = np.array([[n_00,n_01],[n_10,n_11]])
    _ ,p_value = fisher_exact(data, alternative = 'greater')
    return n_00,n_01,n_10,n_11 , p_value


def generate_matrix_pr(df,df_exp,k,s):

    new_total_ex = set(df_exp['exp_cond_id'].tolist())
    print(len(new_total_ex))

    sub = df['gene_site'].unique().tolist()
    
    df_comb = get_combinations(k+"_"+s,sub)

    df.set_index("gene_site", inplace = True) 

    df_comb[['n_00','n_01','n_10','n_11','p-Value']] = df_comb.apply(lambda x:get_n(x['site1'],x['site2'],
                                                                                                          df , new_total_ex),  axis = 1, result_type='expand')

    df_comb.drop_duplicates(inplace=True)
    
    df_comb.sort_values(by=['p-Value'], inplace=True)
    return df_comb
