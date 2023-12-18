import pandas as pd 
from scipy.stats import fisher_exact
import numpy as np 
import os 
import pymysql
# from dotenv import read_dotenv

# read_dotenv()
# DB_HOST = os.getenv('DB_HOST')
# DB_USER = os.getenv('DB_USER')
# DB_PASSWORD = os.getenv('DB_PASSWORD')
# DB_NAME = os.getenv('DB_NAME')

def get_query(k,s):
    q = f'''SELECT mapped_gene, Exp_cond_id , mapped_phosphosite FROM phospodb_nisar_profile_data WHERE Exp_cond_id IN(
        select DISTINCT Exp_cond_id FROM phospodb_nisar_profile_data where mapped_gene ='{k}' AND mapped_phosphosite='{s}') '''
    return q

def get_query_pr(k,s):

    q = f'''SELECT mapped_gene, Exp_cond_id , mapped_phosphosite FROM phospodb_nisar_profile_data WHERE Exp_cond_id IN(
        select DISTINCT Exp_cond_id FROM phospodb_nisar_profile_data where mapped_gene ='{k}' AND mapped_phosphosite='{s}') '''
    return q



def get_query_exp(k,s):
    q = f'''select DISTINCT Exp_cond_id FROM phospodb_nisar_differential_data where mapped_genesymbol ='{k}' '''
    return q
 

def get_d(query):
    # connection = pymysql.connect(
    #     host=DB_HOST,
    #     port=3306,
    #     user=DB_USER,
    #     password=DB_PASSWORD,
    #     database=DB_NAME,
    # )
    # cursor = connection.cursor()
    # cursor.execute(query)
    # result = cursor.fetchall()
    # df_result = pd.DataFrame(result, columns=[desc[0] for desc in cursor.description])
    # cursor.close()
    # connection.close()
    # return df_result
    pass

def get_combinations(k,s):
    li = [[k,i] for i in s]
    df = pd.DataFrame(li , columns=['site1', 'site2'])
    df = df.loc[df["site2"] != k]

    return df

def get_n00(s1,s2,t):
    s1 = {i.split('--&&--')[0] for i in s1}
    s2 = {i.split('--&&--')[0] for i in s2}
    
    return len(t - s1.union(s2))


def get_nud(s1,s2):

    s1 = [i.split('--&&--')[0] for i in s1]
    s2 = [i.split('--&&--')[0] for i in s2]

    s1_f = [i for i in s1 if i not in s2 ]
    s2_f = [i for i in s2 if i not in s1]

    for i in s2_f:
        s1_f.append(i)

    return len(set(s1_f))

def get_ndu(s1,s2):
    list_of_dicts = []

    for s1s in s1:
        c1,e1 = s1s.split('--&&--')
        for s2s in s2:
            c2,e2 = s2s.split('--&&--')
            if (c1 == c2) and (e1 == "Up-regulated") and (e2 == "down-regulated"):
                    list_of_dicts.append(c1)

    for s1s in s1:
        c1,e1 = s1s.split('--&&--')
        for s2s in s2:
            c2,e2 = s2s.split('--&&--')
            if (c1 == c2) and (e1 == "down-regulated") and (e2 == "Up-regulated"):
                list_of_dicts.append(c1)

    return len(set(list_of_dicts))


def get_uudd(s1,s2):
    total_matches = sum(pair in s2 for pair in s1)
    return total_matches

def get_n(s1,s2,df, total):

    s1_con = set(df.at[s1,'condition_exp'])
    s2_con = set(df.at[s2,'condition_exp'])

    n_00 = get_n00(s1_con,s2_con,total)
    n_ud = get_nud(s1_con,s2_con)
    n_du = get_ndu(s1_con,s2_con)
    n_uudd = get_uudd(s1_con,s2_con)

    data1 = np.array([[n_00,n_ud],[n_du,n_uudd]])
    _ ,p_value = fisher_exact(data1, alternative = 'greater')    

    data2 = np.array([[n_00,n_ud],[n_uudd,n_du]])
    _ ,p_value_r = fisher_exact(data2, alternative = 'greater')

    return n_00, n_ud, n_du, n_uudd, p_value , p_value_r


def generate_matrix(df,k,s):

    new_total_ex = set(df['Exp_cond_id'].tolist())
    
    print(len(new_total_ex))

    df['condition_exp'] = df[['Exp_cond_id', 'expression']].apply(lambda x: '--&&--'.join(map(str, x)), axis=1)

    df["gene_site"] = df["mapped_genesymbol"] + '_' +df['mapped_phosphosite']

    sub = df['gene_site'].unique().tolist()
    
    df_comb = get_combinations(k+"_"+s,sub)

    df = df.groupby('gene_site').agg(pd.Series.tolist).reset_index()

    df.set_index("gene_site", inplace = True) 

    df_comb[['n_00','n_uddu_nd','n_uddu','n_uudd','p-Value','p-Value_r']] = df_comb.apply(lambda x:get_n(x['site1'],x['site2'],
                                                                                                          df , new_total_ex),  axis = 1, result_type='expand')

    df_comb.drop_duplicates(inplace=True)
    
    df1 = df_comb[['site1','site2','n_00','n_uddu_nd','n_uddu','n_uudd','p-Value']]
    df2 = df_comb[['site1','site2','n_00','n_uddu_nd','n_uudd','n_uddu','p-Value_r']]
    
    df1.sort_values(by=['p-Value'], inplace=True)
    df2.rename({'p-Value_r':'p-Value'},axis=1, inplace = True)
    df2.sort_values(by=['p-Value'], inplace=True)

    return df1 , df2
