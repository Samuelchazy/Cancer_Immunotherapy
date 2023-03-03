import scanpy as sc

from warnings import filterwarnings
filterwarnings('ignore')

sc.settings.verbosity = 3
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white')

import anndata
import scvi
import numpy as np
import pandas as pd
pd.options.display.max_rows = 999


def get_state_distr_Q_prop(adata, g_name):
    """
    naive cell state distribution is calculated based on 'unperturbed' cell population when 'g_name' has zero expression value
    adata: anndata
    g_name: target gene name, ex) 'Tcf7'
    """
    df_tmp=pd.DataFrame([], columns=['progenitor', 'effector', 'terminal exhausted', 'cycling', 'other' ]) # empty df for handling cell states < 5

    # boolean index for adata slicing in 'unperturbed' cell population
    b2_idx = (adata[:, g_name][adata[:, g_name].obs.condition == "Unperturbed"].to_df() == 0.0)[g_name].values 
    len_zero_expr2 = len(adata[adata.obs.condition == "Unperturbed"][b2_idx, g_name].obs)
    return (pd.concat([df_tmp.T,
                       adata[adata.obs.condition == "Unperturbed"][b2_idx, g_name].obs.groupby("state").size()/len_zero_expr2], axis=1)
                     .fillna(0.0).T
                     .values[0])
  

def get_KO_state_distr2(adata, name_gKO):
    """ state vector of KO gene """
    num_cells = len(adata[adata.obs.condition == name_gKO].obs)
    ser_tmp = adata[adata.obs.condition == name_gKO].obs.groupby("state").size()/num_cells
    try: 
        return ser_tmp[[ 'progenitor', 'effector','terminal exhausted', 'cycling', 'other']], len(ser_tmp), num_cells
    except:
        return ser_tmp, len(ser_tmp), num_cells
  
  
def f_obj_part_a(row):
    """ part-a scoring function """
    obj = np.array([1,0,0,0,0])
    return (row*obj).sum()
    
    
def f_obj_part_b(row):
    """ part-b scoring function """
    svec_unpert = np.array([0.0675, 0.2097, -0.3134, 0.3921])
    return (row/svec_unpert).sum()
  
  
def handle_zero_div_scalar(nom, denom, scale_f): # assume both div_factor, pred_state_vec are scalars
    if denom == 0:
        if nom == 0:
            return 0
        else:
            print(nom, nom*scale_f) # DEBUGGING: to trace values in this case
            return nom*scale_f # scale but avoid too large value => scale_f = smallest non-zero score of part-a/-b object
    else:
        return nom/denom   


def f_obj_ch3_Q_checkpoint_v3(row): # assuming row = ['gene', 'a_i', 'b_i', 'c_i', 'd_i', 'e_i']
    """ proposed scoring function for immune checkpoint blockade therapy """
    div_factor = get_state_distr_Q_prop(adata, row[0])[0]
    #scale_f is arbitrarily chosen, could be hyper-para; here, set as the part-a non-zero smallest number = 0.011241 of 64KO genes
    scale_f = 1/0.011241
    return handle_zero_div_scalar(row[1], div_factor, scale_f)
    

def f_obj_ch3_Q_cart_v2(row): # assuming row = ['gene', 'a_i', 'b_i', 'c_i', 'd_i']
    """ proposed scoring function for car t-cell therapy """
    div_factor = get_state_distr_Q_prop(adata, row[0])[:4]*np.array([1,1,-1,1])
    #scale_f is arbitrarily chosen, could be hyper-para; here, set as the part-b non-zero smallest number = 0.371409 of 64 KO genes
    scale_f = 1/0.371409
    return sum([handle_zero_div_scalar(n_i, d_i, scale_f) for n_i, d_i in zip(row[1:5], div_factor)])


def main():
    
    global adata, ser_64KO_cases, df_64KO_4scoring
    # read data
    adata = sc.read_h5ad('../data/sc_training.h5ad')
    
    # 64 knock-out genes
    ser_all_cases = adata.obs.groupby('condition').size().sort_values(ascending=False)
    ser_64KO_cases = ser_all_cases.drop(['Unperturbed','Fzd1','P2rx7'])

    # state vectors of 64KO genes
    df_64KO_4scoring = []
    for name_gKO in ser_64KO_cases.index:
        df_64KO_4scoring.append(get_KO_state_distr2(adata, name_gKO)[0])
        
    df_64KO_4scoring = pd.concat(df_64KO_4scoring, axis=1).T.fillna(0.0)
    df_64KO_4scoring.index = ser_64KO_cases.index

    # 1) adding cycling constraint
    df_64KO_4scoring.columns = [ 'a_i', 'b_i', 'c_i', 'd_i', 'e_i' ]
    df_64KO_4scoring['cycling_constraint']=df_64KO_4scoring.d_i.apply(lambda x: 1 if x > 0.05 else 0) 
    df_64KO_4scoring = df_64KO_4scoring.reset_index(names=['gene'])

    # 2) apply part-a objective and add its ranking
    df_64KO_4scoring['score_part_a'] = df_64KO_4scoring[['a_i', 'b_i', 'c_i', 'd_i', 'e_i']].apply(f_obj_part_a, axis=1)
    df_64KO_4scoring.sort_values(by=['score_part_a'], ascending=[False], inplace=True)
    df_64KO_4scoring['part_a_rank'] = list(map(lambda x: x+1, list(range(len(df_64KO_4scoring)))))
    print("part-a done")
    
    # 3) apply part-b objective and add its ranking
    df_64KO_4scoring['score_part_b'] = df_64KO_4scoring[['a_i', 'b_i', 'c_i', 'd_i']].apply(f_obj_part_b, axis=1)
    df_64KO_4scoring.sort_values(by=['score_part_b'], ascending=[False], inplace=True)
    df_64KO_4scoring['part_b_rank'] = list(map(lambda x: x+1, list(range(len(df_64KO_4scoring)))))
    print("part-b done")
    
    # 4) apply challenge3 Q-checkpoint objective and add its ranking
    df_64KO_4scoring['score_Q_chkpnt'] = df_64KO_4scoring[['gene', 'a_i', 'b_i', 'c_i', 'd_i', 'e_i']].apply(f_obj_ch3_Q_checkpoint_v3, axis=1)
    df_64KO_4scoring.sort_values(by=['score_Q_chkpnt'], ascending=[False], inplace=True)
    df_64KO_4scoring['Q_chkpnt_rank'] = list(map(lambda x: x+1, list(range(len(df_64KO_4scoring)))))
    print("Q-checkpoint done")
    
    # 5) apply challenge3 Q-cart objective and add its ranking
    df_64KO_4scoring['score_Q_cart'] = df_64KO_4scoring[['gene', 'a_i', 'b_i', 'c_i', 'd_i']].apply(f_obj_ch3_Q_cart_v2, axis=1)
    df_64KO_4scoring.sort_values(by=['score_Q_cart'], ascending=[False], inplace=True)
    df_64KO_4scoring['Q_cart_rank'] = list(map(lambda x: x+1, list(range(len(df_64KO_4scoring)))))
    print("Q-cart done")
    
    return df_64KO_4scoring
    
if __name__ == '__main__':
    df_64KO_4scoring = main()    
    print(df_64KO_4scoring[['gene', 'cycling_constraint', 
                            'score_part_a', 'part_a_rank', 
                            'score_part_b', 'part_b_rank',
                            'score_Q_chkpnt', 'Q_chkpnt_rank', 
                            'score_Q_cart', 'Q_cart_rank']].reset_index(drop=True))
