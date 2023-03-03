print("This is a model based on correlated KO genes with target gene")

import scanpy as sc
import numpy as np
import pandas as pd
import anndata
import scvi

import scipy.stats as stats
from collections import defaultdict
import os

from warnings import filterwarnings
filterwarnings('ignore')

sc.settings.verbosity = 3
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white')


col_name  = ['progenitor','effector','terminal exhausted','cycling','other'] 
tgenes = ['Aqr', 'Bach2', 'Bhlhe40', 'Ets1', 'Fosb', 'Mafk', 'Stat3']
g_pert66 = ['Arid4b','Arid5b', 'Atf2','Batf','Crem','Ctnnb1','Dkk3','Dvl1','Dvl2','Dvl3',
            'Eef2','Egr1','Elf1','Eomes','Ep300','Ezh2','Foxm1','Foxo1','Foxp1',
            #'Fzd1', 
            'Fzd3','Fzd6','Gsk3b','Hif1a','Hmgb1','Hmgb2','Id2','Id3','Ikzf3','Il12rb1',
            'Il12rb2','Irf2','Irf9','Klf2','Ldhb','Lef1','Litaf','Lrp1','Myb','Nr3c1',
            'Nr4a1','Nr4a2','Nr4a3','Oxnad1',
            #'P2rx7', 
            'Prdm1','Rad21','Rela','Rps6','Runx2',
            'Runx3','Satb1','Sox4','Sp100','Sp140','Stat4','Sub1','Tbx21','Tcf3','Tcf7',
            'Tox','Tox2','Tpt1','Yy1','Zeb2','Zfp292']

def get_state_distr(adata, g_name):
    """
    naive cell state distribution is calculated based on 'unperturbed' cell population when 'g_name' has zero expression value
    adata: anndata
    g_name: target gene name, ex) 'Tcf7'
    """
    # boolean index for adata slicing in 'unperturbed' cell population
    b2_idx = (adata[:, g_name][adata[:, g_name].obs.condition == "Unperturbed"].to_df() == 0.0)[g_name].values 
    len_zero_expr2 = len(adata[adata.obs.condition == "Unperturbed"][b2_idx, g_name].obs)
    return (adata[adata.obs.condition == "Unperturbed"][b2_idx, g_name].obs.groupby("state").size()/len_zero_expr2 )[[ 
                                                 'progenitor', 'effector', 'terminal exhausted', 'cycling', 'other' ]].values

                     
def get_KO_state_cell_cnt(adata, name_gKO):
    # KO case
    num_cells = len(adata[adata.obs.condition == name_gKO].obs)
    ser_tmp = adata[adata.obs.condition == name_gKO].obs.groupby("state").size()
    try: 
        return ser_tmp[[ 'progenitor', 'effector','terminal exhausted', 'cycling', 'other']], len(ser_tmp), num_cells
    except:
        return ser_tmp, len(ser_tmp), num_cells
        
def build_sol(df_reg, df_cor, g_list, tbl_cor_genes, alpha=0.7):
    res_final = []
    for gname in g_list:
        if len(tbl_cor_genes[gname]) == 0:
            res_final.append(df_reg.loc[[gname]]) # fallback when no-correlated genes detected
        else:
            res_final.append(df_cor.loc[[gname]]*alpha + df_reg.loc[[gname]]*(1-alpha))

    return pd.concat(res_final, axis=0).reset_index(names=['gene'])        

def main():
    # Read data
    adata = sc.read_h5ad('../../../data/sc_training.h5ad') # path should be actual data location of yours
    
    ##### linear decoder variational autoencoder model from scvi-tools
    scvi.model.LinearSCVI.setup_anndata(adata, layer="rawcounts") # setup the anndata for scvi-tools
    mdl_ldvae = scvi.model.LinearSCVI(adata, n_latent=5) # initialize LinearSCVI model
    
    # load the pre-trained model
    mdl_name = 'TcellKO_mdl_ldvae_all_cond_Z5_submit'
    mdl_ldvae = scvi.model.LinearSCVI.load(os.path.join('./', mdl_name), adata)    
    
    ##### uncomment below 3 lines if you want to repeat to train the model 
    #os.system(f'cp -r {os.path.join("./", mdl_name)} {mdl_name + "_old"}') # backup the pre-loaded model before training
    #mdl_ldvae.train(max_epochs=250, plan_kwargs={'lr':5e-3}, check_val_every_n_epoch=10) # train for 250 epochs, compute metrics every 10 epochs
    #mdl_ldvae.save(os.path.join('./', mdl_name, overwrite=True) # save the freshly trained model

    Z_hat = mdl_ldvae.get_latent_representation()
    for i, z in enumerate(Z_hat.T):
        adata.obs[f'Z_{i}'] = z
    
    ##### weights of linear decoder
    loadings = mdl_ldvae.get_loadings() # weights

    df_tgenes = []
    df_top10pos = []
    df_top10neg = []
    df_all = []
    df_all_KO = []
    for clmn_ in loadings:
        loading_ = loadings[[clmn_]].sort_values(by=clmn_, ascending=False)
        loading_ ['gRank'] = list(map(lambda x: x+1, list(range(len(loading_)))))
        df_top10pos.append(loading_.head(10))
        df_tgenes.append(loading_.loc[[i for i in tgenes]])
        df_top10neg.append(loading_.tail(10))
        df_all.append(pd.concat([loading_.loc[[i for i in tgenes]].reset_index(),
                                 loading_.head(10).reset_index(),
                                 loading_.tail(10).reset_index()], 
                                axis=0))
        df_all_KO.append(loading_.loc[[i for i in g_pert66]])
    
    df_tg_lv = pd.concat(df_all, axis=1).head(7)[[f'Z_{i}' for i in range(5)]]
    df_tg_lv.index = tgenes
    df_all_KO_lv = pd.concat(df_all_KO, axis=1)[[f'Z_{i}' for i in range(5)]]

    ##### pearson-r b/w target gene and 66 KO genes 
    tbl_tg_kog = defaultdict(list)
    for tg in df_tg_lv.T:
        tbl_tg_kog[tg]
        for kog in df_all_KO_lv.T: 
            pr_pv = stats.stats.pearsonr(df_tg_lv.loc[tg], df_all_KO_lv.loc[kog])
            if pr_pv[0] > 0.7 and pr_pv[1] < 0.05: 
                print(tg,kog,pr_pv)
                tbl_tg_kog[tg].append(kog)
                
    print(tbl_tg_kog) # table of correlated KO genes
                   
    #####
    g_val = tgenes[:3] #['Aqr', 'Bach2', 'Bhlhe40']
    g_test = tgenes[3:] #['Ets1', 'Fosb', 'Mafk', 'Stat3']
    ser_KO_cases = adata.obs.groupby('condition').size().sort_values(ascending=False)

    df_all_cell_cnt = []
    for name_gKO in ser_KO_cases.index:
        df_all_cell_cnt.append(get_KO_state_cell_cnt(adata, name_gKO)[0])  # 66 KO
        
    df_all_cell_cnt = pd.concat(df_all_cell_cnt, axis=1).T
    df_all_cell_cnt.index = ser_KO_cases.index

    df_res_cnt = pd.concat([ser_KO_cases, df_all_cell_cnt.fillna(0.0)], axis=1)
    df_res_cnt.columns = ['num_cells']+list(df_all_cell_cnt.columns)
    
    ##### cell state by counted cell number
    l_all = []

    for k, v in tbl_tg_kog.items():
        l_tmp = []
        for item in v:
            l_tmp.append(df_res_cnt.loc[item][1:].values)
        df_tmp = pd.DataFrame(l_tmp)
        l_all.append(df_tmp.sum()/df_tmp.sum().sum())
        
    pd.DataFrame(l_all) # distribution by sum of cell number of all correlated KO genes per a target gene
    
    pred_cor_g_cnt = pd.DataFrame(data=l_all, index = tbl_tg_kog.keys())
    pred_cor_g_cnt.columns = col_name             

    #####
    # for validation genes
    result_val_reg = []
    for g_name in g_val:
        result_val_reg.append(get_state_distr(adata, g_name))
    
    df_val_reg = pd.DataFrame(data=result_val_reg, index=g_val, columns = col_name)
    df_val = build_sol(df_val_reg, pred_cor_g_cnt, g_val, tbl_tg_kog, alpha=1)
    df_val.columns = ['gene','a_i', 'b_i', 'c_i', 'd_i', 'e_i']
    print(df_val)
    df_val[['gene','a_i', 'b_i', 'c_i', 'd_i', 'e_i']].to_csv('../solution/validation_output.csv', index=False)


    # for test genes
    result_test_reg = []
    for g_name in g_test:
        result_test_reg.append(get_state_distr(adata, g_name))
        
    df_test_reg = pd.DataFrame(data=result_test_reg, index=g_test, columns = col_name)
    df_test = build_sol(df_test_reg, pred_cor_g_cnt, g_test, tbl_tg_kog, alpha=1)
    df_test.columns = ['gene','a_i', 'b_i', 'c_i', 'd_i', 'e_i']
    print(df_test)
    df_test[['gene','a_i', 'b_i', 'c_i', 'd_i', 'e_i']].to_csv('../solution/test_output.csv', index=False)
    
    print("Done")
if __name__ == '__main__':     
    main()
