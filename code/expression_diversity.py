import os, sys
import numpy as np
import pandas as pd
import networkx as nx
from scipy.stats import mannwhitneyu
from sklearn.svm import SVC
from sklearn.metrics import accuracy_score
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import SelectFromModel
from sklearn.cluster import AgglomerativeClustering
import multiprocessing as mp



#--------------------------------------------------------------------
def load_ubiquitin():
    """
    Function to load auxiliary data on Ub ligases. 
    Combines data and annotations from two databases - hUbiquitome and UbiNet2.0 -
    and filters for genes in the nonredundant human Uniprot reference proteome
    Paths of input files hardcoded (below)
    """

    ub = pd.read_csv("../data/auxilliary/ubiquitin/hUbiquitome.csv", header=0, index_col=False, sep=',')
    ub1 = pd.read_csv("../data/auxilliary/ubiquitin/Categorized_human_E3_ligases.csv", header=0, index_col=False, sep=',')
    match_genes = pd.read_csv("../data/processed/id_human.txt", header=0, index_col=None, sep='\t')

    output = pd.DataFrame(columns=['Gene', 'Class', 'Category'])

    # E2: only one DB
    list_e2 = list(set(ub['E2']))

    # E3: merge 3 DBs and convert IDs
    list_e3 = []
    for i in list(ub['E3']): # + list_dub:
        if i in list(match_genes['uniprot']):
            current_gene = list(match_genes[match_genes['uniprot']==i]['gene'])[0]
            list_e3.append(current_gene)
    for i in list(set(ub1['E3'])):
        if i in list(match_genes['gene']):
            list_e3.append(i)
    list_e3 = list(set(list_e3))
    
    # DUB: only one DB, convert IDs
    list_dub = []
    for i in list(set(ub['DUB'])):
        if i in list(match_genes['uniprot']):
            current_gene = list(match_genes[match_genes['uniprot']==i]['gene'])[0]
            list_dub.append(current_gene)

    for i in list_e2:
        if str(i) != str('nan'):
            output.loc[len(output)]=(i, 'E2', "E2")

    for i in list_e3: # + list_dub:
        if i in list(ub1['E3']):
            current_domain = list( ub1[ub1['E3']==i]['Category'])[0]
            current_domain = current_domain.split('-')
            if current_domain[0] == "HECT" or len(current_domain) == 1:
                current_domain = current_domain[0]
            elif len(current_domain) > 1:
                current_domain = "_".join(current_domain[0:2])
        else:
            current_domain = "Uncategorized"
        output.loc[len(output)]=( i, "E3", current_domain)

    for i in list_dub:
        output.loc[len(output)]=( i, "DUB" , "DUB")

    return output



#--------------------------------------------------------------------
def parse_meta(meta, human_meta, mca_meta):
    """
    Function to extract info on clusters from the metadata files for coloring
    meta_types_clusters: for each cluster one DF entry
    meta_types_cells: for each cell one DF entry
    Output written to text files for plotting etc. 
    """

    meta_types = pd.DataFrame(columns=['cluster', 'type'])
    for i in list(set(list(meta['cluster']))):
        current_samples = list(meta[meta['cluster']==i]['sample'])
        current_cat = []
        cluster_type = 'None'

        for j in current_samples:
            if j in list(human_meta['sample_name']):
                current_type = str(list(human_meta[human_meta['sample_name']==j]['old.cluster'])[0])
                current_type = current_type.split()[-2]
            elif j in list(mca_meta['exp_component_name']):
                current_type = str(list(mca_meta[mca_meta['exp_component_name']==j]['cell_type_alias_label'])[0])
                current_type = current_type.split()[0]
            current_cat.append(current_type)

        current_cat2 = list(set(current_cat))


        if len(current_cat2) == 1:
            cluster_type = current_cat[0]

        elif len(current_cat2) > 1: 		# choice by qualified majority
            res = np.zeros(( len(current_cat2) ), dtype=int)
            for kx, k in enumerate(current_cat2):
                res[kx] = np.sum(np.array(current_cat)==k)
            best = np.argmax(res)
            cluster_type = current_cat2[best]

        meta_types.loc[len(meta_types)] = (i, cluster_type)

    meta_types.to_csv("../data/processed/meta_types_clusters.txt", header=True, index=False, sep='\t')
  
    meta_types_cell = pd.DataFrame(columns=['sample', 'cluster', 'type', 'subtype'])
    for i in list(meta['sample']):

        current_cluster = str( list(meta[meta['sample']==i]['cluster'])[0] )
        if len(current_cluster) == 0:
            current_cluster = "NA"

        current_type = list(meta_types[meta_types['cluster']==str(current_cluster)]['type'])
        if len(current_type) > 0:
            current_type = current_type[0]
        else:
            current_type = "NA"

        meta_types_cell.loc[len(meta_types_cell)] = (i, current_cluster, current_type, None)

    inh = 1
    exc = 1
    nnn = 1
    for i in list(set(meta_types_cell['cluster'])):
        current_type = str( list(meta_types_cell[meta_types_cell['cluster']==i]['type'])[0] )
        if len(current_type) > 0 and current_type == "Inh":
            current_subtype = current_type + str(inh)
            inh += 1
        elif len(current_type) > 0 and current_type == "Exc":
            current_subtype = current_type + str(exc)
            exc += 1
        else:
            current_subtype = "non" + str(nnn)
            nnn += 1

        current_sel = meta_types_cell['cluster'] == i
        meta_types_cell['subtype'][current_sel] = current_subtype

    meta_types_cell.to_csv("../data/processed/meta_types_cell.txt", header=True, index=False, sep='\t')

    return meta_types



#--------------------------------------------------------------------
def filterbylists(expMAT, varMAT, li_channels, li_synapse, li_chap, li_ub):
    """
    Load auxiliary gene lists of interest, match with genes in expression data
    Data sources (see paper for references):
        Ion channels: ModelDB
        Synapse proteome: SynaptomeDB minus channels, chaperones, and ubiquitin genes
        Chaperones: Proteostasis Consortium
        Ub genes: hUbiquitome and UbiNet2.0
    """

    def get_subset(inDF, inList):
        list_new = []
        for i in inList:
            if i in list(inDF.index):
                list_new.append(i)
        list_new = sorted(list_new)
        newDF = inDF.loc[list_new]
        return newDF

    # make synapse exclusive of other stuff
    list_synapse2 = []
    for i in li_synapse:
        if i not in li_channels and i not in li_chap and i not in li_ub:
            list_synapse2.append(i)

    expression_channels = get_subset(expMAT, li_channels) 
    expression_channels.to_csv("../data/processed/expression_channels.txt", header=True, index=True, sep='\t')
    variability_channels = get_subset(varMAT, li_channels)
    variability_channels.to_csv("../data/processed/variability_channels.txt", header=True, index=True, sep='\t')

    expression_synapse = get_subset(expMAT, list_synapse2)
    expression_synapse.to_csv("../data/processed/expression_synapse.txt", header=True, index=True, sep='\t')
    variability_synapse = get_subset(varMAT, list_synapse2)
    variability_synapse.to_csv("../data/processed/variability_synapse.txt", header=True, index=True, sep='\t')

    expression_chap = get_subset(expMAT, li_chap)
    expression_chap.to_csv("../data/processed/expression_chap.txt", header=True, index=True, sep='\t')
    variability_chap = get_subset(varMAT, li_chap)
    variability_chap.to_csv("../data/processed/variability_chap.txt", header=True, index=True, sep='\t')

    expression_ub = get_subset(expMAT, li_ub)
    expression_ub.to_csv("../data/processed/expression_ub.txt", header=True, index=True, sep='\t')
    variability_ub = get_subset(varMAT, li_ub)
    variability_ub.to_csv("../data/processed/variability_ub.txt", header=True, index=True, sep='\t')


    return expression_channels, expression_synapse, expression_chap, expression_ub



#--------------------------------------------------------------------
def diffexp(dataIN, dataIndex, subset, meta1, meta2, colID1, colID2):
    """
    Differential expression (DE) testing by nonparametric MWU test as in Seurat
    Add colIDs: for array indexing of conditions
    """

    list_subset = []
    for i in list(subset.index):
        if i in dataIndex:
            list_subset.append(dataIndex.index(i))
    data = dataIN[list_subset,:]
    data = np.array(data)

    types = []
    for ix, i in enumerate(list(meta1['sample'])):
        current_cluster = list(meta1[meta1['sample']==i]['cluster'])[0]
        current_type = list(meta2[meta2['cluster']==current_cluster]['type'])[0]
        if current_type == "Exc":
            types.append("Exc")
        elif current_type == "Inh":
            types.append("Inh")
        else:
            types.append("non")
    types = np.array(types, dtype=str)

    subset1 = types == colID1
    subset2 = types == colID2

    N = np.shape(data)[0]
    result = np.zeros(( N  ))
    res_fc = np.zeros(( N ))
   
    for i in range(N):
        data_i = data[i, subset1]
        data_j = data[i, subset2]
        if (np.sum(data_i) + np.sum(data_j) > 0) and np.any(data_i != data_j):
            pval = mannwhitneyu(data_i, data_j, alternative='two-sided')[1]
        else:
            pval = np.nan

        result[i] = pval        
        res_fc[i] = np.mean(data_i + 1) / np.mean(data_j + 1)   # pseudo count of 1, like Seurat

        #if np.sum(data_j) > 0:
        #    fc = np.mean(data_i[data_i>0]) / np.mean(data_j[data_j>0])
        #    #fc = np.around(np.log2( np.mean(data_i[data_i>0]) / np.mean(data_j[data_j>0]) ), 4)
        #else:
        #    fc = np.nan
        #res_fc[i] = fc
   
    return result, res_fc



#--------------------------------------------------------------------
def pval_adjust(pval):
    """
    correct for multiple testing
    Benjamini & Hochberg FDR method (R funct p.adjust)
    """
    pval = np.array(pval)
    if np.size(pval) > 1:
        padj = np.nan * np.ones( (np.size(pval) ) )
        nn = np.isnan(pval)
        pval = pval[~nn]
        n = len(pval)
        i = np.array( range(n)[::-1]) + 1
        o = np.array( sorted(range(n), key=lambda k: pval[k])[::-1] )
        ro = np.array( sorted(o, key=lambda k: o[k]) )
        adj = np.minimum.accumulate( float(n)/np.array(i) * pval[o] )
        for i in range(len(adj)):
            adj[i] = min(1, adj[i])
        padj[~nn] = adj[ro]
    else:
        padj = pval

    return padj



#--------------------------------------------------------------------
def diffexp_pairwise(dataIN, indexIN, list_exp, metaIN, metaClusterIN, classlist):
    """
    Pairwise DE testing with FDR ajustment
    """

    queries = [ ['Exc', 'Inh'], ['Exc', 'non'], ['Inh', 'non'] ]

    resultDF = pd.DataFrame(columns=['class', 'label', 'dir', 'count'])
    for i in range(len(list_exp)):
        current_exp = list_exp[i]
        current_class = classlist[i]

        for j in range(len(queries)):
            current_query = queries[j]
            current_q1 = current_query[0]
            current_q2 = current_query[1]
            current_label = "/".join([current_q1, current_q2])

            pval, fc = diffexp(dataIN, indexIN, current_exp, metaIN, metaClusterIN, current_q1, current_q2)
            adj = pval_adjust(pval)
            current_up = np.nansum( fc > 1.2 ) / len(current_exp) * 100
            current_down = np.nansum( fc < 0.8 ) / len(current_exp) * 100

            resultDF.loc[len(resultDF)] = (current_class, current_label, 'up', np.around(current_up,3) )
            resultDF.loc[len(resultDF)] = (current_class, current_label, 'down', np.around(current_down,3) )

    print(resultDF)
    resultDF.to_csv("../data/processed/classes_diffexp.txt", header=True, index=False, sep='\t')



#--------------------------------------------------------------------
def parse_heatmap(avrgEXP, EXP, metaDF, suffixOUT):
    """
    wrangles data into heatmap
    x: cosine similarity, clustered by full expression profile within Exc, Inh, non
    y: cosine similarity, clustered by genes
    values: log2FC relative to overall median expression
    """ 

    # rel to median
    rel_data = np.array(avrgEXP)
    for i in range(np.shape(rel_data)[0]):
        current_median = np.nanmedian(rel_data[i,:])
        if current_median > 0:
            rel_data[i, :] /= current_median         # FC rel to median
    resultDF = pd.DataFrame(data=rel_data, columns=avrgEXP.columns, index=avrgEXP.index)


    avrg_data = np.array(avrgEXP) 
    list_avrg = list(avrgEXP.columns)

    sel_exc = np.zeros(( len(list_avrg) ), dtype=bool)
    sel_inh = np.zeros(( len(list_avrg) ), dtype=bool)
    sel_non = np.zeros(( len(list_avrg) ), dtype=bool)

    for ix, i in enumerate(list_avrg):
        current_type = list(metaDF[metaDF['cluster']==i]['type'])[0]
        if current_type == "Exc":
            sel_exc[ix] = True
        elif current_type == "Inh":
            sel_inh[ix] = True
        else:
            sel_non[ix] = True

    list_columns = []
    rel_data = np.array(resultDF)
    # cluster columns based on full expression data
    for ix, current_sel in enumerate([sel_exc, sel_inh, sel_non]):
        dat = rel_data[:,current_sel]
        dat = np.transpose(dat)
        model = AgglomerativeClustering(n_clusters=np.shape(dat)[0], affinity='cosine', linkage='average', compute_full_tree=True)
        model.fit(dat)
        labels = model.labels_
        current_list = np.array(list_avrg)[current_sel]
        list_columns += list(current_list[labels])

    rel_dataDF = resultDF.loc[list(EXP.index)]
    data = np.array(rel_dataDF)
    # cluster input genes
    model = AgglomerativeClustering(n_clusters=np.shape(data)[0], affinity='cosine', linkage='average', compute_full_tree=True)
    model.fit(data)
    gene_labels = model.labels_
    rel_dataDF = rel_dataDF.iloc[gene_labels][list_columns]

    # add some empty columns for plotting
    spacer = 10
    pos1 = int( np.sum(sel_exc) )
    pos2 = int( pos1 + spacer + np.sum(sel_inh) )

    for i in range(spacer):
        rel_dataDF.insert( pos1+i, "space1"+str(i), list(np.repeat("NA", np.shape(data)[0] )) )
    for i in range(spacer):
        rel_dataDF.insert( pos2+i, "space2"+str(i), list(np.repeat("NA", np.shape(data)[0] )) )
   
    rel_dataDF.to_csv("../data/processed/heatmap_"+str(suffixOUT)+".txt", header=True, index=True, sep='\t')

    return rel_dataDF


#--------------------------------------------------------------------
def summary_heatmap(HMAP, metaDF, suffixOUT):
    """
    bla
    """

    list_columns = list(HMAP.columns)

    sel_exc = np.zeros(( len(list_columns) ), dtype=bool)
    sel_inh = np.zeros(( len(list_columns) ), dtype=bool)
    sel_non = np.zeros(( len(list_columns) ), dtype=bool)

    for ix, i in enumerate(list_columns):
        if i in list(metaDF['cluster']):
            current_type = list(metaDF[metaDF['cluster']==i]['type'])[0]
            if current_type == "Exc":
                sel_exc[ix] = True
            elif current_type == "Inh":
                sel_inh[ix] = True
            else:
                sel_non[ix] = True

    data_exc = np.array(HMAP)[:,sel_exc]
    data_inh = np.array(HMAP)[:,sel_inh]
    data_non = np.array(HMAP)[:,sel_non]

 

    exc_up = np.around(np.nansum(data_exc > 2, 1)/ np.shape(data_exc)[1] * 100, 2)       # pct without nan vals
    exc_do = np.around(np.nansum(data_exc < 0.5, 1)/np.shape(data_exc)[1] * 100, 2)
    inh_up = np.around(np.nansum(data_inh > 2, 1)/np.shape(data_inh)[1] * 100, 2)
    inh_do = np.around(np.nansum(data_inh < 0.5, 1)/np.shape(data_inh)[1] * 100, 2)
    non_up = np.around(np.nansum(data_non > 2, 1)/np.shape(data_non)[1] * 100, 2)
    non_do = np.around(np.nansum(data_non < 0.5, 1)/np.shape(data_non)[1] * 100, 2)


    res = pd.DataFrame({'exc_up': list(exc_up), 'exc_down': list(exc_do), 
        'inh_up': list(inh_up), 'inh_down': list(inh_do), 
        'non_up': list(non_up), 'non_do': list(non_do)}, index=list(HMAP.index))
    print(res)
    
    res.to_csv("../data/processed/log2fc_sig_"+str(suffixOUT)+".txt", header=True, index=True, sep='\t')


#--------------------------------------------------------------------
def total_abundance_category(expIN, annoIN, selColumn, suffixOUT):
    """
    Computes cumulative expression of genes in input categories
        expIN: expression dasta
        annoIN: annotation of protein homeostasis genes
        selColumn: workaround for different column names in different annoIN dataframes ... 
    """

    list_index = []
    result_data = np.zeros(( len(list(set(annoIN[selColumn]))), np.shape(np.array(expIN))[1]  ))
    for ix, i in enumerate( list(set(annoIN[selColumn])) ): 
        current_genes = list(annoIN[annoIN[selColumn]==i]['Gene'])
        list_genes = []
        for j in current_genes:
            if j in list(expIN.index):
                list_genes.append(j)
        list_index.append(i)
        result_data[ix,:] = np.sum( np.array(expIN.loc[list_genes]), 0) #/ len(list_genes)

    resultDF = pd.DataFrame(data=result_data, index=list_index, columns=expIN.columns)
    resultDF.to_csv("../data/processed/abundance_bycategories_"+str(suffixOUT)+".txt", header=True, index=True, sep='\t')

    print(resultDF)
    






if __name__ == '__main__':


    ## LOAD META DATA
    meta_human = pd.read_csv("../data/AllenBrainAtlas/human_LGN_2021_metadata.csv", header=0, index_col=False, sep=',')
    meta_mca = pd.read_csv("../data/AllenBrainAtlas/MCA_metadata.csv", header=0, index_col=False, sep=',')
    meta_merged = pd.read_csv("../data/processed/merged_meta.txt", header=0, index_col=False, sep='\t')
    data_merged = np.loadtxt("../data/processed/data_merged.txt")

    ## LOAD AUX DATA
    channels = pd.read_csv("../data/processed/list_channels.txt", header=None, index_col=False)
    syn = pd.read_csv("../data/auxilliary/synapse/synapse_proteome.txt", header=0, index_col=False, sep='\t')
    chap = pd.read_csv("../data/auxilliary/list_chaperones.txt", header=0, index_col=False, sep=' ')

    if os.path.exists("../data/auxilliary/ub.csv"):
        ub = pd.read_csv("../data/auxilliary/ub.csv", header=0, index_col=False, sep='\t')
    else:
        ub = load_ubiquitin()
        ub.to_csv("../data/auxilliary/ub.csv", header=True, index=False, sep='\t')

    list_channels = list(set(channels[0]))
    list_chap = list(set(chap['Gene']))
    list_ub = list(set(ub['Gene']))
    list_synapse = []						# make synapse stuff all other than chap and ub
    for i in list(set(syn['symbol'])):
        if i not in list_channels and i not in list_chap and i not in list_ub:
            list_synapse.append(i)

    ## LOAD EXP DATA
    avrg_exp = pd.read_csv("../data/processed/cluster_average_norm.txt", header=0, index_col=0, sep='\t')
    variab = pd.read_csv("../data/processed/cluster_variable.txt", header=0, index_col=0, sep='\t')


    print("computing now")
    ## COMPUTE FUNCTIONS
    meta_clusters = parse_meta(meta_merged, meta_human, meta_mca)
    exp_channels, exp_synapse, exp_chap, exp_ub = filterbylists(avrg_exp, variab, list_channels, list_synapse, list_chap, list_ub)

    diffexp_pairwise(data_merged, list(avrg_exp.index), [exp_channels, exp_synapse, exp_chap, exp_ub], meta_merged, meta_clusters, ['Ch', 'Syn', 'Chap', 'Ub'])

    hmap_chap = parse_heatmap(avrg_exp, exp_chap, meta_clusters, "chap")
    hmap_ub = parse_heatmap(avrg_exp, exp_ub, meta_clusters, "ub")
    
    summary_heatmap(hmap_chap, meta_clusters, "chap")
    summary_heatmap(hmap_ub, meta_clusters, "ub")

    total_abundance_category(exp_chap, chap, "System", "chap")  
    total_abundance_category(exp_ub, ub, "Category", "ub")