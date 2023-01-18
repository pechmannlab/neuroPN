import os, sys
import numpy as np
import pandas as pd
import networkx as nx
from Bio import SeqIO
from scipy.stats import mannwhitneyu



#--------------------------------------------------------------------
def extract_cluster(fileIN, metaIN):
    """
    Extract the mapping between experiment id and cluster label
    """

    counter = 0
    for line in open(fileIN, 'r'):
        current_line = line.strip('\n').replace('"','').split(',')
        current_name = current_line[0]

        if counter == 0:                    # header row
            current_sel = list(metaIN['exp_component_name'])

            list_cells = []
            list_names = []
            list_genes = []

            for i in current_sel:
                if i in current_line:
                    list_cells.append(current_line.index(i))
                    list_names.append(i)
            list_cells = np.array(list_cells, dtype=int)

            result = np.zeros(( len(match_human), len(list_cells) ), dtype=int)
            counter += 1

        elif counter > 0:                   # data rows
            if current_name in list(match_human['gene']):
                list_genes.append(current_line[0])
                current_line = np.array(current_line)[list_cells]
                result[counter-1,:] = current_line          # header row not needed
                counter += 1
   
    list_cluster = []
    list_region = []
    for i in list_names:
        current_meta = metaIN[metaIN['exp_component_name']==i]
        current_cluster = str(list(current_meta['cluster_label'])[0]).replace(" ", "_")
        if str(current_cluster) == "nan":
            current_cluster = 'NA'
        current_cluster = current_cluster.replace('-', '_') # problem with '-' in colname in R
        current_region = list(current_meta['region_label'])[0]

        list_cluster.append(current_cluster)
        list_region.append(current_region)

    metaDF = pd.DataFrame({'sample': list_names, 'cluster': list_cluster, 'region': list_region})


    return metaDF, list_genes, result



#--------------------------------------------------------------------
def normalize_countmat(dataIN):
    """
    Normalize to CPM (counts per million)
    """

    dataOUT = np.zeros(( np.shape(dataIN) ))


    for i in range(np.shape(dataIN)[1]):
        dataOUT[:,i] = np.around( dataIN[:,i]/np.sum(dataIN[:,i]) * 1000000, 4)

    
    return dataOUT



#--------------------------------------------------------------------
def average_expression(dataIN):
    """
    Compute average expression of cluster
    follows Seurat implementation: nonzero entries, no pseudo-count
    """
    
    dataIN = np.array(dataIN) 

    dat = np.sum(np.transpose(dataIN), 0)
    nz  = np.sum(np.transpose(dataIN) > 0, 0)

    if np.all( nz > 0):
        avrg = dat/nz
    else:
        avrg = np.zeros(( len(dat) ))
        avrg[nz>0] = dat[nz>0] / nz[nz>0]

    avrg_norm = np.around( avrg/np.nansum(avrg) * 1000000, 4)

    return avrg, avrg_norm



#--------------------------------------------------------------------
def diffexp(dataIN, clustIDs):
    """
    DE (differential expression) testing by nonparametric MWU test
    FC computed with pseudo-count of 1; as implemented in Seurat; 
    """

    list_clusters = sorted(list(set(clustIDs)))
    L = len(list_clusters)
    N = np.shape(dataIN)[0]
    result = np.zeros(( N, int((L*(L-1)/2.)) ))
    res_fc = np.zeros(( N, int((L*(L-1)/2.)) ))


    idx = 0
    for i in range(L):
        print(i, list_clusters[i])
        sel_i = np.array(clustIDs) == list_clusters[i]
        for j in range(i+1,L):
            sel_j = np.array(clustIDs) == list_clusters[j]
 
            for k in range(N):
                data_i = dataIN[k, sel_i]
                data_j = dataIN[k, sel_j]

                if (np.nansum(data_i) > 0) and (np.nansum(data_j) > 0) and np.any(data_i != data_j):
                    pval = mannwhitneyu(data_i, data_j, alternative='two-sided')[1]
                    #fc = np.mean(data_i[data_i>0]) / np.mean(data_j[data_j>0])
                    fc = np.nanmean(data_i + 1) / np.nanmean(data_j + 1)
                else:
                    pval = np.nan
                    fc = np.nan
                                  
                result[k, idx] = pval
                res_fc[k,idx] = fc

            idx += 1

    return result, res_fc



#--------------------------------------------------------------------
def match_genenames(seqlist, uniDF, genesDF):
    """
    Best match of gene names between annotations from different sources
    """

    resultDF = pd.DataFrame(columns=['uniprot', 'unigene', 'gene', 'index', 'order'])

    list_matched = []

    for i in range(len(uniDF)):
        current_data = uniDF.iloc[i]
        current_entry = current_data['Entry name']
        if current_entry in seqlist:
            current_names = list( str(current_data['Gene names']).split() )

            for j in current_names:
                current_gene = list(genesDF[genesDF['gene']==j]['gene'])
                current_idx = list(genesDF[genesDF['gene']==j].index)
                current_order = current_names.index(j)
                if len(current_gene) > 0 and current_gene[0] not in list_matched:
                    resultDF.loc[len(resultDF)] = (current_entry, j, current_gene[0], current_idx[0], current_order)
                    list_matched.append(current_gene[0]) # only one entry per gene
                    break

    print(len(uniDF), len(genesDF), len(resultDF), len(set(resultDF['gene'])))

    return resultDF

 

#--------------------------------------------------------------------
def variable_genes(dataIN, genesIN):
    """
    Detect variable genes within clusters
    Testing 20 bins of equal size; following the Seurat implementation
    """

    list_names = list(genesIN)
    data = np.array(dataIN)

    result = np.zeros(( len(list_names) ), dtype=int)

    sel = np.sum(data,1) > 0
    data = data[sel,:]
    names_sel = np.array(list_names)[sel]

    x = np.sum(data, axis=1) / np.sum(data>0, axis=1)
    y = np.var(data, axis=1) / x

    x_min = np.nanmin(x)
    x_max = np.nanmax(x)
    bin_width = (x_max - x_min) / 20.

    for i in range(20):
        current_start = x_min + i * bin_width
        current_end = x_min + (i+1) * bin_width

        current_sel = (x >= current_start) * (x < current_end)
        current_x = x[current_sel]
        current_y = y[current_sel]
        current_names = names_sel[current_sel]

        current_y2 = current_y[np.isnan(current_y)==False]
        current_names2 = current_names[np.isnan(current_y)==False]
        if len(current_y2) > 0:
            current_ycut = np.mean(current_y2) + 2* np.std(current_y2)

            N = np.sum(current_y2 > current_ycut)
            R = current_names2[current_y2 > current_ycut]
            for j in list(R):
                current_idx = list_names.index(j)
                result[current_idx] = 1

    return result




#--------------------------------------------------------------------
def cocorrelation(coexpMAT, theta, list1, list2, preOUT):
    """
    Extracts cocorrelation network from coexpression matrix
    Filtered based on threshold theta and input lists list1 and list2
    """

    N,M = np.shape(coexpMAT)

    bins = np.arange(1, 41)/20 -1
    dat = coexpMAT[np.isnan(coexpMAT)==False]
    hist, be = np.histogram(dat, bins=bins, density=False)
    np.savetxt("../data/processed/histogram_"+preOUT+".txt", hist)

    G = nx.Graph()
    for i in range(N):
        current_query = list1[i]
        for j in range(i+1, M):
            current_test = list2[j]
            current_score = coexpMAT[i,j]     
            if current_score >= theta:
                #G.add_edge(i,j, weight=current_score)
                G.add_edge(current_query, current_test, weight=current_score)
             
    nx.write_gpickle(G, "../data/processed/coexp_"+preOUT+".nx")
    


#--------------------------------------------------------------------
def cluster_expression(metaDF, dataMAT, list_genes):
    """
    Function to export:
    - average expression per cluster
    - variable genes per cluster
    Removes all-zero genes that are never expressed
    """

    data_avrg = np.zeros(( np.shape(dataMAT)[0], len(set(metaDF['cluster'])) ))
    data_avrg_norm = np.zeros(( np.shape(dataMAT)[0], len(set(metaDF['cluster'])) ))

    list_avrg = list(set(metaDF['cluster']))
    data_var = np.zeros(( np.shape(dataMAT)[0], len(set(metaDF['cluster'])) ))

    for ix, i in enumerate(list_avrg):
        if i != "nan":
            current_idx = metaDF['cluster'] == str(i)
            current_data = dataMAT[:,current_idx] 
            current_var = variable_genes(current_data, list_genes)
            data_var[:,ix] = current_var 
            current_avrg, current_avrg_norm = average_expression(current_data)
            data_avrg[:,ix] = current_avrg
            data_avrg_norm[:,ix] = current_avrg_norm

    df_avrg = pd.DataFrame(data=data_avrg, columns=list_avrg, index=list_genes)
    df_avrg.to_csv("../data/processed/cluster_average.txt", header=True, index=True, sep='\t')
    
    df_avrg_norm = pd.DataFrame(data=data_avrg_norm, columns=list_avrg, index=list_genes)
    df_avrg_norm.to_csv("../data/processed/cluster_average_norm.txt", header=True, index=True, sep='\t')
    
    df_var = pd.DataFrame(data=data_var, columns=list_avrg, index=list_genes)
    df_var.to_csv("../data/processed/cluster_variable.txt", header=True, index=True, sep='\t')

 

#--------------------------------------------------------------------
def average_foldchange(fcMAT):
    """
    Computes average foldchange for each gene to identify variable genes between clusters
    """

    fc_mean = np.zeros(( np.shape(fcMAT)[1] ))
    for i in range(np.shape(fcMAT)[1]):
        current_data = fcMAT[i,:]
        if np.any(np.isnan(current_data)==False):
            current_data = current_data[np.isnan(current_data)==False]
            current_data = np.log2(current_data)        # necess bc setting stuff to zero above (check!)
            fc_mean[i] = np.mean(np.abs(current_data))

    np.savetxt("../data/processed/average_foldchange.txt", fc_mean, fmt='%.4f')



#--------------------------------------------------------------------
def coexpression(fcMAT):
    """
    Computes coexpression network from transformed expression data
    Pairwise correlation coefficient between gene expression vectors
        With a large number of clusters, the input array fcMAT has many columns; 
        to yield 'high-quality' correlation coefficients, at least 100 non nan 
        elements were required here. However, there is no problem with missing values
        once pseudo-counts are included in the fold-changes and this becomes irrelevant
    """


    N = np.shape(fcMAT)[0]
    result = np.zeros(( N, N )) * np.nan
    for i in range(N):
        current_i = fcMAT[i,:]
        for j in range(i+1, N):
            current_j = fcMAT[j,:]
            sel = (np.isnan(current_i) == False) * (np.isnan(current_j)==False)
            if np.sum(sel) > 100:       # min non-nan entries
                result[i,j] = result[j,i] = np.corrcoef(current_i[sel], current_j[sel])[0][1]

    return result







if __name__ == '__main__':


    ## LOAD DATA
    # Uniprot reference proteomes
    human_seq = SeqIO.to_dict(SeqIO.parse('../data/fasta/UP000005640_9606.fasta', "fasta"))

    # list of proteins in Uniprot reference proteomes
    list_human = [entry.split('|')[2] for entry in list(human_seq.keys())]

    # Uniprot mapping between entry name and gene name, from Uniprot.org
    uni_human = pd.read_csv("../data/uniprot/uniprot-human.tab", header=0, sep='\t', index_col=False)

    # data from Allen Brain Map
    genes_human = pd.read_csv("../data/AllenBrainAtlas/human_LGN_2021_genes-rows.csv", header=0, sep=',', index_col=False)
    exon_human =pd.read_csv("../data/AllenBrainAtlas/human_LGN_2021_exon-matrix.csv", header=0, sep=',', index_col=0)

    match_human = match_genenames(list_human, uni_human, genes_human)
    match_human.to_csv("../data/processed/id_human.txt", header=True, index=False, sep='\t')
    el_human = exon_human.loc[list(match_human['gene'])]
    el_human.to_csv("../data/processed/exons_human.txt", header=True, index=True, sep='\t')

    human_meta = pd.read_csv("../data/AllenBrainAtlas/human_LGN_2021_metadata.csv", header=0, index_col=False, sep=',')
    mca_meta = pd.read_csv("../data/AllenBrainAtlas/MCA_metadata.csv", header=0, index_col=False, sep=',')
    match_human = pd.read_csv("../data/processed/id_human.txt", header=0, index_col=False, sep='\t')


    ### MERGE DATA
    print("merge data")
    list_LGN = []
    list_LGN_cluster = []
    list_region = []
    for i in list(exon_human.columns):
        if i in list(human_meta['sample_name']):
            list_LGN.append(i)
            current_cluster = list(human_meta[human_meta['sample_name']==i]['cluster_label'])[0]
            list_LGN_cluster.append(current_cluster)
            list_region.append("LGN")
    df_LGN = exon_human[list_LGN]
    DFmeta_LGN = pd.DataFrame({'sample': list_LGN, 'cluster': list_LGN_cluster, 'region': list_region})

    # 1. load both data, normalize, merge
    list_clusters = list(set(mca_meta['cluster_label']))
    DFmeta_MCA, genes, data_MCA = extract_cluster("../data/AllenBrainAtlas/MCA_exon.csv", mca_meta)
    data_LGN = np.array(exon_human.loc[genes], dtype=int)

    DFmeta = DFmeta_MCA.append(DFmeta_LGN, ignore_index=True)
    DFmeta.to_csv("../data/processed/merged_meta.txt", header=True, index=False, sep='\t')

    data = np.concatenate((data_MCA, data_LGN), axis=1)

    genes = list( np.array(genes)[np.sum(data,1)>0] )
    data = data[np.sum(data,1)>0,:]             # omit all-zero genes

    data = normalize_countmat(data)
    
    # remove the cells that did not cluster, as the following analyses are based on the clusters
    sel_keep = DFmeta['cluster']!="NA"
    DFmeta = DFmeta.loc[sel_keep] 
    data = data[:,sel_keep]
    DFmeta.to_csv("../data/processed/merged_meta.txt", header=True, index=False, sep='\t')
    np.savetxt("../data/processed/data_merged.txt", data, fmt='%f')


    ## COMPUTE EXPRESSION AVERAGES AND NETWORKS
    # average and variable expression
    print("computer cluster averages")
    cluster_expression(DFmeta, data, genes)

    # differential expression transformation & coexpression networks
    print("compute DE transformation")
    res, fc = diffexp(data, list(DFmeta['cluster']))
    average_foldchange(fc)
    coexp = coexpression(fc)
    coexp_5pct = np.around(np.percentile(coexp[np.isnan(coexp)==False], 95), 2)
    coexp_3pct = np.around(np.percentile(coexp[np.isnan(coexp)==False], 97), 2)
    coexp_1pct = np.around(np.percentile(coexp[np.isnan(coexp)==False], 99), 2)

    print("compute coexpression networks")
    cocorrelation(coexp, coexp_5pct, genes, genes, "all_top5")
    cocorrelation(coexp, coexp_3pct, genes, genes, "all_top3")
    cocorrelation(coexp, coexp_1pct, genes, genes, "all_top1")

  
