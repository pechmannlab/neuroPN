import os, sys
import numpy as np
import pandas as pd
import networkx as nx
import itertools as it
import subprocess

from scipy.stats import mannwhitneyu

from sklearn.svm import SVC
from sklearn.metrics import accuracy_score
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import SelectFromModel

from scipy.cluster.hierarchy import dendrogram
from sklearn.cluster import AgglomerativeClustering



#--------------------------------------------------------------------
def orthologs():
    """
    Load and parse ortholog assignments
    """

    orthomat = pd.read_csv("../OMA/Output/OrthologousMatrix.txt", header=4, index_col=None, sep='\t')

    mapseqid= pd.read_csv("../OMA/Output/Map-SeqNum-ID.txt", header=None, index_col=None, skiprows=1, sep='\t')
    mapseqid.columns = ['species', 'id', 'description']

    geneDF = pd.DataFrame(columns=['human', 'macaque', 'mouse'])
    uniprotDF = pd.DataFrame(columns=['human', 'macaque', 'mouse'])

    for i in range(len(orthomat)):
        current_data = orthomat.iloc[i]
        current_human = "None"
        current_macaque = "None"
        current_mouse = "None"
        for j in range(len(current_data)):
            current_species = current_data.index[j]
            current_id = list(current_data)[j]

            current_map = mapseqid[mapseqid['species']==current_species]
            if current_id > 0 and current_id in list(current_map['id']):
                current_map = current_map[current_map['id']==current_id]
                current_description = list(current_map['description'])[0] #.split()[0]
                current_description = current_description.split()[0].split('|')[2]
                current_spec = current_description.split('_')[1]

                if current_spec == "HUMAN":
                    current_uni_human = current_description
                    if current_description in list(id_human['uniprot']):
                        current_human = list(id_human[id_human['uniprot']==current_description]['gene'])[0]
                    else:
                        current_human = "None"
                elif current_spec == "MACMU":
                    current_uni_macaque = current_description
                    if current_description in list(id_macaque['uniprot']):
                        current_macaque = list(id_macaque[id_macaque['uniprot']==current_description]['gene'])[0]
                    else:
                        current_macaque = "None"
                elif current_spec == "MOUSE":
                    current_uni_mouse = current_description
                    if current_description in list(id_mouse['uniprot']):
                        current_mouse = list(id_mouse[id_mouse['uniprot']==current_description]['gene'])[0]
                    else:
                        current_mouse = "None"
                
        geneDF.loc[len(geneDF)] = ([current_human, current_macaque, current_mouse])
        uniprotDF.loc[len(uniprotDF)] = ( [current_uni_human, current_uni_macaque, current_uni_mouse])

    geneDF.to_csv("../data/processed/ortho_genes.txt", header=True, index=False, sep='\t')
    uniprotDF.to_csv("../data/processed/orth_uniprot.txt", header=True, index=False, sep='\t')



#--------------------------------------------------------------------
def normalize_countmat(dataIN):
    """
    Normalize to CPM (counts per million)
    lazy loop through columns for large matrices
    """

    dataOUT = np.zeros(( np.shape(dataIN) ))

    for i in range(np.shape(dataIN)[1]):
        if np.nansum(dataIN[:,i]) > 0:
            dataOUT[:,i] = np.around( dataIN[:,i]/np.nansum(dataIN[:,i]) * 1000000, 4)
        else:
            dataOUT[:,i] = np.copy(dataIN[:,i])
    
    return dataOUT




#--------------------------------------------------------------------
def average_expression(dataIN):
    """
    average expression by cluster as in Seurat
    returns average and normalized average
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
def eval_clusters(expIN, metaIN, idIN):
    """
    Evaluate clustersa and perform DE transformation for species expression data
    Returns coexpression correlation matrix for input expression and annotation data
    """
    
    list_uniprot = []
    for i in list(idIN['gene']):
        if i in list(expIN.index):
            list_uniprot.append(i)
    data = expIN.loc[list_uniprot]
    #data = data.loc[ np.sum(np.array(data),1) > 0 ]        # this is dealt with during ortho assignment and consensus
    
    list_clusters = []
    for i in list(expIN.columns):
        current_cluster = list(metaIN[metaIN['sample_name']==i]['cluster_label'])[0]
        list_clusters.append(current_cluster)
    list_clusters = sorted(list(set(list_clusters)))

    list_data = []
    for i in list_clusters:
        current_samples = list(metaIN[metaIN['cluster_label']==i]['sample_name'])
        list_samples = []
        for j in current_samples:
            if j in list(data.columns):
                list_samples.append(j)
        current_data = np.array( data[list_samples] )
        #current_data_norm = normalize_countmat(current_data) # not needed
        current_data_norm = np.copy(current_data)
        list_data.append(current_data_norm)

    print( list(set(metaIN['old.class'])))

    L = len(list_clusters)
    N = np.shape(data)[0]
    result = np.zeros(( N, int((L*(L-1)/2.)) ))
    res_fc = np.zeros(( N, int((L*(L-1)/2.)) ))

    idx = 0
    for i in range(L):
        data_i = list_data[i]
        for j in range(i+1,L):
            data_j = list_data[j]
            print(i,j)
            for k in range(N):
                data_ik = data_i[k, :]
                data_jk = data_j[k, :]
                if (np.nansum(data_ik) > 0) and (np.nansum(data_jk) > 0) and np.any(data_ik != data_jk):
                    pval = mannwhitneyu(data_ik, data_jk, alternative='two-sided')[1]
                    #fc = np.nanmean(data_ik[data_ik>0]) / np.nanmean(data_jk[data_jk>0])
                    fc = np.nanmean(data_ik + 1) / np.nanmean(data_jk + 1)
                else:
                    pval = np.nan
                    fc = np.nan

                result[k, idx] = pval
                res_fc[k,idx] = np.around(fc, 4)
                
            idx += 1

    print(np.shape(res_fc))
    N = np.shape(res_fc)[0]
    coexp = np.zeros(( N, N )) * np.nan  
    for i in range(N):
        current_i = res_fc[i,:]
        for j in range(i+1, N):
            current_j = res_fc[j,:]
            sel = (np.isnan(current_i) == False) * (np.isnan(current_j)==False)
            if np.sum(sel) > 10:       # min non-nan entries for correlation coefficient, lower as smaller dataset
                coexp[i,j] = coexp[j,i] = np.around( np.corrcoef(current_i[sel], current_j[sel])[0][1], 4 )

    print(np.nansum(coexp))
    list_idx = pd.DataFrame({'gene': list(data.index)})

    return coexp, list_idx





#--------------------------------------------------------------------
def ortho_interactions(ortho, idx_hu, idx_ma, idx_mo, co_hu, co_ma, co_mo):
    """
    Function to call orthologous interactions between species coexpresssion networks
    Returns binary NxNx3 interaction matrix for N genes in 3 species
    
    Threshold for interactions is more lenient here, i.e. top 20%. This is because the 
    input data is smaller and too few interactions are conserved if too stringend here
    """

    def parse_graph(selIN, li_hu, id_hu):
        G = nx.Graph()
        for i in range(np.shape(selIN)[0]):     
            for j in range(i+1, np.shape(selIN)[0]):   
                if selIN[i,j] == 1: 
                    gene1 = li_hu[id_hu[i]]
                    gene2 = li_hu[id_hu[j]]         
                    G.add_edge(gene1, gene2)
        return G

    list_human = list(idx_hu['gene'])
    list_macaque = list(idx_ma['gene'])
    list_mouse = list(idx_mo['gene'])

    sel_none = np.sum( np.array(ortho) == 'None', 1 ) == 0
    ortho2 = ortho.iloc[sel_none]

    N = len(ortho2)
    idxDF = pd.DataFrame(columns=['human', 'macaque', 'mouse'])
    for i in range(N):
        
        i_human = ortho2.iloc[i]['human']
        if i_human in list_human:
            ix_human = list_human.index(i_human)
        else:
            ix_human = -1

        i_macaque = ortho2.iloc[i]['macaque']
        if i_macaque in list_macaque:
            ix_macaque = list_macaque.index(i_macaque)
        else: 
            ix_macaque = -1

        i_mouse = ortho2.iloc[i]['mouse']
        if i_mouse in list_mouse:
            ix_mouse = list_mouse.index(i_mouse)
        else:
            ix_mouse = -1

        if ix_human >=0 and ix_macaque >= 0 and ix_mouse >= 0:
            idxDF.loc[len(idxDF)] = (ix_human, ix_macaque, ix_mouse)

    idxMAT = np.array(idxDF)

    index_human = np.array(idxMAT[:,0], dtype=int)
    huMAT = co_hu[index_human,:]
    huMAT = huMAT[:,index_human]
    theta_human = np.percentile(huMAT[np.isnan(huMAT)==False], 80.)
    huMAT[np.isnan(huMAT)] = 0

    index_macaque = np.array(idxMAT[:,1], dtype=int)
    maMAT = co_ma[index_macaque,:]
    maMAT = maMAT[:,index_macaque]
    theta_macaque = np.percentile(maMAT[np.isnan(maMAT)==False], 80.)
    maMAT[np.isnan(maMAT)] = 0

    index_mouse = np.array(idxMAT[:,2], dtype=int)
    moMAT = co_mo[index_mouse,:]
    moMAT = moMAT[:,index_mouse]
    theta_mouse = np.percentile(moMAT[np.isnan(moMAT)==False], 80.)
    moMAT[np.isnan(moMAT)] = 0

    coexp2 = np.zeros(( np.shape(huMAT)[0], np.shape(huMAT)[1], 3))
    coexp2[:,:,0] = np.array(huMAT > theta_human, dtype=int)
    coexp2[:,:,1] = np.array(moMAT > theta_mouse, dtype=int)
    coexp2[:,:,2] = np.array(maMAT > theta_macaque, dtype=int)

    coexp_cons = (np.sum(coexp2, 2) >= 2)

    return coexp2



#--------------------------------------------------------------------
def match_clusters(exp_human, human_meta, exp_mouse, mouse_meta, exp_macaque, macaque_meta, genes_ortho):
    """
    Function to match clusters from different species based on expression profile. 
        The idea was to only use matching clusters for testing the conservation of
        interactions between species. However, using all available data and thus computing
        coexpression networks from the full data was much more informative
        - best cluster orthology by Monte Carlo (params are overkill but computatrionally cheap, so whatever)
        - by-pass option to keep all clusters
    In publication: match clusters in numbers but keep all, then DE transform
    Returns: filtered expression data for human, mouse, and macaquee for computing coexp networks
    """

    cluster_selection = False #True        # move as input var! 


    def normalize_data(dataIN):
        """
        normalize to CPM (counts per million)
        """
        dataOUT = dataIN.copy(deep=True)
        for i in list(dataOUT.columns):
            dataOUT[i] = np.around( np.array(dataOUT[i])/np.sum(np.array(dataOUT[i])) * 1000000, 4)
        return dataOUT

    def filter_clusters(metaIN):
        """
        restrict to 'Glutamatergic', 'GABAergic' b/c macaque only has those
        """
        list_output = []
        for i in list(metaIN['cluster_label']):
            current_class = list(metaIN[metaIN['cluster_label']==i]['old.class'])[0]
            if current_class in list(['Glutamatergic', 'GABAergic']):
                list_output.append(i)
        list_output = sorted(list(set(list_output)))
        return list_output


    def extract_clusters(assgnL, mta, dta, outPRE):
        sel = []
        for i in list(assgnL):
            current_list = list(mta[mta['cluster_label']==i]['sample_name'])
            if len(current_list) > 0:
                sel += current_list
        sel2 = []
        for i in sel:
            if i in list(dta.columns):
                sel2.append(i)
        output = dta[sel2]
        output.to_csv("tmp."+outPRE, header=True, index=True, sep='\t')
        return output


    def extract_averages(assgnL, mta, dta, outPRE):
        output = np.zeros(( np.shape(dta)[0], len(list(assgnL)) )) 
        list_columns = []
        for ix, i in enumerate(list(assgnL)):
            current_list = list(mta[mta['cluster_label']==i]['sample_name'])
            current_list2 = []
            for j in current_list:
                if j in list(dta.columns):
                    current_list2.append(j)
            if len(current_list2) > 0:
                list_columns.append(i)
                current_data = dta[current_list2]
                current_avrg, current_norm = average_expression(current_data)
                output[:,ix] = current_avrg
            
        output = pd.DataFrame(data=output, columns=list_columns)
        output.index = dta.index
        #output.to_csv("tmp."+outPRE, header=True, index=True, sep='\t')
        return output


    #ortho_cons = genes_ortho.iloc[ np.sum( np.array(genes_ortho)=='None', 1) == 0 ]
    # dbl-check that all genes are in indeces
    list_sel = []
    for i in range(len(genes_ortho)):
        current_cons = genes_ortho.iloc[i]
        if 'None' not in list(current_cons):
            list_sel.append(i)
    list_sel = list(set(list_sel))
    ortho_cons = genes_ortho.iloc[list_sel]

    # NORMALIZE
    norm_human = normalize_data(exp_human)
    norm_mouse = normalize_data(exp_mouse)
    norm_macaque = normalize_data(exp_macaque)

    cons_human = norm_human.loc[list(ortho_cons['human'])]
    cons_mouse = norm_mouse.loc[list(ortho_cons['mouse'])]
    cons_macaque = norm_macaque.loc[list(ortho_cons['macaque'])]

    if cluster_selection:
        list_clusters_human = filter_clusters(human_meta)
        list_clusters_mouse = filter_clusters(mouse_meta)
        list_clusters_macaque = filter_clusters(macaque_meta)
    else:
        list_clusters_human = sorted(list(set(human_meta['cluster_label'])))
        list_clusters_mouse = sorted(list(set(mouse_meta['cluster_label'])))
        list_clusters_macaque = sorted(list(set(macaque_meta['cluster_label'])))


    # average linkage
    def avrg_linkage(meta1, clusters1, exp1, meta2, clusters2, exp2):

        result = np.zeros((len(clusters1), len(clusters2) ))
        for m in range(len(clusters1)):
            idx1 = list(meta1[meta1['cluster_label']==clusters1[m]]['sample_name'])
            list_idx1 = []
            for i in idx1:
                if i in list(exp1.columns):
                    list_idx1.append(i)
            data1 = np.array(exp1[list_idx1])

            for n in range(len(clusters2)):
                idx2 = list(meta2[meta2['cluster_label']==clusters2[n]]['sample_name'])
                list_idx2 = []
                for j in idx2:
                    if j in list(exp2.columns):
                        list_idx2.append(j)
                data2 = np.array(exp2[list_idx2])

                res = []
                for ix in range(np.shape(data1)[1]):
                    for jx in range(np.shape(data2)[1]):
                        sel = (np.isnan(data1[:,ix]) == False) * (np.isnan(data2[:,jx]) == False)
                        corr = np.corrcoef(data1[sel,ix], data2[sel,jx])
                        res.append(corr[0,1])
                res = np.array(res)   
                result[m,n] = np.mean(res)

        return result

    res1 = avrg_linkage(human_meta, list_clusters_human, cons_human, mouse_meta, list_clusters_mouse, cons_mouse)
    res2 = avrg_linkage(human_meta, list_clusters_human, cons_human, macaque_meta, list_clusters_macaque, cons_macaque)
    res3 = avrg_linkage(macaque_meta, list_clusters_macaque, cons_macaque, mouse_meta, list_clusters_mouse, cons_mouse)

    N = int( np.min([len(list_clusters_human), len(list_clusters_mouse), len(list_clusters_macaque)]) )

    trails_max = 10000
    current_trail = 0
    best_score = 0
    best_idx = None
    while current_trail < trails_max:
        current_trail += 1

        idxMAT = np.zeros(( 3, N ), dtype=int)
        idxMAT[0,:] = np.random.choice( np.arange(len(list_clusters_human)), N, replace=False)
        idxMAT[1,:] = np.random.choice( np.arange(len(list_clusters_mouse)), N, replace=False)
        idxMAT[2,:] = np.random.choice( np.arange(len(list_clusters_macaque)), N, replace=False)

        score = np.nansum(res1[idxMAT[0,:], idxMAT[1,:]]) + np.nansum(res2[idxMAT[0,:], idxMAT[2,:]]) + np.nansum(res3[idxMAT[2,:], idxMAT[1,:]])

        i_max = 1000
        i = 0
        while i < i_max:
            i += 1
            current_idxMAT = np.copy(idxMAT)
            current_species = np.random.choice( np.arange(3), 1)
            current_swap = np.random.choice( np.arange(N), 2, replace=False)
            current_swap0 = current_idxMAT[current_species, current_swap[0]]
            current_swap1 = current_idxMAT[current_species, current_swap[1]]
            current_idxMAT[current_species, current_swap[0]] = current_swap1
            current_idxMAT[current_species, current_swap[1]] = current_swap0
            current_score = np.nansum(res1[current_idxMAT[0,:], current_idxMAT[1,:]]) + np.nansum(res2[current_idxMAT[0,:], current_idxMAT[2,:]]) + np.nansum(res3[current_idxMAT[2,:], current_idxMAT[1,:]])

            if current_score > score:        
                idxMAT = np.copy(current_idxMAT)
                score = np.copy(current_score)

        if score > best_score:
            best_score = np.copy(score)
            best_idx = np.copy(idxMAT)

    print(best_score)
    print(best_idx)
    best_ind = np.zeros(( N ))
    for i in range(N):
        best_ind[i] =  np.nansum([res1[best_idx[0,i], best_idx[1,i]], res2[best_idx[0,i], best_idx[2,i]], res3[best_idx[2,i], best_idx[1,i]] ])
    print(best_ind)

    # only keep clusters that actually match, cumulative correlation above threshold
    # all for DE, selection for cluster EL
    sel_score = best_ind > 0      
    sel_el = best_ind > 0.8
    best_idx = best_idx[:,sel_score]
    el_idx = best_idx[:,sel_el]
    print(el_idx)

    # all clusters for DE-transfor coexpression 
    assign_human = list( np.array(list_clusters_human)[ np.array(best_idx[0,:], dtype=int)] )
    assign_mouse = list( np.array(list_clusters_mouse)[ np.array(best_idx[1,:], dtype=int)] )
    assign_macaque = list( np.array(list_clusters_macaque)[ np.array(best_idx[2,:], dtype=int)] )
    assignment = pd.DataFrame({'human': list(assign_human), 'mouse': list(assign_mouse), 'macaque': list(assign_macaque)})

    ortho_human = extract_clusters(list(assignment['human']), human_meta, cons_human, "human")
    ortho_mouse = extract_clusters(list(assignment['mouse']), mouse_meta, cons_mouse, "mouse")
    ortho_macaque = extract_clusters(list(assignment['macaque']), macaque_meta, cons_macaque, "macaque")

    # matching clusters for EL conservation
    assel_human = list( np.array(list_clusters_human)[ np.array(el_idx[0,:], dtype=int)] )
    assel_mouse = list( np.array(list_clusters_mouse)[ np.array(el_idx[1,:], dtype=int)] )
    assel_macaque = list( np.array(list_clusters_macaque)[ np.array(el_idx[2,:], dtype=int)] )
    assignment_el = pd.DataFrame({'human': list(assel_human), 'mouse': list(assel_mouse), 'macaque': list(assel_macaque)})

    orthel_human = extract_averages(list(assignment_el['human']), human_meta, cons_human, "human")
    orthel_mouse = extract_averages(list(assignment_el['mouse']), mouse_meta, cons_mouse, "mouse")
    orthel_macaque = extract_averages(list(assignment_el['macaque']), macaque_meta, cons_macaque, "macaque")

    orthel_human.to_csv("../data/processed/coevo/ortho_el_human.txt", header=True, index=True, sep='\t')
    orthel_mouse .to_csv("../data/processed/coevo/ortho_el_mouse.txt", header=True, index=True, sep='\t')
    orthel_macaque.to_csv("../data/processed/coevo/ortho_el_macaque.txt", header=True, index=True, sep='\t')

    return ortho_human, ortho_mouse, ortho_macaque



#--------------------------------------------------------------------
def get_indices(chapIN, ubIN, synIN, riboIN, listID):
    """
    Generate index DF of input genes of interest corresponding to
    indices in orthologous expression matrices 
    """

    indexDF = pd.DataFrame(columns=['Gene', 'index', 'cat', 'system'])
    
    for i in list(chapIN['Gene']):
        if i in listID:
            current_idx = listID.index(i)
            current_cat = "chap"
            current_sys = list(chapIN[chapIN['Gene']==i]['System'])[0]
            indexDF.loc[len(indexDF)] = (i, current_idx, current_cat, current_sys)

    for i in list(ubIN['Gene']):
        if i in listID:
            current_idx = listID.index(i)
            current_cat = "ub"
            current_sys = list(ubIN[ubIN['Gene']==i]['Category'])[0]
            indexDF.loc[len(indexDF)] = (i, current_idx, current_cat, current_sys)

    for i in list(synIN):
        if i in listID:
            current_idx = listID.index(i)
            current_cat = "syn"
            current_sys = "None" 
            indexDF.loc[len(indexDF)] = (i, current_idx, current_cat, current_sys)

    for i in list(riboIN):
        if i in listID:
            current_idx = listID.index(i)
            current_cat = "ribo"
            current_sys = "None" 
            indexDF.loc[len(indexDF)] = (i, current_idx, current_cat, current_sys)

    return indexDF



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def randomized_network(data1, data2, data3):
    """
    Wrapper function to call the R script 'randomize_network.R':
        generates randomized network with same degree distribution
        wrapper for BiRewire R package
        output is written to file by RScript
    """

    np.savetxt("tmp.inp1", data1, fmt='%i')
    np.savetxt("tmp.inp2", data2, fmt='%i')
    np.savetxt("tmp.inp3", data3, fmt='%i')

    cmd_rand = "Rscript randomize_network.R"
    output = subprocess.run(cmd_rand, shell=True)


    os.remove("tmp.inp1")
    os.remove("tmp.inp2")
    os.remove("tmp.inp3")

    return None



#--------------------------------------------------------------------
def eval_density(consHU, consMO, consMA, listIDX, POST, N=1000):
    """
    Compares network density (true interactions vs. possible interactions) in 
    species and conserved networks
    """

    cons3 = (consHU + cons_mo + cons_ma) == 3

    rand_hu = np.zeros(( N ))
    rand_mo = np.zeros(( N ))
    rand_ma = np.zeros(( N ))
    rand_cons = np.zeros(( N ))

    for i in range(N):
        idx_rand = np.random.choice(np.arange(len(list_id)), len(listIDX))

        rand_hu[i] = np.sum( consHU[idx_rand,:][:,idx_rand])/np.prod(np.shape(consHU[idx_rand,:][:,idx_rand])) 
        rand_mo[i] = np.sum( consMO[idx_rand,:][:,idx_rand])/np.prod(np.shape(consMO[idx_rand,:][:,idx_rand])) 
        rand_ma[i] = np.sum( consMA[idx_rand,:][:,idx_rand])/np.prod(np.shape(consMA[idx_rand,:][:,idx_rand])) 
        rand_cons[i] = np.sum( cons3[idx_rand,:][:,idx_rand])/np.prod(np.shape(cons3[idx_rand,:][:,idx_rand])) 

    true_hu = np.sum( consHU[listIDX,:][:,listIDX])/np.prod(np.shape(consHU[listIDX,:][:,listIDX])) 
    true_mo = np.sum( consMO[listIDX,:][:,listIDX])/np.prod(np.shape(consMO[listIDX,:][:,listIDX])) 
    true_ma = np.sum( consMA[listIDX,:][:,listIDX])/np.prod(np.shape(consMA[listIDX,:][:,listIDX])) 
    true_cons = np.sum( cons3[listIDX,:][:,listIDX])/np.prod(np.shape(cons3[listIDX,:][:,listIDX])) 
   
    np.savetxt("../data/processed/coevo/cons_rand_"+POST+".txt", rand_cons)

    return true_cons, np.mean(rand_cons), true_cons/np.mean(rand_cons)   


 
#--------------------------------------------------------------------
def parse_degree(consHU, consMO, consMA, listID, chapIN, ubIN):
    """
    Compares network degree per gene in true conserved and randomized conserved networks
    """
   
    cons3 = (consHU + cons_mo + cons_ma) == 3

    rand_degree = np.loadtxt("../data/processed/coevo/cons_rand_degree.txt")
    rand_95 = np.zeros(( np.shape(rand_degree)[0] ))
    colo = np.zeros(( np.shape(rand_degree)[0] ))
    cat = []
    for i in range( np.shape(rand_degree)[0] ):
        rand_95[i] = np.percentile(rand_degree[i,:], 95)

        if listID[i] in chapIN:
            current_cat = "chap"
        elif listID[i] in ubIN:
            current_cat = "ub"
        else:
            current_cat = "None"
        cat.append(current_cat)
    

    degDF = pd.DataFrame({'cons': list(np.sum(cons3,1)), 
                          'rand': list(np.mean(rand_degree,1)), 
                          'q95': list(rand_95),
                          'cat': list(cat)}, index=list(listID) )
    degDF.to_csv("../data/processed/coevo/degree.txt", header=True, index=True, sep='\t')





if __name__ == '__main__':


    # LOAD DATA
    # expression data
    lgn_human = pd.read_csv("../data/AllenBrainAtlas/human_LGN_2021_exon-matrix.csv", header=0, index_col=0, sep=',')
    lgn_macaque = pd.read_csv("../data/AllenBrainAtlas/macaque_LGN_2021_exon-matrix.csv", header=0, index_col=0, sep=',')
    lgn_mouse = pd.read_csv("../data/AllenBrainAtlas/mouse_LGN_2021_exon-matrix.csv", header=0, index_col=0, sep=',')
    # meta data
    meta_human = pd.read_csv("../data/AllenBrainAtlas/human_LGN_2021_metadata.csv", header=0, index_col=0, sep=',')
    meta_macaque =  pd.read_csv("../data/AllenBrainAtlas/macaque_LGN_2021_metadata.csv", header=0, index_col=0, sep=',')
    meta_mouse = pd.read_csv("../data/AllenBrainAtlas/mouse_LGN_2021_metadata.csv", header=0, index_col=0, sep=',')
    # gene ids
    id_human = pd.read_csv("../data/processed/id_human.txt", header=0, index_col=False, sep='\t')
    id_macaque = pd.read_csv("../data/processed/id_macaque.txt", header=0, index_col=False, sep='\t')
    id_mouse = pd.read_csv("../data/processed/id_mouse.txt", header=0, index_col=False, sep='\t')


    # extract sets of orthologous genes
    if not os.path.exists("../data/processed/ortho_genes.txt"):
        orthologs()
    ortho_genes = pd.read_csv("../data/processed/ortho_genes.txt", header=0, index_col=False, sep='\t')

    # identify best sets of orthologous clusters
    ortho_human, ortho_mouse, ortho_macaque = match_clusters(lgn_human, meta_human, lgn_mouse, meta_mouse, lgn_macaque, meta_macaque, ortho_genes)
 
    # compute species specific coexpression networks and find conserved interactions
    coexp_human, idx_human = eval_clusters(ortho_human, meta_human, id_human)
    coexp_mouse, idx_mouse = eval_clusters(ortho_mouse, meta_mouse, id_mouse)
    coexp_macaque, idx_macaque = eval_clusters(ortho_macaque, meta_macaque, id_macaque)
  
    coexp = ortho_interactions(ortho_genes, idx_human, idx_macaque, idx_mouse, coexp_human, coexp_macaque, coexp_mouse)
    coexp2 = np.sum( coexp, 2)

    cons_hu = coexp[:,:,0]      # can keep condensed matrix but more cryptic
    cons_mo = coexp[:,:,1]
    cons_ma = coexp[:,:,2]

    anno_chap = pd.read_csv("../data/auxilliary/list_chaperones.txt", header=0, index_col=False, sep=' ')
    anno_ub = pd.read_csv("../data/ub/ub.csv", header=0, index_col=False, sep='\t')
    anno_syn = pd.read_csv("../data/auxilliary/synapse/synapse_proteome.txt", header=0, index_col=False, sep='\t')
    anno_channels = pd.read_csv("../data/processed/list_channels.txt", header=None, index_col=False)

    idx_human = pd.read_csv("tmp.idx.human", header=0, index_col=None)

    ribo = pd.read_csv("human.ribosome", header=None, index_col=None, sep=',')
    list_ribo = list(ribo[0])

    list_chap = list(anno_chap['Gene'])
    list_ub = list(anno_ub['Gene'])
    list_channels = list(anno_channels[0])
    list_syn = list(set(anno_syn['symbol']))
    list_synapse = []                       # make synapse stuff all other than chap and ub
    for i in list_syn:
        if i not in list_channels and i not in list_chap and i not in list_ub:
            list_synapse.append(i)
    list_synapse += list_channels # synapse to include channels

    list_id = list(idx_human['gene'])

    idxDF = get_indices(anno_chap, anno_ub, list_synapse, list_ribo, list_id)
    idxDF.to_csv("tmp.idx", header=True, index=False, sep='\t')



    #randomized_network(cons_hu, cons_mo, cons_ma)
    parse_degree(cons_hu, cons_mo, cons_ma, list_id, list_chap, list_ub)

    resDF = pd.DataFrame(columns=['true', 'rand', 'ratio'])
    resDF.loc[len(resDF)] = eval_density(cons_hu, cons_mo, cons_ma, list(idxDF[idxDF['cat']=="chap"]['index']), "chap")
    resDF.loc[len(resDF)] = eval_density(cons_hu, cons_mo, cons_ma, list(idxDF[idxDF['cat']=="ub"]['index']), "ub")
    resDF.loc[len(resDF)] = eval_density(cons_hu, cons_mo, cons_ma, list(idxDF[idxDF['cat']=="syn"]['index']), "syn")

    res_ctrl = np.zeros(( 100, 3 ))
    for i in range(100):
        idx_ctrl = np.random.choice(np.arange(len(list_id)), 1000)
        res_ctrl[i,:] = eval_density(cons_hu, cons_mo, cons_ma, idx_ctrl, "control", N=1000)

    resDF.loc[len(resDF)] = np.mean(res_ctrl, 0)        #eval_density(idx_ctrl)
    resDF.index = list(['chap', 'ub', 'syn', 'ctrl'])
    resDF.to_csv("../data/processed/coevo/density.txt", header=True, index=True, sep='\t')
    print(resDF)