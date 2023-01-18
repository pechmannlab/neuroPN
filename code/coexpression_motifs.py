import os, sys
import numpy as np
import pandas as pd
import networkx as nx
import itertools as it
import time
import multiprocessing as mp
import glob
import subprocess
from scipy.stats import mannwhitneyu



# Some code is recycled/modified from [Sahoo & Pechmann, PeerJ, 2022]
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def partition(n):
    """
    partition of integer n
    """

    result = []
    result.append((n, ))
    for x in range(1, n):
        for y in partition(n - x):
            result.append( tuple(sorted((x, ) + y)))

    return result

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def composition(n):
    """
    computes combinations of partition of integer n
    """

    composition = []
    pat = partition(n)
    for p in pat:
        for comp in list(set(list(it.permutations(p)))):
            if comp not in composition:
                composition.append(comp)

    return composition

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def parse_motif(motif_list):
    """
    parses motif, flattens potentially nested list
    """

    motif = []
    for i in motif_list:
        if np.size(i) == 1:
            motif.append(i)
        elif np.size(i) > 1:
            for j in i:
                motif.append(j)
    #motif = list(set(motif))

    return motif 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def to_g6(motif_array):
    """
    extract motif and convert to graph6 format
    """
    
    g_motif = nx.from_numpy_array(motif_array)
    g6 = nx.to_graph6_bytes(g_motif, nodes=None, header=False)

    g6 = g6.rstrip().decode('utf-8')

    return g6

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def rename_graph(G, labels=[]):
    """
    renames nodes in graph to integers
    returns renamed graph and dictionary with node mappings
    """

    if len(labels) == 0:
        label_dict = {}
        label_rev = {}
        label_num = 0 
        for node in G.nodes():
            label_dict[node] = label_num
            label_rev[label_num] = node
            label_num += 1
    else:
        label_dict = labels 
        label_rev = {}
        for i in list(label_dict.keys()):
            label_rev[label_dict[i]] = i

    Gnew = nx.relabel_nodes(G, label_dict, copy=True)

    return Gnew, label_dict, label_rev
 


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def enumerate_kmotif(G_mat, S0, k, labels_rev, TMPDIR, writeORF=True):
    """
    Implementation of Kavosh motif enumeration (Kashani et al. BMC Bioinfo, 2009) with own modifications
    G_mat: input graph (PPI) in np.array form
    S0: starting 'source' node as index for graph arrays
    k: size of motifs to be enumerated
    labels_rev: dictionary to link matrix indices to ORF names
    """

    if k > 6:
        print("do you really want this? check your hardware first")
        return None

    list_graphlet = composition(k-1)

    discovered_graphlets = {} 
    num_discovered = 0

    resultcols = ['ORF'+str(i+1) for i in range(k)]   
    resultDF = pd.DataFrame(columns=resultcols)

    num = 0
    fnum = 0

    for g in range(len(list_graphlet)):
        current_graphlet = list_graphlet[g]
        num_levels = np.size(current_graphlet) 
        current_level = 0

        visited = np.zeros(( np.shape(G_mat)[0]  ), dtype=bool)
        visited[S0] = True

        queue = []
        queue.append([S0])
   
        for i in range(np.size(current_graphlet)):
            queue.append([])

        read_pointer = np.zeros(( num_levels + 1 ), dtype=int)
        len_queue = [len(queue[x]) for x in range(len(queue))]

        current_motif = []
        s = 1
     
        while np.any(read_pointer < (len_queue) ): 
            if np.size(current_graphlet) > 1: 
                if current_level < np.size(current_graphlet):
                    s = current_graphlet[current_level]
                else:
                    s = 1 #placeholder, though not needed
            elif np.size(current_graphlet) == 1:
                s = int(current_graphlet[0])

            if read_pointer[current_level] < len(queue[current_level]):  
                current_node = queue[current_level][read_pointer[current_level]]
                read_pointer[current_level] += 1
                current_motif.append(current_node)

                if current_level < num_levels:
                    if s == 1:
                        if np.size(current_node) == 1:
                            current_parent = current_node
                        elif np.size(current_node ) > 1:
                            current_parent = current_node[0] 
                        current_children = np.where(G_mat[current_parent,:] > 0)[0] 
                        current_children = current_children[ visited[current_children] == False ] 
                        visited[current_children] = True

                        queue[int(current_level+1)] += list(current_children)
                        current_level += 1               

                    elif s > 1:
                        if np.size(current_node) == 1:
                            current_parent = current_node
                        elif np.size(current_node ) > 1:
                            current_parent = current_node[0] 
                        current_children = np.where( G_mat[current_parent,:] > 0)[0] 
                        current_children = current_children[ visited[current_children] == False ] 
                        #visited[current_children] = True ## see below

                        if current_level < num_levels - 1:
                            current_combs = [list(x) for x in it.permutations(current_children, int(s) ) ]
                        elif current_level == num_levels -1:
                            current_combs = list(it.combinations(current_children,int(s) ))

                     
                        current_combs = list(current_combs)
                        queue[int(current_level+1)] += current_combs
                        visited[current_children] = True 
                     
                        current_level += 1               

                elif current_level == num_levels:
                    result_motif = parse_motif(current_motif)
                    current_motif.pop()

                    result_orf = []
                    for node in result_motif:
                        result_orf.append(labels_rev[node])
                        output_motif = [result_orf[x] for x in range(len(result_orf))] 

                    if writeORF:
                        resultDF.loc[len(resultDF)] = sorted(output_motif) 
                    else:
                        resultDF.loc[len(resultDF)] = sorted(result_motif) 
   
                    g_array = G_mat[result_motif,:][:,result_motif]

                    current_g6 = to_g6(g_array)
                    if current_g6 in list(discovered_graphlets.keys()):
                        discovered_graphlets[current_g6] += 1
                    else:
                        discovered_graphlets[current_g6] = 1

                    num += 1
           
            elif read_pointer[current_level] == len(queue[current_level]):
                list_current_level = np.array(queue[current_level], dtype=int)
                visited[list_current_level] = False 
                queue[current_level] = []
                read_pointer[current_level] = 0
                current_motif.pop()
                current_level -= 1 
            
            len_queue = [len(queue[x]) for x in range(len(queue))]

            if np.all(read_pointer == len_queue):   
                break

    print(S0, num)
    resultDF.to_csv(TMPDIR+str(S0)+'_motif.txt', header=True, index=False, sep='\t')

    return None



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def find_motifs(G, k, listINP, Fout_pre, writeORF):
    """
    wrapper function to enumerate motifs
    goes through each node as source (i.e. takes time to compute)
    G:  PPI network in nx format
    GI: GI network in nx format
    k:  size of motif
    Fout_pre: prefix for output files
    randomized: logical whether PPI network is randomized, using R package Birewire
    randomized_gi: logical whether GI network is randomized
    theta_func: threshold for fraction of GI edges in motif to consider 'functional'
    writeDF: write extended output to disk
    """

    TMPDIR = 'TMP_' + Fout_pre + '/'
    if not os.path.exists(TMPDIR):
        os.mkdir(TMPDIR)

    G, labels, labels_rev = rename_graph(G)

    G_mat = nx.to_numpy_array(G, nodelist=list(np.arange(len(labels))) )
    G_mat = np.array(G_mat > 0, dtype=int)      # force unweighted
    np.fill_diagonal(G_mat, 0)


    iter_max = np.shape(G_mat)[0]
    # for restart omit what already computed
    #iter_list = list(np.arange(iter_max))

    iter_list = []
    for i in listINP:
        if i in list(labels.keys()):
            iter_list.append( labels.get(i) )

    pool = mp.Pool(processes=20)        # parallelize computation
    for s in iter_list:
        pool.apply_async(motif_parallel, [(s, G_mat, k, labels_rev, TMPDIR, writeORF)] ) 
    pool.close()
    pool.join()

    # COLLECT DATA FROM INDIVIDUAL FILES (SIMPLER FOR RESTART)  
    results = []
    for f in glob.glob(TMPDIR + '*_motif.txt'):
        df = pd.read_csv(f, header=0, index_col=False, sep='\t')
        results.append(df)        
    resultDF = pd.concat(results, ignore_index=True, sort=True)
    resultDF = resultDF.drop_duplicates()
    resultDF.to_csv('motif_' + Fout_pre +'.txt', header=True, index=False, sep='\t')

    return None



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def motif_parallel(packaged_input):
    """
    import/export wrapper to run motifs on multi-cores
    """

    # pack/unpack vars for muliprocessing
    source = packaged_input[0]
    G_mat = packaged_input[1]
    k = packaged_input[2]
    labels_rev = packaged_input[3]
    TMPDIR = packaged_input[4]
    writeORF = packaged_input[5]

    enumerate_kmotif(G_mat, source, k, labels_rev, TMPDIR, writeORF)
 
    return None



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def enriched_selection(dataIN, sel1, sel2, listINP, listIN):
    """
    Identify genes expressed higher in sel1 than sel2
    Hardcoded threshold 20%, i.e. fc > 1.2
    """
    N = np.shape(avrg_data)[0]
    result = np.zeros(( N ))
    fc = np.zeros(( N ))
    for i in range(N):
        data_i = avrg_data[i, sel1]
        data_j = avrg_data[i, sel2]
        if (np.sum(data_i) + np.sum(data_j) > 0) and np.any(data_i != data_j):
            pval = mannwhitneyu(data_i, data_j, alternative='two-sided')[1]
        else:
            pval = np.nan
        result[i] = pval
        fc[i] = np.mean(data_i + 1) / np.mean(data_j + 1)      # pseudo count of 1, like Seurat
         
    sel_orf = avrg_exp.index[ fc > 1.2 ]

    sel_orf2 = []
    sel_pn = []
    for i in sel_orf:
        if i in listIN:
            sel_orf2.append(i)
        if i in listINP:
            sel_pn.append(i)

    return sel_pn, sel_orf2



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def cluster_motifs(G_IN, motifIN, DIR):
    """
    Parses network of motifs
    Finds modules in network with 
    nx.algorithms.community.greedy_modularity_communities(G)
    Output is lists of genes in network modules written to file
    """

    CLUSTERDIR = 'clusters_' + DIR
 
    # LOAD DATA
    list_orfs = []
    N = 0
    motif_data = pd.read_csv(motifIN, header=0, index_col=False, sep='\t')
    current_columns = motif_data.columns
    current_orfs = list(current_columns)
    N += len(motif_data)
    for i in current_orfs:
        list_orfs += list(set(list(motif_data[i])))
    list_orfs = sorted(list(set(list_orfs)))           # make ordered so that following results are deterministic

    G = nx.subgraph(G_IN, list_orfs)  
    modules = nx.algorithms.community.greedy_modularity_communities(G)

    for ix, i in enumerate(modules):
        if len(i) > 20:                      # no orphan 'modules'
            gsub = nx.subgraph(G_IN, i)
            gmat = nx.to_numpy_array(gsub, nodelist=list(i) )
            gmat = np.array(gmat > 0, dtype=int)
            np.savetxt("../data/processed/clusters/"+DIR+"/cluster."+str(ix)+".txt", gmat, fmt='%i')
            nodesDF = pd.DataFrame({'ORF': list(i)})
            nodesDF.to_csv("../data/processed/clusters/"+DIR+"/cluster."+str(ix)+"_orfs.txt", header=True, index=False)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def parse_networks(G_IN, dirIN, listCHAP, listUB, prefixOUT, chapDF, ubDF):
    """
    join and parse modules into network viz:
    - load modules and parse into graph, edges labeled by module
    - load the rest of the base network, edges labeled as such
    - write files for plotting with ggraph
    """


    # load clusters
    networkDF = pd.DataFrame(columns=['from', 'to', 'cluster'])
    nodesDF = pd.DataFrame(columns=['Gene', 'cat', 'cluster'])
    G_clust = nx.Graph()
    G_pn = nx.Graph()
    n = int( len( list(glob.glob(dirIN + 'cluster.*.txt') )) / 2. )
    
    for i in range(n):
        current_data = np.loadtxt(dirIN+"cluster."+str(i)+".txt")
        current_orfs = pd.read_csv(dirIN+"cluster."+str(i)+"_orfs.txt")

        for jx, j in enumerate(list(current_orfs['ORF'])):
            for kx, k in enumerate(list(current_orfs['ORF'])):
                if kx > jx and current_data[jx, kx] == 1:
                    G_clust.add_edge(j, k)
                    if j in listCHAP:
                        j_anno = "chap"
                    elif j in listUB:
                        j_anno = "ub"
                    else:
                        j_anno = "syn"

                    if k in listCHAP:
                        k_anno = "chap"
                    elif k in listUB:
                        k_anno = "ub"
                    else:
                        k_anno = "syn"

                    if j_anno != "syn" and k_anno == "syn":
                        networkDF.loc[len(networkDF)] = (j,  k,  "cluster"+str(i))
                        G_pn.add_edge(j,k)
                    elif j_anno == "syn" and k_anno != "syn":
                        networkDF.loc[len(networkDF)] = (k,  j,  "cluster"+str(i))
                        G_pn.add_edge(j,k)
                    else:
                        networkDF.loc[len(networkDF)] = (j, k, "cluster"+str(i))

                    if j not in list(nodesDF['Gene']):
                        nodesDF.loc[len(nodesDF)] = (j, j_anno, "cluster"+str(i))
                    if k not in list(nodesDF['Gene']):
                        nodesDF.loc[len(nodesDF)] = (k, k_anno, "cluster"+str(i))

    for i,j in G_IN.edges():
        if i in list(nodesDF['Gene']) and j in list(nodesDF['Gene']):    
            networkDF.loc[len(networkDF)] = (i, j, "nocluster")

    networkDF.to_csv("../data/processed/clusters/"+prefixOUT+"_network.txt", header=True, index=False, sep='\t')

    deg = dict(G_pn.degree())
    clust = nx.clustering(G_pn)
    list_degree = []
    list_clust = []
    list_sys = []
    for i in list(nodesDF['Gene']):
        current_degree = deg.get(i, 0)
        list_degree.append(current_degree)
        current_clust = clust.get(i,0)
        list_clust.append(current_clust)
        if i in list(chapDF['Gene']):
            current_sys = list(chapDF[chapDF['Gene']==i]['System'])[0]
        elif i in list(ub['Gene']):
            current_sys = list(ub[ub['Gene']==i]['Category'])[0]
        else:
            current_sys = "None"
        list_sys.append(current_sys)
    nodesDF['degree'] = list_degree
    nodesDF['PNsyst'] = list_sys
    #nodesDF['clustcoeff'] = list_clust

    nodesDF.to_csv("../data/processed/clusters/"+prefixOUT+"_annotation.txt", header=True, index=False, sep='\t')
  
    return nodesDF



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def pn_network(G_IN, chapDF, ubDF):
    """
    Subgraph of interactions that involve at least one protein homeostasis network (PN) component
    Parses the degree of PN genes to other stuff into DF
    """

    list_pn = list(chapDF['Gene']) + list(ubDF['Gene'])

    G_pn = nx.Graph()
    for i in G_IN.edges():
        n1, n2 = i[0], i[1]
        if n1 in list_pn and n2 not in list_pn:
            G_pn.add_edge(n1, n2)
        elif n1 not in list_pn and n2 in list_pn:
            G_pn.add_edge(n1, n2)

    list_nodes = sorted(list(G_pn.nodes()))
    deg = dict(G_pn.degree())
    list_degree = []
    list_sys = []
    for i in list_nodes:
        current_degree = deg.get(i, 0)
        list_degree.append(current_degree)
       
        if i in list(chapDF['Gene']):
            current_sys = list(chapDF[chapDF['Gene']==i]['System'])[0]
        elif i in list(ubDF['Gene']):
            current_sys = list(ubDF[ubDF['Gene']==i]['Category'])[0]
        else:
            current_sys = "None"
        list_sys.append(current_sys)


    nodesDF = pd.DataFrame({'Gene': list_nodes, 'degree': list_degree, 'PNsyst': list_sys})

    return nodesDF



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def reldegree(G_IN, G_ref, chapDF, ubDF):
    """
    Relative connectivity of PN components in coexpression motif networks
    """

    resDF_true = pn_network(G_IN, chap, ub)
    print(resDF_true)

    list_cat = list(set( list(chapDF['System']) + list(ubDF['Category']) ))
    dict_cat = {}

    N = 1000
    res = np.zeros(( N ))
    list_res = []
    for i in range( N ):

        list_sample = list( np.random.choice( np.array(list(G_ref.nodes())), len(list(G_IN.nodes())) ) )
        G_sample = nx.subgraph(G_ref, list_sample)
        resDF = pn_network(G_sample,  chap, ub)
       
        #print(resDF)
        for j in list_cat:
            if j in list(resDF['PNsyst']):
                current_sys = np.mean( np.array(resDF[resDF['PNsyst']==j]['degree']) )
                if j in dict_cat:
                    dict_cat[j].append(current_sys)
                else:
                    dict_cat[j] = [current_sys]

        res[i] = np.mean(np.array(resDF_true['degree'])) / np.mean(np.array(resDF['degree']))
        list_res += list(resDF['degree'])

    pval = mannwhitneyu(np.array(resDF_true['degree']), np.array(list_res), alternative='two-sided')[1]
    result_degree = np.array([np.mean(res), np.std(res), pval])

    catDF = pd.DataFrame(columns=['System', 'enrich', 'pval'])
    for i in list(set(list(resDF_true['PNsyst']))):
        if i != "None":
            current_sys = np.mean( np.array(resDF_true[resDF_true['PNsyst']==i]['degree']) )
            current_ref = np.mean( np.array(dict_cat[i]) )
            current_ratio = current_sys/current_ref
            current_pval = mannwhitneyu(np.array(resDF_true[resDF_true['PNsyst']==i]['degree']), np.array(dict_cat[i]), alternative='two-sided')[1]

            catDF.loc[len(catDF)] = (i, current_ratio, current_pval)
  
    return result_degree, catDF








if __name__ == '__main__':


    ## LOAD DATA
    avrg_exp = pd.read_csv("../data/processed/cluster_average_norm.txt", header=0, index_col=0, sep='\t')
    meta = pd.read_csv("../data/processed/meta_types_clusters.txt", header=0, index_col=False, sep='\t')

    channels = pd.read_csv("../data/processed/list_channels.txt", header=None, index_col=False)
    syn = pd.read_csv("../data/auxilliary/synapse/synapse_proteome.txt", header=0, index_col=False, sep='\t')
    chap = pd.read_csv("../data/auxilliary/list_chaperones.txt", header=0, index_col=False, sep=' ')
    ub = pd.read_csv("../data/auxilliary/ub.csv", header=0, index_col=False, sep='\t')
 
    list_channels = list(set(channels[0]))
    list_chap = list(set(chap['Gene']))
    list_ub = list(set(ub['Gene']))
    list_synapse = []                       # make synapse stuff all other than chap and ub
    for i in list(set(syn['symbol'])):
        if i not in list_channels and i not in list_chap and i not in list_ub:
            list_synapse.append(i)
    list_synapse += list_channels # synapse to include channels
    list_all = list_synapse + list_channels + list_chap + list_ub
    list_pn = list_chap + list_ub


    # cell types for Exc, Inh, and non
    list_avrg = list(avrg_exp.columns)
    sel_exc = np.zeros(( len(list_avrg) ), dtype=bool)
    sel_inh = np.zeros(( len(list_avrg) ), dtype=bool)
    sel_non = np.zeros(( len(list_avrg) ), dtype=bool)

    for ix, i in enumerate(list_avrg):
        current_type = list(meta[meta['cluster']==i]['type'])[0]
        if current_type == "Exc":
            sel_exc[ix] = True
        elif current_type == "Inh":
            sel_inh[ix] = True
        else:
            sel_non[ix] = True


    # expression data and expression enriched genes between sets
    avrg_data = np.array(avrg_exp)

    inp_exc_inh, sel_exc_inh = enriched_selection(avrg_exp, sel_exc, sel_inh, list_pn, list_all)
    inp_inh_exc, sel_inh_exc = enriched_selection(avrg_exp, sel_inh, sel_exc, list_pn, list_all)

    inp_exc_non, sel_exc_non = enriched_selection(avrg_exp, sel_exc, sel_non, list_pn, list_all)
    inp_inh_non, sel_inh_non = enriched_selection(avrg_exp, sel_inh, sel_non, list_pn, list_all)
    sel_neuron_non = []
    inp_neuron_non = []
    list_nn = list(set(list(sel_exc_non) + list(sel_inh_non)))
    for i in list_nn:
        if i in sel_exc_non and i in sel_inh_non:
            sel_neuron_non.append(i)
        if i in list_pn:
            inp_neuron_non.append(i)



    # load network data and find motifs in subgraphs

    print("loading network data")
    G = nx.read_gpickle("../data/processed/coexp_all_top5.nx")

    G_syn = nx.subgraph(G, list_all)
    G_exc_inh = nx.subgraph(G, sel_exc_inh)
    G_inh_exc = nx.subgraph(G, sel_inh_exc)
    G_neu_non = nx.subgraph(G, sel_neuron_non)

 
    start_time = time.time()
    print("enumerating motifs")
    for k in [3]:
        print("now computing k =", k)
        find_motifs(G_syn, k, list_pn, "syn_"+str(k), True )
        find_motifs(G_exc_inh, k, inp_exc_inh, "exc_inh_"+str(k), True )
        find_motifs(G_inh_exc, k, inp_inh_exc, "inh_exc_"+str(k), True )
        find_motifs(G_neu_non, k, inp_neuron_non, "neu_non_"+str(k), True )
    print("--- %s seconds --- " % (time.time() - start_time))

    cluster_motifs(G_syn, 'motif_syn_3.txt', 'syn')
    cluster_motifs(G_exc_inh, 'motif_exc_inh_3.txt', 'exc_inh')
    cluster_motifs(G_inh_exc, 'motif_inh_exc_3.txt', 'inh_exc')
    cluster_motifs(G_neu_non, 'motif_neu_non_3.txt', 'neu_non')


    parse_networks(G_syn, "../data/processed/clusters/syn/", list_chap, list_ub, "syn", chap, ub)

    degreeDF = pd.DataFrame(columns=['OR', 'std', 'pval'])
    degreeDF.loc[len(degreeDF)], res_cat_exc_inh = reldegree(G_exc_inh, G_syn, chap, ub)
    degreeDF.loc[len(degreeDF)], res_cat_inh_exc = reldegree(G_inh_exc, G_syn, chap, ub)
    degreeDF.loc[len(degreeDF)], res_cat_neu_non = reldegree(G_neu_non, G_syn, chap, ub)
    degreeDF.index = ['Exc_Inh', 'Inh_Exc', 'Neuron_non']
    degreeDF.to_csv("../data/processed/clusters/degree.txt", header=True, index=True, sep='\t')

    categoryDF = pd.DataFrame(columns=['System', 'enrich_EI', 'pval_EI', 'enrich_IE', 'pval_IE', 'enrich_Nn', 'pval_Nn'])
    categoryDF.loc[len(categoryDF)] = np.concatenate( (np.array(res_cat_exc_inh[res_cat_exc_inh['System']=="HSP70_system"]), 
                         np.array(res_cat_inh_exc[res_cat_inh_exc['System']=="HSP70_system"][['enrich','pval']]),
                         np.array(res_cat_neu_non[res_cat_neu_non['System']=="HSP70_system"][['enrich','pval']])), axis=1 ).flatten()
    
    categoryDF.loc[len(categoryDF)] = np.concatenate( (np.array(res_cat_exc_inh[res_cat_exc_inh['System']=="HSP90_system"]), 
                         np.array(res_cat_inh_exc[res_cat_inh_exc['System']=="HSP90_system"][['enrich','pval']]),
                         np.array(res_cat_neu_non[res_cat_neu_non['System']=="HSP90_system"][['enrich','pval']])), axis=1 ).flatten()
    
    categoryDF.loc[len(categoryDF)] = np.concatenate( (np.array(res_cat_exc_inh[res_cat_exc_inh['System']=="RING_Single"]), 
                         np.array(res_cat_inh_exc[res_cat_inh_exc['System']=="RING_Single"][['enrich','pval']]),
                         np.array(res_cat_neu_non[res_cat_neu_non['System']=="RING_Single"][['enrich','pval']])), axis=1 ).flatten()
    
    categoryDF.loc[len(categoryDF)] = np.concatenate( (np.array(res_cat_exc_inh[res_cat_exc_inh['System']=="CCT/TRiC_system"]), 
                         np.array(res_cat_inh_exc[res_cat_inh_exc['System']=="CCT/TRiC_system"][['enrich','pval']]),
                         np.array(res_cat_neu_non[res_cat_neu_non['System']=="CCT/TRiC_system"][['enrich','pval']])), axis=1 ).flatten()
 
    categoryDF.loc[len(categoryDF)] = np.concatenate( (np.array(res_cat_exc_inh[res_cat_exc_inh['System']=="sHSP_system"]), 
                         np.array(res_cat_inh_exc[res_cat_inh_exc['System']=="sHSP_system"][['enrich','pval']]),
                         np.array(res_cat_neu_non[res_cat_neu_non['System']=="sHSP_system"][['enrich','pval']])), axis=1 ).flatten()
  
    categoryDF.loc[len(categoryDF)] = np.concatenate( (np.array(res_cat_exc_inh[res_cat_exc_inh['System']=="RING_Complex"]), 
                         np.array(res_cat_inh_exc[res_cat_inh_exc['System']=="RING_Complex"][['enrich','pval']]),
                         np.array(res_cat_neu_non[res_cat_neu_non['System']=="RING_Complex"][['enrich','pval']])), axis=1 ).flatten()
    
    categoryDF.loc[len(categoryDF)] = np.concatenate( (np.array(res_cat_exc_inh[res_cat_exc_inh['System']=="HECT"]), 
                         np.array(res_cat_inh_exc[res_cat_inh_exc['System']=="HECT"][['enrich','pval']]),
                         np.array(res_cat_neu_non[res_cat_neu_non['System']=="HECT"][['enrich','pval']])), axis=1 ).flatten()
    

    categoryDF.to_csv("../data/processed/clusters/degree_system.txt", header=True, index=True, sep='\t')

    resDF = parse_networks(G_exc_inh, "../data/processed/clusters/exc_inh/", list_chap, list_ub, "exc_inh", chap, ub)
    resDF = parse_networks(G_inh_exc, "../data/processed/clusters/inh_exc/", list_chap, list_ub, "inh_exc", chap, ub)
    resDF = parse_networks(G_neu_non, "../data/processed/clusters/neu_non/", list_chap, list_ub, "neu_non", chap, ub)