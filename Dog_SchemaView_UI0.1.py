import pickle, requests, json, mne
import matplotlib
from mne.viz import circular_layout, plot_connectivity_circle
import numpy as np
import matplotlib.pyplot as plt
import ctypes
import mne.viz.utils as utils
from functools import partial

def export(item, name):
    file = open('C:/Users/webaccus/OneDrive - Agilent Technologies/Desktop/DGAE/%s.pkl' % name, 'wb')
    pickle.dump(item, file)
    file.close()

def get(name):
    file = open('C:/Users/webaccus/OneDrive - Agilent Technologies/Desktop/DGAE/%s.pkl' % name, 'rb')
    item = pickle.load(file)
    file.close()
    return item

class formatData():
    def __init__(self):
        self.sym2ensg = "/xrefs/symbol/canis_lupus_familiaris/"
        self.orthology = "/homology/symbol/canis_lupus_familiaris/"
        self.translate = "/map/translation/%s/1..100000?"
        self.server = "http://rest.ensembl.org"
    
    def buildCorrelationMap(listA, listB, length):
        print('finding correlations for genes in slice with other slices %s' % listA)
        cor = []
        string = 'https://string-db.org'
        interactions = '/api/tsv/network?identifiers='
        for b, y in enumerate(listB):
            if listA == y[0]:
                continue
            total = 0.0
            correlations = requests.get(string+interactions+listA+y[0]+'&species=9612')
            correlations = correlations.text
            correlations = correlations.split('\n')
            for c, line in enumerate(correlations):
                line = line.split('\t')
                if len(line) > 5:
                    if line[5] != 'score':
                        #print(line[2],line[3],line[5])
                        if line[2] in listA or line[3] in listA:
                            total += float(line[5])
            cor.append(float(total)/c)
        if len(cor) < length:
            for m in range(length-len(cor)):
                lst = []
                for i in range(length):
                    lst.append(0.0)
                cor.append(lst)
        #print(len(cor))
        return cor

    def mainCorrelationMap(geneList, slices, location):
        length = len(slices)
        corMap = []
        labelNames = []
        labelColors = []
        for i in range(length-1):
            label = slices[i][2]
            #labelNames.append(str(label[0]) + '@' + location[i])
            labelNames.append(str(label[0])+str(i))
            a = label[1]
            labelColors.append((0.75,0.3,0.1,a))
        for i in range(length-1):
            corMap.append(formatData.buildCorrelationMap(geneList[i][0], geneList, length-1))
        #for i in range(length-1):    
        #    if i <5:
        #        corMap.append(formatData.buildCorrelationMap(geneList[i][0], geneList, length-1))
        #    else:
        #        lst = []
        #        for m in range(length-1):
        #            lst.append(0.0)
        #        corMap.append(lst)
        print(corMap)
        return corMap, labelNames, labelColors

    def prepareHomologyGraph(self, labelNames, node, geneList, sliceCorrelations):
        genesInSlice = geneList[node]
        descriptions = genesInSlice[1].split('\n')
        genesInSlice = genesInSlice[0].split('%0d')
        homologs = {}
        filterOrthologs = input('type True to filter orthologs or False to filter paralogs')
        if json.loads(filterOrthologs.lower()) == True:
            typeFilter = 'ENSEMBL_ORTHOLOGUES'
        else:
            typeFilter = 'ENSEMBL_PARALOGUES'
        print(typeFilter)
        for i, a in enumerate(genesInSlice):
            b = descriptions[i]
            ENSG = requests.get(self.server+self.sym2ensg+a, headers={ "Content-Type" : "application/json"})
            ENSG = ENSG.json()
            #print(ENSG['id'])
            print(a, b)
            #if type(ENSG) is list:
                #homology = requests.get(server+self.orthology+ENSG[0]["id"], headers={ "Content-Type" : "application/json"})
            homology = requests.get(self.server+self.orthology+a, headers={ "Content-Type" : "application/json"})
            if homology.ok:
                homology = homology.json()
                #print(str(homology))
                homologs[a] = []
                homologs[a].append(['id', 'protein_id', 'perc_id', 'per_pos', 'cigar_line', 'align_seq'])
                for x in homology:
                    if x != 'error':
                        homologies = homology['data'][0]['homologies']
                        for A, data in enumerate(homologies):
                            #print(data)
                            #print('\n\nmethod_link_type:\t', data['method_link_type'])
                            #print(typeFilter)
                            if data['method_link_type'] != typeFilter:
                                print('\n\nmethod_link_type:\t', data['method_link_type'])
                                for mapping in homologies[A]:
                                    #print(mapping)
                                    if mapping == 'source':
                                        if float(homologies[A][mapping]['perc_id']) > 40 and float(homologies[A][mapping]['perc_pos']) > 50:
                                            print(mapping)
                                            print('id:\t', homologies[A][mapping]['id'])
                                            print('protein_id:\t', homologies[A][mapping]['protein_id'])
                                            print("perc_id:\t", homologies[A][mapping]['perc_id'])
                                            print("perc_pos:\t", homologies[A][mapping]['perc_pos'])
                                            print('cigar_line:\t', homologies[A][mapping]['cigar_line'])
                                            print('align_seq:\t', homologies[A][mapping]['align_seq'])
                                            homologs[a].append([mapping, homologies[A][mapping]['id'], homologies[A][mapping]['protein_id'], homologies[A][mapping]['perc_id'], homologies[A][mapping]['perc_pos'], homologies[A][mapping]['cigar_line'], homologies[A][mapping]['align_seq']])
                                    elif mapping == 'target':
                                        if float(homologies[A][mapping]['perc_id']) > 40 and float(homologies[A][mapping]['perc_pos']) > 50:
                                            homologs[a].append([mapping, homologies[A][mapping]['id'], homologies[A][mapping]['protein_id'], homologies[A][mapping]['perc_id'], homologies[A][mapping]['perc_pos'], homologies[A][mapping]['cigar_line'], homologies[A][mapping]['align_seq']])
                                            
                        break
        print(genesInSlice)
        index = input('type Gene symbol to print homologies')
        print(homologs[index])
        for i, data in enumerate(homologs[index]):
            if i > 1:
                print(data, self.server+self.translate % data[2])
                transcriptCoords = requests.get(self.server+self.translate % data[2], headers={ "Content-Type" : "application/json"}) 
                if transcriptCoords.ok:
                    transcriptCoords = transcriptCoords.json()
                    print(str(transcriptCoords))
        thisLabel = labelNames[node]
        print(sliceCorrelations)

    def prepareCorrelationGraph(self, labelNames, node, geneList, sliceCorrelations):
        print('pickling')
        genesInSlice = geneList[node]
        print(genesInSlice)
        string = 'https://string-db.org'
        interactions = '/api/tsv/network?identifiers='
        for i, percent in enumerate(sliceCorrelations):
            print(percent)
            if percent > 0:
                correlations = requests.get(string+interactions+geneList[i][0]+genesInSlice[0]+'&species=9612')
                correlations = correlations.text
                print(correlations)
            #correlations = correlations.split('\n')
        descriptions = genesInSlice[1].split('\n')
        genesInSlice = genesInSlice[0].split('%0d')
        for i, gene in enumerate(genesInSlice):
            print(gene)
            print(descriptions[i])
        correlations = {'genesInSlice':  open('C:/Users/webaccus/OneDrive - Agilent Technologies/Desktop/DGAE/genesInSlice.pkl', 'wb'), 'geneList':  open('C:/Users/webaccus/OneDrive - Agilent Technologies/Desktop/DGAE/geneList.pkl', 'wb'), 'sliceCorrelations': open('C:/Users/webaccus/OneDrive - Agilent Technologies/Desktop/DGAE/sliceCorrelations.pkl', 'wb')}
        pickle.dump(genesInSlice, correlations['genesInSlice'])
        correlations['genesInSlice'].close()
        pickle.dump(geneList, correlations['geneList'])
        correlations['geneList'].close()
        pickle.dump(sliceCorrelations, correlations['sliceCorrelations'])
        correlations['sliceCorrelations'].close()

class displayGraphics():
    def __init__():
        import matplotlib.pyplot as plt

    def grabData(event, fig=None, axes=None, indices=None,
                                         n_nodes=0, node_angles=None,
                                         ylim=[9, 10], labelNames=None, geneList=None, corMap=None):
        if event.inaxes != axes:
            return

        #event buttons: 1 = left click, 2 = double click, 3 = right click
        if event.button == 3:  # right click
            # click must be near node radius
            if not ylim[0] <= event.ydata <= ylim[1]:
                return
            node_angles = node_angles % (np.pi * 2)
            node = np.argmin(np.abs(event.xdata - node_angles))
            return formatData().prepareHomologyGraph(labelNames, node, geneList, corMap[node])
        elif event.button == 1:  # left click
            # click must be near node radius
            if not ylim[0] <= event.ydata <= ylim[1]:
                return
            node_angles = node_angles % (np.pi * 2)
            node = np.argmin(np.abs(event.xdata - node_angles))
            return formatData().prepareCorrelationGraph(labelNames, node, geneList, corMap[node])

    def buildConnectivityMap(slices, length, geneList, corMap, labelNames, labelColors):
        node_order = list()
        node_order.extend(labelNames)
        node_angles = circular_layout(labelNames, node_order, start_pos=90,
                                  group_boundaries=None)
        label_names = np.array(labelNames)
        label_colors = np.array(labelColors)
        cor_map = np.array(corMap)
        #instance = x+1
        print('Almost Done! Building connectivity circle and UI')
        fig, axes, indices, n_nodes, node_angles, ID, plt = plot_connectivity_circle(cor_map, label_names, n_lines=300, node_angles=node_angles, 
                             node_colors=label_colors,
                             title='Chromosome', show=False)
        print(fig, axes)
        callback = partial(displayGraphics.grabData, fig=fig,
                               axes=axes, indices=indices, n_nodes=n_nodes,
                               node_angles=node_angles, labelNames=labelNames, geneList=geneList, corMap=corMap)
        fig.canvas.mpl_connect('button_press_event', callback)
        utils.plt_show(True)
        #plt.show()

    def plotSubChromosome():
        print('nope')
        
corMap, labelNames, labelColors = formatData.mainCorrelationMap(get('geneList'), get('alignedSlices'), get('location'))
displayGraphics.buildConnectivityMap(get('alignedSlices'), get('length'), get('geneList'), corMap, labelNames, labelColors)    
