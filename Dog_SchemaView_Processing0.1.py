import requests, sys, re, json, mne, pickle
from mne.viz import circular_layout, plot_connectivity_circle
import numpy as np
import matplotlib.pyplot as plt
import ctypes
import mne.viz.utils as utils
from functools import partial

#source files
#MASHA sequence
file = open('C:/Users/webaccus/OneDrive - Agilent Technologies/Desktop/DGAE/Plink/plink.bim', 'r')
ref1 = open('C:/Users/webaccus/OneDrive - Agilent Technologies/Desktop/DGAE/canfam3/chr1_allsnps.txt', 'r')
file.seek(0)
ref1.seek(0)
start = 0
end = 0
q = False
found = []
SNPs = []
mutations = {"rs0000": [1, 2, 3]}
genes = {"xxxx": [1, 2, 3]}
#apis
server = "http://rest.ensembl.org"
variant = "/variation/canis_lupus_familiaris?phenotypes=1;"
overlap = "/overlap/region/canis_lupus_familiaris/"
phenotype = "/phenotype/gene/canis_lupus_familiaris/"
homology = "/homology/id/ENSG00000157764?content-type=application/json"
sym2Ensg = "/xrefs/symbol/canis_lupus_familiaris/BRCA2?content-type=application/json"
LD = "/ld/canis_lupus_familiaris/pairwise/"
headers={ "Content-Type" : "application/json", "Accept" : "application/json"}
significance = "/eqtl/variant_name/canis_lupus_familiaris/"

class analyseInput():
    def __init__():
        import requests, sys, re, json

    def getFeatureOverlap(loc):
        text = requests.get(server+overlap+loc+'?feature=gene;', headers={ "Content-Type" : "application/json"})
        text = text.json()
        geneNames = ['', '']
        for i in text:
            #print(i)
            try:
                geneNames[0] = geneNames[0] + i['external_name'] + '%0d'
                geneNames[1] = geneNames[1] + i['description'] + '\n'
            except:
                geneNames = geneNames
        #print(str(geneNames))
        return geneNames

    def sliceChromosome(inputChromosome, refrenceChromosome, slices):
        print('dividing chromosome into smaller slices for analysis')
        Alist = []
        Blist = []
        Ase = []
        Bse = []
        cor = []
        sliver = slices
        endA = 0
        endB = 0
        i = 0
        for a in range(sliver):
            Alist.append([])
            for i in range(endA, int(endA+len(inputChromosome)/sliver)):
                Alist[a].append(inputChromosome[i])
            start = inputChromosome[endA].split('\t')
            start = start[3]
            end = inputChromosome[i].split('\t')
            end = end[3]
            Ase.append((int(start)+3002782, int(end)+3002782))
            endA = i
            Blist.append([])
            for i in range(endB, int(endB+len(refrenceChromosome)/sliver)):
                Blist[a].append(refrenceChromosome[i])
            start = refrenceChromosome[endB].split('\t')
            start = start[2]
            end = refrenceChromosome[i].split('\t')
            end = end[2]
            Bse.append((start, end))
            endB = i
        inputSet, refrenceSet = Alist, Blist
        inputSlices, refrenceSlices = Ase, Bse
        return inputSet, refrenceSet, inputSlices, refrenceSlices

    def alignSlices(inputSlices, refrenceSlices, inputSliceBounds, refrenceSliceBounds, slices):
        print('aligning input slices and refrence slices to optimize speed')
        location = []
        alignedSlices = []
        for i in range(slices):
            alignedSlices.append([i, []])
            for t in range(1, len(refrenceSliceBounds)):
                if int(refrenceSliceBounds[t][0]) >= int(inputSliceBounds[i][0]) and int(refrenceSliceBounds[t][0]) <= int(inputSliceBounds[i][1]) or int(refrenceSliceBounds[t][1]) <= int(inputSliceBounds[i][1]) and int(refrenceSliceBounds[t][1]) >= int(inputSliceBounds[i][0]):
                    alignedSlices[i][1].append(t)
                else:
                    continue
            if alignedSlices[i][1]:
                location.append(str(inputSliceBounds[i][0])+'-'+str(inputSliceBounds[i][1]))
        print(location)
        return alignedSlices, location

    def getBreedDominanceBySlice(inputSlices, refrenceSlices, sliceAlignments):
        print('analysing slices for breed dominance and traits')
        geneList = []
        labels = []
        thisLabel = []
        slices = {}
        for Slice, c in enumerate(sliceAlignments):
            print(Slice)
            if not c[0] or not c[1]:
                continue
            else:
                inputLine = inputSlices[c[0]][0].split('\t')
                subStart = str(inputLine[0] + ':' + inputLine[3])
                inputLine = inputSlices[c[0]][len(inputSlices[c[0]])-1].split('\t')
                loc = str(inputLine[3])
                loc = subStart + '-' + loc
                print('getting feature overlap for %s from ensembl database' % loc)
                geneList.append(analyseInput.getFeatureOverlap(loc))
                for i, inputLine in enumerate(inputSlices[c[0]]):
                    inputLine = inputLine.split('\t')
                    for s in c[1]:
                        breed = {}
                        refrenceLine = refrenceSlices[s][0].split('\t')
                        for refrenceLine in refrenceSlices[s]:
                            refrenceLine = refrenceLine.split('\t')
                            if refrenceLine[5] in breed:
                                breed[refrenceLine[5]][0] += 1
                                if str(inputLine[1]) == str(refrenceLine[0]):
                                    breed[refrenceLine[5]][1] += 1
                            else:
                                breed[refrenceLine[5]] = [1, 0]
                                if str(inputLine[1]) == str(refrenceLine[0]):
                                    breed[refrenceLine[5]][1] += 1
                        subEnd = int(refrenceLine[2])
                    try:
                        slices[Slice].append(breed)
                    except:
                        slices[Slice] = []
                        slices[Slice].append(breed)
                    inputLine = inputSlices[c[0]][len(inputSlices[c[0]])-1].split('\t')
                highest = []
                matches = {}
                for dictionary in slices[Slice]:
                    for i in dictionary:
                        if dictionary[i][1] > 0:
                            if i in matches:
                                matches[i] += dictionary[i][1]
                            else:
                                matches[i] = dictionary[i][1]
                        if len(highest) < 10:
                            highest.append((i, dictionary[i][0]))
                        else:
                            for x, cont in enumerate(highest):
                                if cont[1] < dictionary[i][0]:
                                    if i in highest != True:
                                        highest[x] = (i, dictionary[i][0])
                highest = tuple(highest)
                highest = dict((y, x) for y, x in highest)
                slices[Slice] = [matches, highest]
                highestRelative = (0, 0)
                for lst in slices[Slice][0]:
                    for subLst in slices[Slice][1]:
                        if lst == subLst:
                            if float(slices[Slice][0][lst]/slices[Slice][1][subLst]) > highestRelative[1]:
                                highestRelative = (lst, float(slices[Slice][0][lst]/slices[Slice][1][subLst]))
                slices[Slice].append(highestRelative)
        for i in slices:
            labels.append(slices[i])
        length = len(slices)
        return labels, length, geneList

class formatData():
    def __init__(self):
        import requests, json, mne
        from mne.viz import circular_layout, plot_connectivity_circle
        import numpy as np
        self.sym2ensg = "/xrefs/symbol/canis_lupus_familiaris/"
        self.orthology = "/homology/symbol/canis_lupus_familiaris/"
        self.translate = "/map/canis_lupus_familiaris/%s?"
    
    def buildCorrelationMap(listA, listB, length):
        print('finding correlations for genes in slice with other slices %s' % listA)
        cor = []
        string = 'https://string-db.org'
        interactions = '/api/tsv/network?identifiers='
        for b, y in enumerate(listB):
            if listA == y[0]:
                continue
            total = 0.0
            correlations = requests.get(string+interactions+listA+y[0])
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
        corMap = []
        labelNames = []
        labelColors = []
        for i in range(length-1):
            label = slices[i][2]
            labelNames.append(str(label[0]) + '@' + location[i])
            a = label[1]
            labelColors.append((0.75,0.3,0.1,a))
        #for i in range(length-1):
            #corMap.append(formatData.buildCorrelationMap(geneList[i][0], geneList, length-1))
        for i in range(length-1):    
            if i <5:
                corMap.append(formatData.buildCorrelationMap(geneList[i][0], geneList, length-1))
            else:
                lst = []
                for m in range(length-1):
                    lst.append(0.0)
                corMap.append(lst)
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
            ENSG = requests.get(server+self.sym2ensg+a, headers={ "Content-Type" : "application/json"})
            ENSG = ENSG.json()
            #print(ENSG['id'])
            print(a, b)
            #if type(ENSG) is list:
                #homology = requests.get(server+self.orthology+ENSG[0]["id"], headers={ "Content-Type" : "application/json"})
            homology = requests.get(server+self.orthology+a, headers={ "Content-Type" : "application/json"})
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
                                        if float(homologies[A][mapping]['perc_id']) > 60 and float(homologies[A][mapping]['perc_pos']) > 50:
                                            #print(mapping)
                                            #print('id:\t', homologies[A][mapping]['id'])
                                            #print('protein_id:\t', homologies[A][mapping]['protein_id'])
                                            #print("perc_id:\t", homologies[A][mapping]['perc_id'])
                                            #print("perc_pos:\t", homologies[A][mapping]['perc_pos'])
                                            #print('cigar_line:\t', homologies[A][mapping]['cigar_line'])
                                            #print('align_seq:\t', homologies[A][mapping]['align_seq'])
                                            homologs[a].append([mapping, homologies[A][mapping]['id'], homologies[A][mapping]['protein_id'], homologies[A][mapping]['perc_id'], homologies[A][mapping]['perc_pos'], homologies[A][mapping]['cigar_line'], homologies[A][mapping]['align_seq']])
                                    elif mapping == 'target':
                                        if float(homologies[A][mapping]['perc_id']) > 60 and float(homologies[A][mapping]['perc_pos']) > 50:
                                            homologs[a].append([mapping, homologies[A][mapping]['id'], homologies[A][mapping]['protein_id'], homologies[A][mapping]['perc_id'], homologies[A][mapping]['perc_pos'], homologies[A][mapping]['cigar_line'], homologies[A][mapping]['align_seq']])
                                            
                        break
        print(genesInSlice)
        index = input('type Gene symbol to print homologies')
        print(homologs[index])
        for i, data in enumerate(homologs[index]):
            if i > 1:
                print(data, server+self.translate % data[2])
                transcriptCoords = requests.get(server+self.translate % data[2], headers={ "Content-Type" : "application/json"}) 
                if transcriptCoords.ok:
                    transcriptCoords.json()
                    print(str(transcriptCoords))
        thisLabel = labelNames[node]
        print(sliceCorrelations)

    def prepareCorrelationGraph(self, labelNames, node, geneList, sliceCorrelations):
        print('pickling')
        genesInSlice = geneList[node]
        print(genesInSlice)
        descriptions = genesInSlice[1].split('\n')
        genesInSlice = genesInSlice[0].split('%0d')
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
        import ctypes
        import mne.viz.utils as utils
        from functools import partial

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
        instance = x+1
        print('Almost Done! Building connectivity circle and UI')
        fig, axes, indices, n_nodes, node_angles, ID, plt = plot_connectivity_circle(cor_map, label_names, n_lines=300, node_angles=node_angles, 
                             node_colors=label_colors,
                             title='Chromosome%d'% instance, show=False)
        print(fig, axes)
        callback = partial(displayGraphics.grabData, fig=fig,
                               axes=axes, indices=indices, n_nodes=n_nodes,
                               node_angles=node_angles, labelNames=labelNames, geneList=geneList, corMap=corMap)
        fig.canvas.mpl_connect('button_press_event', callback)
        utils.plt_show(True)
        #plt.show()

    def plotSubChromosome():
        print('nope')

def export(item, name):
    file = open('C:/Users/webaccus/OneDrive - Agilent Technologies/Desktop/DGAE/%s.pkl' % name, 'wb')
    pickle.dump(item, file)
    file.close()

def get(name):
    file = open('C:/Users/webaccus/OneDrive - Agilent Technologies/Desktop/DGAE/%s.pkl' % name, 'rb')
    item = pickle.load(file)
    file.close()
    return item

with open('C:/Users/webaccus/OneDrive - Agilent Technologies/Desktop/DGAE/Plink/plink.bim', 'r') as f:
    data = f.readlines()
data = data[755:]

#print chromosome cutoff lines
chromosome = 1
cutoff = []
for i, line in enumerate(data):
    line = line.split('\t')
    try:
        replaceSpace = line[5].replace(' ', '')
        data[i].replace('%s' % line[5], '%s' % replaceSpace)
    except Exception as e: print(repr(e))
    try:
        replaceComma = replaceSpace.replace(',', '')
        data[i].replace('%s' % replaceSpace, '%s' % replaceComma)
    except Exception as e: print(repr(e))
    if chromosome != line[0]:
        cutoff.append(i+1)
        chromosome = line[0]
print(cutoff)

chromosome = []
masha = []
for i in range(38):
    chromosome.append([])
    with open('C:/Users/webaccus/OneDrive - Agilent Technologies/Desktop/DGAE/canfam3/chr%d_allsnps.txt' % (i+1), 'r') as f:
        chromosome[i] = f.readlines()
    masha.append([])
    masha[i] = data[cutoff[i]:cutoff[i+1]]

for x in range(38):
    slices = 200
    
    #contents
    inputSlices, refrenceSlices, inputSliceBounds, refrenceSliceBounds = analyseInput.sliceChromosome(masha[x], chromosome[x], slices)
    
    sliceAlignments, location = analyseInput.alignSlices(inputSlices, refrenceSlices, inputSliceBounds, refrenceSliceBounds, slices)                  

    print('\n\n\n\n\n\n\n\n\n\n\n\n', x+1)
    slices, length, geneList = analyseInput.getBreedDominanceBySlice(inputSlices, refrenceSlices, sliceAlignments)    
    export(slices, 'alignedSlices')
    export(length, 'length')
    export(geneList, 'geneList')
    export(location, 'location')
    
    #print(slices, length, geneList)
    #corMap, labelNames, labelColors = formatData.mainCorrelationMap(get('geneList'), get('alignedSlices'), get('location'))
    #displayGraphics.buildConnectivityMap(get('alignedSlices'), get('length'), get('geneList'), corMap, labelNames, labelColors)    

