import requests, sys, re, json, mne, pickle, os, ctypes, threading
from mne.viz import circular_layout, plot_connectivity_circle
import numpy as np
import matplotlib.pyplot as plt
import PySimpleGUIQt as sg
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

def export(item, name):
    file = open('C:/Users/webaccus/OneDrive - Agilent Technologies/Desktop/DGAE/%s.pkl' % name, 'wb')
    pickle.dump(item, file)
    file.close()

def get(name):
    file = open('C:/Users/webaccus/OneDrive - Agilent Technologies/Desktop/DGAE/%s.pkl' % name, 'rb')
    item = pickle.load(file)
    file.close()
    return item

#simpleGUI user interface
layout = [
    [sg.Text('select a .bim file containing the sequenced dna of your dog')],
    [sg.In(), sg.FileBrowse()],
    [sg.Text('Type the name of your dog', size=(15, 1)), sg.InputText()],
    [sg.Submit(), sg.Cancel()]
]
window = sg.Window('Simple data entry window', layout)
event, values = window.read()
window.close()
print(event, values)
fname = values[0]
dogName = values[1]

with open(fname, 'r') as f:
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

layout = [  [sg.Text('What you print will display below:')],
            [sg.Multiline(key='-OUTPUT-', size=(50, 10), autoscroll=True)]  ]

window = sg.Window('Window Title', layout, finalize=True)

#dogName = input('type the name of your dog to create a new project file')
progress = '##'

for x in range(38):
    slices = 200
    progress += '##'

    window['-OUTPUT-'].update('analyzing chromosome %d \nprogress: %s' % (x+1, progress))
    window.Refresh()
    
    #contents
    inputSlices, refrenceSlices, inputSliceBounds, refrenceSliceBounds = analyseInput.sliceChromosome(masha[x], chromosome[x], slices)
    
    sliceAlignments, location = analyseInput.alignSlices(inputSlices, refrenceSlices, inputSliceBounds, refrenceSliceBounds, slices)                  

    print('\n\n\n\n\n\n\n\n\n\n\n\n', x+1)
    slices, length, geneList = analyseInput.getBreedDominanceBySlice(inputSlices, refrenceSlices, sliceAlignments)    

    newpath = r'C:/Users/webaccus/OneDrive - Agilent Technologies/Desktop/DGAE/%s/%d' % (dogName, x+1) 
    if not os.path.exists(newpath):
        os.makedirs(newpath)
    
    export(slices, '%s/%d/alignedSlices' % (dogName, x+1))
    export(length, '%s/%d/length' % (dogName, x+1))
    export(geneList, '%s/%d/geneList' % (dogName, x+1))
    export(location, '%s/%d/location' % (dogName, x+1))
    
    #print(slices, length, geneList)
    #corMap, labelNames, labelColors = formatData.mainCorrelationMap(get('geneList'), get('alignedSlices'), get('location'))
    #displayGraphics.buildConnectivityMap(get('alignedSlices'), get('length'), get('geneList'), corMap, labelNames, labelColors)    

