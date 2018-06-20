#!/usr/bin/python3
import pandas as pd
import matplotlib.pyplot as plt
import re
import os
import sys
from matplotlib.colors import ListedColormap
from os import walk
import matplotlib.pyplot as plt
plt.switch_backend('agg')

def plot_bar(x, tax_level):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ylabel = list(Species.transpose().index)
    barplot = x.transpose().plot(kind='barh', stacked=True, cmap='gist_ncar', edgecolor='black', linewidth=1, width = 0.8, legend=False, figsize = (20,(len(x.transpose())*2)), ax = ax)

    barplot.axes.set_title("Relative abundance at the {} level".format(tax_level),fontsize=22)
    barplot.set_xlabel("Percentage",fontsize=20)
    barplot.set_ylabel("Samples",fontsize=20)
    barplot.spines['right'].set_visible(False)
    barplot.spines['top'].set_visible(False)
    barplot.spines['left'].set_visible(False)
    barplot.spines['bottom'].set_visible(False)
    ax.set_yticklabels(ylabel)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    plt.xticks(fontsize = 16)
    plt.yticks(fontsize=14)
    plt.subplots_adjust(left=0, bottom=0.5, right=1, top=1, wspace=0, hspace=0)
    barplot.figure.savefig('Barplot_{}.pdf'.format(tax_level), bbox_inches='tight')

    legend = barplot.legend(loc=9, bbox_to_anchor=(1.5, 1), prop={'size': 15})

    def export_legend(legend, tax_level):
        fig  = legend.figure
        fig.canvas.draw()
        bbox  = legend.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
        fig.savefig('legend_{}.pdf'.format(tax_level), bbox_inches=bbox)
    export_legend(legend, tax_level)
    ax.legend().set_visible(False)




path = sys.argv[1]

f=[]

#Read all tables inside a folder
for (dirpath, dirnames, filenames) in walk(path):
    f.extend(filenames)
    break




#Join all tables in a single matrix
dfs = [pd.read_csv(path +filename, index_col = 0, sep = '\t').rename(index=str, columns={"Metaphlan2_Analysis": re.sub(".{1,}_","",re.sub("_abundance.{1,}","", filename))}) for filename in f]
Df_final = dfs[0].join(dfs[1:], how = 'outer').fillna(0.0).reset_index()
Df_final['#SampleID'] = Df_final['#SampleID'].replace(to_replace = '.{1,}\|', value = '', regex = True)
Df_final = Df_final.set_index('#SampleID')

#Split tables in 7 tables based on tax_level
Kingdom = Df_final[Df_final.index.str.match('k__')].reset_index().replace(to_replace = 'k__', value = '', regex = True).rename(columns={'#SampleID': 'Kingdom'}).set_index('Kingdom')
Phyllum = Df_final[Df_final.index.str.match('p__')].reset_index().replace(to_replace = 'p__', value = '', regex = True).rename(columns={'#SampleID': 'Phyllum'}).set_index('Phyllum')
Class = Df_final[Df_final.index.str.match('c__')].reset_index().replace(to_replace = 'c__', value = '', regex = True).rename(columns={'#SampleID': 'Class'}).set_index('Class')
Order = Df_final[Df_final.index.str.match('o__')].reset_index().replace(to_replace = 'o__', value = '', regex = True).rename(columns={'#SampleID': 'Order'}).set_index('Order')
Family = Df_final[Df_final.index.str.match('f__')].reset_index().replace(to_replace = 'f__', value = '', regex = True).rename(columns={'#SampleID': 'Family'}).set_index('Family')
Genus = Df_final[Df_final.index.str.match('g__')].reset_index().replace(to_replace = 'g__', value = '', regex = True).rename(columns={'#SampleID': 'Genus'}).set_index('Genus')
Species = Df_final[Df_final.index.str.match('s__')].reset_index().replace(to_replace = 's__', value = '', regex = True).replace(to_replace = '_', value = ' ', regex = True).rename(columns={'#SampleID': 'Species'}).set_index('Species')

#Plot tax_levels
plot_bar(Kingdom, Kingdom.index.name)
plot_bar(Phyllum, Phyllum.index.name)
plot_bar(Class, Class.index.name)
plot_bar(Order, Order.index.name)
plot_bar(Family, Family.index.name)
plot_bar(Genus, Genus.index.name)
plot_bar(Species, Species.index.name)
































