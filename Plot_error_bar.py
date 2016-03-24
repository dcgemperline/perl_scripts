import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42

import numpy as np
import matplotlib.pyplot as plt
import pandas
from scipy import stats


font = {'family' : 'Arial',
        'weight' : 'bold',
        'size' : 6.5}
mpl.rc('font',**font)
def getdf(file):
    df = pandas.read_csv(file)
    return df

def setupErrorBarPlot(ax0, x, y, yerror = None):
    #ax0.set_title(test)
    ##SetScale
    ax0.set_yscale('linear')
    ##ax.setyscale('log')

    ## Get rid of extraneous features for later processing in illustrator
    ax0.spines['right'].set_visible(False)
    ax0.spines['bottom'].set_visible(False)
    ax0.xaxis.set_ticks_position('top')
    ax0.yaxis.set_ticks_position('left')

    #Set X Axis Top
    ax0.xaxis.tick_top()

    #SetAxes Labels
    ax0.set_xlabel("Log10(Average[fmol])",fontsize=7)
    ax0.set_ylabel("Log10(Average(dNSAF))",fontsize=7)

    #AxesTicks and Label Rotation sometimes useful for line graphs
    ##ax.minorticks_on()
    ##ax.axes.set_xticklabels(df.Log10Fmol, rotation=45)

    ## Get rid of top and right bounding axes

    #DisplayGrid
    ##ax.grid(True,linestyle='-',color='0.75')

    ##Set Limits - U are goign to change this tomorrow
    ## this will also cause the ticks to be uneven, so you will have to find
    ##Some what to specify this
    
    ax0.set_xlim(1,5)
    ax0.set_ylim(-5,0)

    
    #Setup Line Thickness
    for axis in ['left','bottom', 'top','right']:
        ax0.spines[axis].set_linewidth(1)
    ax0.tick_params(width=1)


    #PlotErrorBarPlot
    (_, caps, _) = ax0.errorbar(x, y, yerr=yerror, fmt='o', ecolor='black',
                                markersize='3.5', markeredgewidth=1.0,
                                markerfacecolor = 'none', markeredgecolor ='black',
                                capthick=0.5, elinewidth=0.5)
    for cap in caps:
        #Sets Errorbar Cap Width in Pixls
        cap.set_markeredgewidth(0.5)
        #Sets Cap Size
        cap.set_markersize(2)
    
   
    ##Generate Linear fit using scipy
    slope, intercept, r_value, p_value, std_error = stats.linregress(x, y)
    
    ##Get R2
    r_value = r_value**2

    if intercept < 0:
        equationstring = "y = {0:.2f}x - {1:.2f}".format(slope, intercept*-1)
    else:
        equationstring = "y = {0:.2f}x + {1:.2f}".format(slope, intercept)

    r_valuestring = "{0:.3f}".format(r_value)
    r_equation = r'R$^2$ ' + "=" + r_valuestring

    #Add Linear fit to graph
    ax0.plot(x, slope*x + intercept, '--', color ='black', linewidth=1.5)

    ##Adjust Ticks
    max_xticks = 5
    xloc = plt.MaxNLocator(max_xticks)
    ax0.xaxis.set_major_locator(xloc)
    ax0.tick_params(direction ='inout')
    ax0.tick_params(length = 4.8)



    ## Add Equation
    #ax0.text(2, -1, 'test')
    TextToAdd = equationstring + '\n' + r_equation
    ax0.text(2, -1, TextToAdd,  fontsize=8)

def OutputErrorBarPlot(inputfile, outputfile):
    df = getdf(inputfile)
    fig, ax0 = plt.subplots()
    fig.set_clip_on=False
    fig.set_size_inches(2.45,2.31)
    #fig.set_size_inches(2.27,2.083)
    setupErrorBarPlot(ax0, df.Log10Fmol, df.Log10Avg, df.Error)
    plt.savefig(outputfile)

def setupScatterPlot(ax0, x, y, yerror = None):
    #ax0.set_title(test)
    ##SetScale
    ax0.set_yscale('linear')
    ##ax.setyscale('log')

    ## Get rid of extraneous features for later processing in illustrator
    ax0.spines['right'].set_visible(False)
    ax0.spines['bottom'].set_visible(False)
    ax0.xaxis.set_ticks_position('top')
    ax0.yaxis.set_ticks_position('left')

    #Set X Axis Top
    ax0.xaxis.tick_top()

    #SetAxes Labels
    ax0.set_xlabel("Log10(Average[fmol])",fontsize=7, fontweight = 'bold')
    ax0.set_ylabel("Log10(Average(dNSAF))",fontsize=7, fontweight = 'bold')

    ##Special Move for Putting Axis at 0,0

    ax0.spines['left'].set_position('zero')

    #AxesTicks and Label Rotation sometimes useful for line graphs
    ##ax.minorticks_on()
    ##ax.axes.set_xticklabels(df.Log10Fmol, rotation=45)

    ## Get rid of top and right bounding axes

    #DisplayGrid
    ##ax.grid(True,linestyle='-',color='0.75')

    ##Set Limits - U are goign to change this tomorrow
    ## this will also cause the ticks to be uneven, so you will have to find
    ##Some what to specify this
    
    ax0.set_xlim(-1,5)
    ax0.set_ylim(-5,0)

    
    #Setup Line Thickness
    for axis in ['left','bottom', 'top','right']:
        ax0.spines[axis].set_linewidth(1)
    ax0.tick_params(width=1)


    ax0.scatter(x, y, s=10, marker='o', facecolor ='none', edgecolor = 'black', lw= 1)
    
   
    ##Generate Linear fit using scipy
    slope, intercept, r_value, p_value, std_error = stats.linregress(x, y)
    
    #Square the Result to get R-squared
    r_value = r_value**2

    if intercept < 0:
        equationstring = "y = {0:.2f}x - {1:.2f}".format(slope, intercept*-1)
    else:
        equationstring = "y = {0:.2f}x + {1:.2f}".format(slope, intercept)

    r_valuestring = "{0:.3f}".format(r_value)
    r_equation = r'R$^2$ ' + "=" + r_valuestring

    #Add Linear fit to graph
    ax0.plot(x, slope*x + intercept, '--', color ='black', linewidth=1.5)

    ##Adjust Ticks
    max_xticks = 7
    xloc = plt.MaxNLocator(max_xticks)
    ax0.xaxis.set_major_locator(xloc)
    ax0.tick_params(direction ='inout')
    ax0.tick_params(length = 4.8)



    ## Add Equation
    #ax0.text(2, -1, 'test')
    TextToAdd = equationstring + '\n' + r_equation
    ax0.text(1, -1, TextToAdd, fontsize=8)

def OutputScatterPlot(inputfile, outputfile):
    df = getdf(inputfile)
    fig, ax0 = plt.subplots()
    fig.set_clip_on=False
    fig.set_size_inches(3.6,2.31)
    setupScatterPlot(ax0, df.Log10Fmol, df.Log10Avg)
    plt.savefig(outputfile)


OutputErrorBarPlot("data_embryo_1pep.csv", "data_embryo_1pep.pdf")
OutputErrorBarPlot("data_egg_1pep.csv", "data_egg_1pep.pdf")
OutputErrorBarPlot("data_embryo_2pep.csv", "data_embryo_2pep.pdf")
OutputErrorBarPlot("data_egg_2pep.csv", "data_egg_2pep.pdf")

OutputScatterPlot("data_embryo_1pep_scatter.csv", "data_embryo_1pep_scatter.pdf")
OutputScatterPlot("data_egg_1pep_scatter.csv", "data_egg_1pep_scatter.pdf")
OutputScatterPlot("data_embryo_2pep_scatter.csv", "data_embryo_2pep_scatter.pdf")
OutputScatterPlot("data_egg_2pep_scatter.csv", "data_egg_2pep_scatter.pdf")





