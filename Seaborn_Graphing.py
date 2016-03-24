import seaborn as sns
import matplotlib as mpl
import numpy as np
## Allows for easy text manipulation in illustrator by changing this
mpl.rcParams['pdf.fonttype'] = 42

import matplotlib.pyplot as plt
import pandas as pd

font = {'family' : 'Arial',
        'weight' : 'bold',
        'size' : 14}
mpl.rc('font',**font)


## Setup figure Style
sns.set(context="paper", font_scale=2, style = "white")
sns.despine()

def GenerateFigure(dataTable, curveTable, outputFileName, labelstyle, plotstyle, shouldLabel):
    ##Read in Data
    df = pd.read_table(dataTable);
    ##Read in Curve
    df2 = pd.read_table(curveTable)
    ## NOTES
    ## replace NaN in subcomplex labels with unlabeled so that it plots correctly, and then setup color
    ## hue so that it turns grey
    ## 
    df = df.replace({'SubComplexAnnotation': { np.NaN : 'Unlabeled'}})


    ## Think about doing a 1D density plot on the axes instead of a contour plot

    # Subset the Volcano Plot by ProteasomeSubComplex

    RP = df.loc[df['SubComplexAnnotation'] == 'RP']
    CP = df.loc[df['SubComplexAnnotation'] == 'CP']
    PIPS_CP = df.loc[df['SubComplexAnnotation'] == 'PIP_CP']
    PIPS_RP = df.loc[df['SubComplexAnnotation'] == 'PIP_RP']
    PIP = df.loc[df['SubComplexAnnotation'] == 'PIP']

    ## Setup XKCD Colors
    colors = ["cool grey", "windows blue", "pale red", "medium green"]
    ## Parse Colors
    colors = sns.xkcd_palette(colors);

    ##SetupCusotmPalleteDictionary
    if labelstyle=="AllPips":
        custompalette=dict(Unlabeled=colors[0], RP=colors[1], CP=colors[2], PIP=colors[3], PIP_CP=colors[3], PIP_RP=colors[3])
    if labelstyle=="CPPips":
        custompalette=dict(Unlabeled=colors[0], RP=colors[1], CP=colors[2], PIP=colors[3], PIP_CP=colors[3], PIP_RP=colors[0])
    if labelstyle=="RPPips":
        custompalette=dict(Unlabeled=colors[0], RP=colors[1], CP=colors[2], PIP=colors[3], PIP_CP=colors[0], PIP_RP=colors[3])
    if labelstyle=="ECM29":
        custompalette=dict(Unlabeled=colors[0], RP=colors[1], CP=colors[2], PIP=colors[3], PIP_CP=colors[0], PIP_RP=colors[0])
    
    



    # Draw the two density plots
    facetgrid = sns.lmplot('Difference','-Log(P-value)', data=df, hue='SubComplexAnnotation', fit_reg=False, palette = custompalette, scatter_kws={'s' : 100}, legend=False, size = 10, aspect = 1)

    ## Force draw of points from facetgrid by axcessing facetgrids ax, and then setting collection zorder
    plt.setp(facetgrid.ax.collections, zorder =100)

    red = sns.color_palette("Reds")[-2]
    blue = sns.color_palette("Blues")[-2]

    sns.kdeplot(CP['Difference'], CP['-Log(P-value)'], cmap="Reds", bw=0.5, shade=True, shade_lowest=False, alpha = 0.5, ax=facetgrid.ax)
    if plotstyle=="kde":
        sns.kdeplot(CP['Difference'], label='CP Density', color =red, bw=0.5, shade=True, shade_lowest=False, alpha = 0.5)
    #ax.collections[0].set_alpha(0)
    sns.kdeplot(RP['Difference'], RP['-Log(P-value)'], cmap="Blues", bw=0.5, shade=True, shade_lowest=False, alpha = 0.5, ax=facetgrid.ax)
    if plotstyle=="kde":
        sns.kdeplot(RP['Difference'], label='RP Density', color =blue, bw=0.5, shade=True, shade_lowest=False, alpha = 0.5)
    #ax.collections[10].set_alpha(0)

    
    plt.plot(df2.iloc[:,0],df2.iloc[:,1], color = "black", label = "SAM Significance")
    
    ## Setup Legend
    facetgrid.fig.get_axes()[0].legend(loc='upper left', title = 'Key')
    #labels = {'Unlabeled','RP','CP','PIPS_CP', 'PIPS_RP', 'PIPs'}

    #facetgrid.add_legend(labels,title='Key')
    #facetgrid.legend(loc='upper left')

    ## Setup Axes
    myaxes = facetgrid.axes
    myaxes[0,0].set_xlim(-8,13)
    myaxes[0,0].set_ylim(0,7.5)




    # Add labels to the plot
    #red = sns.color_palette("Reds")[-2]
    #blue = sns.color_palette("Blues")[-2]
    #ax.text(2.5, 8.2, "virginica", size=16, color=blue)
    #ax.text(3.8, 4.5, "setosa", size=16, color=red)

    #sns.plt.show();
    if shouldLabel:
        if labelstyle=="CPPips":
            for index, row in PIPS_CP.iterrows():
                xposition=row['Difference']
                yposition=row['-Log(P-value)']
                labeltoAdd=row['ProteasomeAnnotation']
                offset=0.1
                facetgrid.ax.text(xposition+offset,yposition+offset, labeltoAdd, size=16, color = colors[3])
            for index, row in PIP.iterrows():
                xposition=row['Difference']
                yposition=row['-Log(P-value)']
                labeltoAdd=row['ProteasomeAnnotation']
                offset=0.1
                facetgrid.ax.text(xposition+offset,yposition+offset, labeltoAdd, size=16, color = colors[3])
        if labelstyle=="RPPips":
            for index, row in PIPS_RP.iterrows():
                xposition=row['Difference']
                yposition=row['-Log(P-value)']
                labeltoAdd=row['ProteasomeAnnotation']
                offset=0.1
                facetgrid.ax.text(xposition+offset,yposition+offset, labeltoAdd, size=16, color = colors[3])
            for index, row in PIP.iterrows():
                xposition=row['Difference']
                yposition=row['-Log(P-value)']
                labeltoAdd=row['ProteasomeAnnotation']
                offset=0.1
                facetgrid.ax.text(xposition+offset,yposition+offset, labeltoAdd, size=16, color = colors[3])
        if labelstyle=="AllPips":
            for index, row in PIPS_CP.iterrows():
                xposition=row['Difference']
                yposition=row['-Log(P-value)']
                labeltoAdd=row['ProteasomeAnnotation']
                offset=0.1
                facetgrid.ax.text(xposition+offset,yposition+offset, labeltoAdd, size=16, color = colors[3])
            for index, row in PIPS_RP.iterrows():
                xposition=row['Difference']
                yposition=row['-Log(P-value)']
                labeltoAdd=row['ProteasomeAnnotation']
                offset=0.1
                facetgrid.ax.text(xposition+offset,yposition+offset, labeltoAdd, size=16, color = colors[3])
            for index, row in PIP.iterrows():
                xposition=row['Difference']
                yposition=row['-Log(P-value)']
                labeltoAdd=row['ProteasomeAnnotation']
                offset=0.1
                facetgrid.ax.text(xposition+offset,yposition+offset, labeltoAdd, size=16, color = colors[3])
        if labelstyle=="ECM29":
            for index, row in PIP.iterrows():
                xposition=row['Difference']
                yposition=row['-Log(P-value)']
                labeltoAdd=row['ProteasomeAnnotation']
                offset=0.1
                facetgrid.ax.text(xposition+offset,yposition+offset, labeltoAdd, size=16, color = colors[3])


    facetgrid.savefig(outputFileName + ".pdf");
    facetgrid.savefig(outputFileName + ".png");


#LabelList = ["PAG1_Minus_vs_Col0_Minus",
#    "PAG1_Minus_vs_PAG1_Plus",
#    "PAG1_Plus_vs_Col0_Plus",
#    "RPT4a_Minus_vs_Col0_Minus",
#    "RPT4a_Minus_vs_RPT4b_Minus",
#    "RPT4a_Plus_vs_Col0_Plus",
#    "RPT4a_Plus_vs_RPT4b_Plus",
#    "RPT4b_Minus_vs_Col0_Minus",
#    "RPT4b_Minus_vs_RPT4b_Plus",
#    "RPT4b_Plus_vs_Col0_Plus",
#    "RPTa_Minus_vs_RPT4a_Plus"]


#TypeList = ["CPPips",
#            "CPPips",
#            "CPPips",
#            "RPPips",
#            "RPPips",
#            "RPPips",
#            "RPPips",
#            "RPPips",
#            "RPPips"] Length here is wrong, make same dimensions as LabelList



GenerateFigure("PAG1_Minus_Col0_Minus.txt", "PAG1_Minus_Col0_Minus_Curve.txt", "PAG1_Minus_Col0_Minus","CPPips","kde", True)
GenerateFigure("PAG1_Plus_Col0_Plus.txt", "PAG1_Plus_Col0_Plus_Curve.txt", "PAG1_Plus_Col0_Plus","CPPips","kde", True)
GenerateFigure("RPT4a_Minus_Col0_Minus.txt", "RPT4a_Minus_Col0_Minus_Curve.txt", "RPT4a_Minus_Col0_Minus","RPPips","kde", True)
GenerateFigure("RPT4a_Plus_Col0_Plus.txt", "RPT4a_Plus_Col0_Plus_Curve.txt", "RPT4a_Plus_Col0_Plus","RPPips","kde", True)
#GenerateFigure("RPT4b_Minus_Col0_Minus.txt", "RPT4b_Minus_Col0_Minus_curve.txt", "RPT4b_Minus_Col0_Minus","RPPips","kde", True)
#GenerateFigure("RPT4b_Plus_Col0_Plus.txt", "RPT4b_Plus_Col0_Plus_curve.txt", "RPT4b_Plus_Col0_Plus","RPPips","kde", True)

