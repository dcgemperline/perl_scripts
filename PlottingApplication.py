import numpy as np
import pandas
from pandas import DataFrame, Series
import statsmodels.api as sm
from sklearn.linear_model import LinearRegression
import scipy, scipy.stats
import matplotlib.pyplot as plt
import matplotlib.markers
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas

df = pandas.read_csv('data.csv')

fig = Figure(figsize=(8,8))
canvas=FigureCanvas(fig)
ax=fig.add_subplot(111);

#Set Title
##ax.set_title("Linear Plot",fontsize=14)

##SetScale
ax.set_yscale('linear')
##ax.setyscale('log')

## Get rid of extraneous features for later processing in illustrator
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')

#SetAxes Labels
ax.set_xlabel("Log10([fmol])",fontsize=14)
ax.set_ylabel("Log10Average",fontsize=14)

#AxesTicks and Label Rotation sometimes useful for line graphs
##ax.minorticks_on()
##ax.axes.set_xticklabels(df.Log10Fmol, rotation=45)

## Get rid of top and right bounding axes

#DisplayGrid
##ax.grid(True,linestyle='-',color='0.75')

##Set Limits
ax.set_ylim(-5,0)
ax.set_xlim(0,5)

#Setup Line Thickness
for axis in ['left','bottom']:
    ax.spines[axis].set_linewidth(2)
ax.tick_params(width=1.5)

#Generate the Scatter Plot
ax.scatter(df.Log10Fmol, df.Log10Avg, s = 20, color = 'black')
##ax.scatter(df.Log10Fmol, df.Log10Avg, s = 20, color = 'blue')

#Generate Linear fit with numpy
m, b = np.polyfit(df.Log10Fmol, df.Log10Avg, 1)

#Add Linear fit to graph
ax.plot(df.Log10Fmol, m*df.Log10Fmol + b, '-', color ='black')

canvas.print_figure('Test.pdf',dpi= 600)
