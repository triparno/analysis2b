#Imported Modules
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator, MultipleLocator,FixedLocator
import matplotlib
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import AutoMinorLocator
from matplotlib import rc
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['pgf.texsystem'] = 'pdflatex'
matplotlib.rcParams['text.usetex']=True
matplotlib.rcParams['pgf.texsystem'] = 'pdflatex'
matplotlib.rc('xtick', labelsize=11)
matplotlib.rc('ytick', labelsize=11)
matplotlib.rcParams['xtick.major.pad'] = 3
matplotlib.rcParams['ytick.major.pad'] = 3
matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{amssymb}"]
matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{bm}"]
plt.gca().xaxis.set_major_locator(MaxNLocator(prune='both'))
plt.gca().yaxis.set_major_locator(MaxNLocator(prune='both'))
matplotlib.rcParams['xtick.direction'] = 'in'
matplotlib.rcParams['ytick.direction'] = 'in'
matplotlib.rcParams['xtick.major.size'] = 6
matplotlib.rcParams['ytick.major.size'] = 6
matplotlib.rcParams['xtick.minor.size'] = 3
matplotlib.rcParams['ytick.minor.size'] = 3
matplotlib.rcParams['xtick.major.width'] = 0.8
matplotlib.rcParams['ytick.major.width'] = 0.8
matplotlib.rcParams['xtick.minor.width'] = 0.5
matplotlib.rcParams['ytick.minor.width'] = 0.5
labpad = 10

matplotlib.rcParams['lines.solid_joinstyle'] = 'round'
sig2 = '#ffee66'
cyan = '#669dff'
sig2 = '#ffee66'
sig2 = '#dab600'
orange = '#ff8000'
'''
sig2 = '#ffee66'
sig2 = '#dab600'
sig1 = '#99dd99'
sig2 = '#dab600
'''
sig1 = '#99dd99'
red = '#92130f'
afblue = '#5D8AA8'
blue = '#0f4d92'
green = '#0f9254'
chrome = '#ffe866'
violet = '#B53389'


grade_red_5 = ['#bd0026', '#f03b20', '#fd8d3c', '#fecc5c', '#ffffb2']
grade_blu_5 = ['#253494', '#2c7fb8', '#41b6c4', '#a1dab4', '#ffffcc']
grade_vlt_5 = ['#7a0177', '#c51b8a', '#f768a1', '#fbb4b9', '#feebe2']
