'''
Script to plot PWMs and first order Matrices
'''

import matplotlib as mpl
from matplotlib.text import TextPath
from matplotlib.patches import PathPatch
from matplotlib.font_manager import FontProperties
import itertools
import operator
import numpy as np
import matplotlib.pyplot as plt
#define global parameters for plotting and the letters  as global variables
fp = FontProperties(family="Arial", weight="bold")
globscale = 1.35
LETTERS = { "T" : TextPath((-0.305, 0), "T", size=1, prop=fp),
            "G" : TextPath((-0.384, 0), "G", size=1, prop=fp),
            "A" : TextPath((-0.35, 0), "A", size=1, prop=fp),
            "C" : TextPath((-0.366, 0), "C", size=1, prop=fp) }
COLOR_SCHEME = {'G': 'orange',
                'A': 'forestgreen',
                'C': 'mediumblue',
                'T': 'crimson'}
COLOR_SCHEME_neg = {'G': 'khaki',
                'A': 'lightgreen',
                'C': 'lightblue',
                'T': 'salmon'}
def letterAt(letter, x, y, yscale=1.0, ax=None,negativ=False):
    '''
    Code taken from the stackoverflow question "sequence logos in matplotlib: aligning xticks"
    https://stackoverflow.com/questions/42615527/sequence-logos-in-matplotlib-aligning-xticks?rq=1
    asked by rightskewed
    answered by ImportanceOfBeingErnest
    modified for the negative functionality
    :param letter (str): A letter to plot in the PWM, has to be part of letters
    :param x (float): x position of the letter to plot
    :param y (float): y position of the letter to plot
    :param yscale (float): scaling parameter in y direction
    :param ax: matplotlib ax handler, usually not use
    :param negativ (bool): should the letter be plotted into the negative direction,
     use a nother color dict in this case
    :return: patch object for matplotlib plot
    '''
    text = LETTERS[letter]
    t = mpl.transforms.Affine2D().scale(1*globscale, yscale*globscale) + \
        mpl.transforms.Affine2D().translate(x,y) + ax.transData
    if negativ:
        p = PathPatch(text, lw=0, fc=COLOR_SCHEME_neg[letter],  transform=t)
    else:
        p = PathPatch(text, lw=0, fc=COLOR_SCHEME[letter],  transform=t)
    if ax != None:
        ax.add_artist(p)
    return p
def IC(row):
    '''

    :param row (array or list): one row of a PWM
    :return (float): information content
    '''
    ICs=[np.log2(x)*-x for x in row]
    IC_tot=np.log2(len(row))-np.sum(ICs)
    IC_return=[x*IC_tot for x in row]
    return IC_return

def mut_info(PWM,DPWM,DPWM_shift=1):
    '''
    for each position: calculate mutual information content by
    multiplying the entries in the PWM of subsequent positions and dividing the entry in the DPWM by this product
    the 2logarithm of this quotient is multiplied with the entry of the DPWM
    see also https://en.wikipedia.org/wiki/Mutual_information
    :param PWM (array): position weight matrix
    :param DPWM (array): dinucleotide position weight matrix (absolute)
    :param DPWM_shift (int): for how many positions in the DPWM shifted relative to the PWM (to the right)
    :return (array): mutial information content between two consequetive positions
    '''
    mut_info=[]
    for i in range(len(PWM)-1):
        mut_info.append([b[1] * np.log2(b[1]/b[0]) for b in zip([a[0]*a[1] for a in itertools.product(PWM[i],PWM[i+1])],DPWM[i+DPWM_shift])])
    return mut_info


class plotter:
    '''
    class to plot PWMs and DPWMs using matplotlib
    '''
    def __init__(self, PWM=None, DPWM=None, DPWM_shift=1):
        '''
        initialize class with PWM and DPWM
        :param PWM (array): position weight matrix
        :param DPWM (array): dinucleotide position weight matrix (absolute)
        :param DPWM_shift (int): for how many positions in the DPWM shifted relative to the PWM (to the right)
        '''
        self.PWM = PWM
        if PWM is not None:
            self.PWM_scores = [sorted(zip(['A', 'C', 'G', 'T'], IC(x)), key=operator.itemgetter(1)) for x in self.PWM]
        self.DPWM = DPWM
        if self.DPWM is not None:
            self.DPWM_scores = [
                sorted(zip(map(''.join, itertools.product('ACGT', 'ACGT')), x), key=operator.itemgetter(1)) for x in
                mut_info(self.PWM, self.DPWM, DPWM_shift=DPWM_shift)]

    def plot_PWM(self, filepath=None, maxi=2.0, space=1.0, ylabel='Information content [bits]'):
        '''
        function to plot the PWM
        :param filepath (str): path to save the resulting figure, including extension (e.g. ".png")
        :param maxi (float): Maximal information content to scale the plot
        :param space (float): how close should the letters be plotted in x direction
        :param ylabel (str): y axis label
        :return : PWM plot
        '''
        fig, ax = plt.subplots(figsize=(10, 3))
        x = 1.0
        for scores in self.PWM_scores:
            y = 0
            for base, score in scores:
                letterAt(base, x, y, score, ax)
                y += score
            x += space
            maxi = max(maxi, y)

        plt.xticks(np.arange(1, (len(self.PWM) + 2 - space) * space, space), np.arange(1, len(self.PWM) + 1))
        plt.xlim((0, x))
        plt.ylim((0, maxi))
        plt.ylabel(ylabel)
        plt.tight_layout()
        plt.show()
        if filepath is not None:
            plt.savefig(filepath, dpi=300)

    def plot_DPWM(self, filepath=None, maxi=1, mini=-1, ylabel='Mututal information [bits]', start=2):
        '''

        :param filepath (str): path to save the resulting figure, including extension (e.g. ".png")
        :param maxi (float): Maximal information content to scale the plot
        :param mini (float): minimal information content to scale the plot
        :param ylabel (str): y axis label
        :param start: start of numbering for positions in x axis, shift relative to PWM
        :return: DPWM plot
        '''
        fig, ax = plt.subplots(figsize=(10, 3))
        x = 1
        for scores in self.DPWM_scores:

            y = np.sum([a[1] for a in filter(lambda b: b[1] < 0, scores)])
            for base, score in scores:
                for z in range(2):
                    letterAt(base[z], x + z, y, abs(score), ax, negativ=score < 0)
                if score < 0:
                    y -= score
                elif score > 0:
                    y += score
                else:
                    y = 0
                # x += 1
            x += 3
        plt.xlim((0, len(self.DPWM_scores) * 3))
        plt.ylim(mini, maxi)
        plt.xticks(np.arange(1.5, len(self.DPWM_scores) * 3, 3), range(start, len(self.DPWM_scores) + start))
        plt.ylabel(ylabel)
        plt.tight_layout()
        plt.show()
        if filepath is not None:
            plt.savefig(filepath, dpi=300)