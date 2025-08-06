from itertools import combinations, product

import matplotlib as mpl
import matplotlib.pyplot as plt 
import numpy as np
from models.plot import Item

class PlotAgent:
    def setScatter(self, axs, x, y, label):
        axs.scatter(x, y, label=label, color="red")
        axs.legend(loc='best')
        pass
    def addPlot(self, axs, x, y, label):
        axs.plot(x, y, label=label, color="red")
        axs.legend()
    def setPlot(self, axs, x, y, label, tags):
        xlabel, ylabel = tags

        axs.plot(x, y, label=ylabel+f"({xlabel})", color="black")
        axs.set_title(label=label, size=10)
        axs.set_ylabel(ylabel+f"({xlabel})", rotation=0, size=10)
        axs.set_xlabel(xlabel, rotation=0, size=10)
        axs.yaxis.set_label_coords(-.13, .95)
        axs.legend(loc='best')
        # axs.yaxis.set_label_coords(-.1, .95)
        # axs.xaxis.set_label_coords(1.05, -0.025)
        pass
    def __init__(self, *args, **kwargs) -> None:

        t, xs, u = args
        
        labels = {  0:u"Координата материальной точки",
                    1:u"Скорость материальной точки",
                    2:""}
        mpl.rcParams.update(mpl.rcParamsDefault)
        plt.style.use('seaborn-v0_8-whitegrid')
        if isinstance(xs[0], np.ndarray):
            self.fig, self.axs = plt.subplots(2,len(xs[0]), figsize=(9, 7))
            for i in range(len(xs[0])):
                self.setPlot(self.axs[0][i], t, xs[0:,i], labels[i], ('t',f'$x_{i+1}$'))
            for i, j in combinations(range(len(xs[0])), 2):
                self.setPlot(self.axs[1][i], xs[0:,i], xs[0:,j], fr"Фазовый портрет $x_{i+1}$ $x_{j+1}$", (f'$x_{i+1}$', f'$x_{j+1}$'))
            self.setPlot(self.axs[-1][-1], t, u, u"Оптимальное управление", ('t', 'u'))
        self.fig.tight_layout()
        plt.style.use("default")
        pass
    def getFig(self):
        return self.fig, self.axs
    pass

class Plot:
    def __init__(self):
        pass
    def post(self,tk=None, xn=None, un=None):
        print(tk)
        plot = PlotAgent(tk, xn, un)
        return plot.fig
        
