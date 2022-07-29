import numpy as np
import matplotlib

def plot_discrete_count_distribution(data, color, plt_handle):
    xlabels, counts = np.unique(data, return_counts=True)
    plt_handle.bar(xlabels, counts, align='center', color=color)
    if issubclass(type(plt_handle), matplotlib.axes.SubplotBase):
        plt_handle.set_xticks(xlabels)
    else:
        plt_handle.xticks(xlabels)
