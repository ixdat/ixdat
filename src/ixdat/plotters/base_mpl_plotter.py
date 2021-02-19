from matplotlib import pyplot as plt


class MPLPlotter:

    def new_ax(self, xlabel=None, ylabel=None):
        fig, ax = plt.subplots()
        if xlabel:
            ax.set_xlabel(xlabel)
        if ylabel:
            ax.set_ylabel(ylabel)
        return ax