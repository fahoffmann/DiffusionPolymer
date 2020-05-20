def prepare_subplots(plt):
    plt.figure(1)
    ax1 = plt.subplot(211)
    ax2 = plt.subplot(212)
    return plt,ax1,ax2

def plot_fit_err(x,y,yfit,yerr,ax1,ax2,label):
    ax1.scatter(x, y, marker='o', label=label)
    ax1.plot(x, yfit, label='')
    ax2.scatter(x, yerr, marker='o', label=label)

def plot_fit(x,y,yfit,ax,label1,label2):
    ax.scatter(x, y, marker='o', label=label1)
    ax.plot(x, yfit, label=label2)

def plot_ax(ax,xmin,xmax,ymin,ymax,xlabel,ylabel):
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.legend()
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    return ax