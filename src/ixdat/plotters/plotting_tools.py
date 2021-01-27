def color_axis(ax, color, lr="right", xy="y"):
    ax.spines[lr].set_color(color)
    ax.tick_params(axis=xy, color=color)
    ax.tick_params(axis=xy, labelcolor=color)
    if xy == "y":
        ax.yaxis.label.set_color(color)
    if xy == "x":
        ax.xaxis.label.set_color(color)
