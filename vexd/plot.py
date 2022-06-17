import matplotlib as mpl
mpl.use('agg')
from matplotlib.cm import ScalarMappable
from matplotlib.colors import BoundaryNorm
from matplotlib.figure import Figure
import pandas as pd
import numpy as np
import io
import roman

def gene_boxplot(gene, gene_results):
    """
    Creates a boxplot showing how a particular gene changed in VExD

    `gene`: A dictionary containing at least 'symbol' and 'ensembl_id' keys,
            used to title the plot
    `gene_results`: The output of geosearch.get_results_by_gene()

    Returns a byte array containing a PNG image of the plot
    """
    # Sort by Baltimore category and then name, extract the data into
    # a format that mpl.axes.boxplot can use
    results = pd.DataFrame(gene_results)\
        .assign(baltimore=lambda df: fromRomanVector(df['virus_baltimore']))\
        .sort_values(by=['baltimore', 'virus'])
    labels = []
    data = []
    for name, subset in results.groupby('virus', sort=False):
        labels.append(name)
        data.append(subset['logfc'].to_list())
    
    # Main figure creation and plotting
    fig = Figure(figsize=(13, 8), layout="constrained")
    ax = fig.add_subplot()
    bplot = ax.boxplot(
        data,
        labels=labels,
        vert=False,
        patch_artist=True,
        widths=0.7,
        medianprops={'color': 'black'},
        showfliers=False,
        showcaps=False,
    )
    # Color the boxplots by genome composition
    for patch, genome_idx in zip(bplot['boxes'], results.drop_duplicates('virus').groupby('baltimore', sort=False).ngroup()):
        patch.set_facecolor(mpl.colormaps['Dark2'](genome_idx))
    # Plot the individual points over the boxplots, with a small amount of jitter
    for y, x in enumerate(data):
        ax.plot(
            x, 
            np.random.uniform(-0.2, 0.2, len(x))+y+1, 
            linestyle="none", 
            color='black', 
            markersize=3, 
            marker='o', 
            fillstyle='none',
        )
    # Set up the axes and labels
    ax.set_xlabel("$log_2$ Fold Change")
    ax.set_ylabel("")
    ax.set_yticklabels(
        labels, 
        fontdict={'fontstyle': 'italic'}
    )
    max_extent = max([abs(x) for x in ax.get_xlim()])
    ax.set_xlim(-max_extent, max_extent)
    ax.spines.top.set_visible(False)
    ax.invert_yaxis()
    # Hide the tick marks on the y axis
    for tick in ax.yaxis.get_major_ticks():
        tick.tick1line.set_markersize(0)
        tick.tick2line.set_markersize(0)
    # Add grid lines and the lines highlighting the 0 point and the change in genome composition
    ax.axvline(x=0, color='k', linestyle='-', linewidth=2, zorder=0)
    class_bounds = _get_class_boundaries(results)
    for bnd in class_bounds:
        ax.axhline(y=bnd+1, color='0.8', linestyle='-')
    ax.xaxis.grid(True)
    ax.grid(axis='y', linestyle=':', color='0.8')
    # Title and colorbar
    fig.suptitle(f"{gene['symbol']} ({gene['ensembl_id']})", fontsize="x-large")
    add_colorbar(fig, class_bounds, results['virus_genome'].drop_duplicates().to_list())
    # Finally, write out the PNG and return
    buffer = io.BytesIO()
    fig.savefig(buffer)
    return buffer.getvalue()

def fromRomanVector(v):
    """
    Translates a list of roman numerals (as strings) into integers

    Input: An iterable (list, pandas Series, etc) that contains strings, 
            each of which is expected to contain roman numerals
    
    Output: A list of integers, corresponding to the input.  Any inputs
            that could not be interpreted will yield None
    """
    result = []
    for val in v:
        try:
            result.append(roman.fromRoman(val))
        except roman.InvalidRomanNumeralError:
            result.append(None)
    return result

def _get_class_boundaries(data):
    """
    Internal function to find the boundaries between Baltimore classes in the data
    """
    results = []
    viruses = data.drop_duplicates('virus', ignore_index=True)
    balt_col = viruses.columns.get_loc('baltimore')
    for i in range(1, viruses.shape[0]):
        if viruses.iat[i, balt_col] != viruses.iat[i-1, balt_col]:
            results.append(i-0.5)
    return results

def add_colorbar(figure, boundaries, names):
    """
    Creates a colorbar that corresponds to the Baltimore categories in the input,
    adds it to the figure.
    """
    augmented_bounds = [-0.5,] + boundaries + [boundaries[-1]+1,]
    cbar = figure.colorbar(
        ScalarMappable(
            BoundaryNorm(augmented_bounds, len(augmented_bounds)-1), 
            mpl.colormaps['Dark2']
        ),
        spacing='proportional',
        ticks=[],
        fraction=0.02,
        aspect=50,
        pad=0
    )
    cbar.ax.invert_yaxis()
    cbar.set_ticks(
        pd.Series(augmented_bounds).rolling(2).mean().tail(-1).to_list(),
        labels=names,
        minor=True
    )
