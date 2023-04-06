import matplotlib as mpl
mpl.use('agg')
from matplotlib.cm import ScalarMappable
from matplotlib.colors import BoundaryNorm, PowerNorm
from matplotlib.figure import Figure
from matplotlib.ticker import AutoLocator
import pandas as pd
import numpy as np
from scipy.stats import gaussian_kde, uniform
import io
import roman

def gene_boxplot(gene, gene_results):
    """
    Creates a boxplot showing how a particular gene changed in VExD.  Viruses
    are sorted and colored by their Baltimore category, and a density plot
    of all viruses is included at the top.

    Parameters
    ----------
    gene: A dictionary containing at least 'symbol' and 'ensembl_id' keys,
        used to title the plot
    gene_results: The output of geosearch.get_results_by_gene()

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
    ax = fig.add_subplot(10, 1, (2,10)) # Bottom 90%
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
    ax.spines.right.set_visible(False)
    ax.invert_yaxis()
    # Hide the tick marks on the y axis
    for tick in ax.yaxis.get_major_ticks():
        tick.tick1line.set_markersize(0)
        tick.tick2line.set_markersize(0)
    # Add grid lines and the lines highlighting the 0 point and the change in genome composition
    ax.axvline(x=0, color='k', linestyle='-', linewidth=2, zorder=0)
    class_bounds = _get_class_boundaries(results)
    for bnd in class_bounds:
        ax.axhline(y=bnd, color='0.8', linestyle='-')
    ax.xaxis.grid(True)
    ax.grid(axis='y', linestyle=':', color='0.8')
    # Title and colorbar
    fig.suptitle(f"{gene['symbol']} ({gene['ensembl_id']})", fontsize="x-large")
    add_colorbar(fig, class_bounds, results['virus_genome'].drop_duplicates().to_list())
    # Density plot of all points
    upper_ax = fig.add_subplot(10, 1, 1, sharex=ax, frame_on=False)
    density = gaussian_kde(results['logfc'].to_numpy())
    density.set_bandwidth(density.factor / 4.)
    density_x = np.linspace(ax.get_xlim()[0], ax.get_xlim()[1], 200)
    upper_ax.fill_between(density_x, density(density_x), facecolor='cadetblue')
    upper_ax.set_axis_off()
    # Finally, write out the PNG and return
    buffer = io.BytesIO()
    fig.savefig(buffer)
    return buffer.getvalue()

def gene_heatmap(gene_list, results, sort_genes=False):
    """
    Creates a pseudo-heatmap showing the response of multiple genes to infection.
    Genes are displayed in rows, and elements in each row are colored based on
    the estimated density of points at that position (from scipy.stats.gaussian_kde)

    Parameters
    ----------
    gene_list: A list of genes, where each gene is represented as a dictionary
        containing (at least) 'ensembl_id' and 'symbol' keys.
    results: A DataFrame containing at least the columns 'ensembl_id' and 'logfc'
    sort_genes: If True, the genes will be sorted such that the genes with the
        least change are last in the heatmap

    Returns a byte array containing a PNG image of the plot
    """
    # Generate density estimates for each gene over the entire range,
    # and stack them into a 2D numpy matrix.  Genes in gene_list that 
    # are not in the results will have a flat density of zero
    overall_density = gaussian_kde(results['logfc'].to_numpy())
    overall_density.set_bandwidth(overall_density.factor / 4.)
    max_extent = results['logfc'].abs().max()
    bygene = results.groupby('ensembl_id')
    xs = np.linspace(-max_extent, max_extent, 400)
    data = None
    for gene in gene_list:
        if gene['ensembl_id'] in bygene.groups:
            values = bygene.get_group(gene['ensembl_id'])['logfc'].to_numpy()
            if len(values) == 1:
                this_gene = np.zeros_like(xs)
                gene_idx = np.searchsorted(xs, values[0])
                this_gene[gene_idx] = 1
            else:
                # Add miniscule amounts of noise to prevent having a singular matrix
                # in the KDE when the values are all 0
                values += uniform.rvs(scale=1E-5, size=len(values))
                density = gaussian_kde(values)
                density.set_bandwidth(density.factor / 4.)
                this_gene = density(xs)
        else:
            this_gene = np.zeros_like(xs)
        if data is None:
            data = this_gene
        else:
            data = np.vstack((data, this_gene))
    if data.ndim == 1: # Only one gene - sorting doesn't matter, make the array 2D
        sorted_data = np.expand_dims(data, axis=0)
        ind = [0,]
    elif sort_genes:
        ind = np.lexsort(([sum(x[190:210]) for x in data],))
        sorted_data = [data[i] for i in ind]
    else:
        ind = range(len(gene_list))
        sorted_data = data
    
    # Create the figure, plot the density estimates, and add a vertical line at x=0
    fig = Figure(figsize=(13, 8), layout="constrained")
    ax = fig.add_subplot(10, 1, (2,10)) # Bottom 90%
    img = ax.imshow(sorted_data, cmap='plasma', aspect='auto', norm=PowerNorm(0.4), interpolation='none')
    ax.axvline(len(xs)/2, color='k', linestyle='--')
    fig.colorbar(img, ax=ax, label='Density')

    # Format the axes.  The x-axis is difficult because the underlying data is
    # indexed from 0 to (# of points density is estimated at), rather than
    # from -max_extent to max_extent (ie, fold change).  Resolve with a custom 
    # formatter (to display the correct value) and a custom locator (to place
    # ticks at nice values in fold-change space)
    ax.set_xlabel("$log_2$ Fold Change")
    fontsize = min([10, 13 - int(len(gene_list)/12)])
    ax.set_yticks(
        range(len(gene_list)), 
        [f'{gene_list[i]["symbol"]} ({gene_list[i]["ensembl_id"]})' for i in ind],
        fontsize=fontsize
    )
    ax.xaxis.set_major_locator(AffineLocator(xs[0], xs[-1], 0, len(xs)))
    get_fc = lambda x, _: "{:.2g}".format(xs[0] + ((xs[-1] - xs[0]) / len(xs)) * x)
    ax.xaxis.set_major_formatter(get_fc)

    # Density plot of the densities (meta-density plot?  density meta-plot?)
    upper_ax = fig.add_subplot(10, 1, 1, sharex=ax, frame_on=False)
    upper_ax.fill_between(range(len(xs)), overall_density(xs), facecolor='cadetblue')
    upper_ax.set_axis_off()

    # Save and return the PNG
    buffer = io.BytesIO()
    fig.savefig(buffer)
    return buffer.getvalue()

class AffineLocator(AutoLocator):
    """
    Tick Locator used when the underlying data is an affine transform of the values
    that should be shown to the user
    """

    def __init__(self, min_display, max_display, min_val, max_val):
        """
        min/max_display are the values that will be shown to the user
        min/max_val are the values that the underlying data will use
        """
        super().__init__()
        display_range = max_display - min_display
        val_range = max_val - min_val
        slope = float(display_range) / val_range
        # transform val coordinates to display coordinates
        self._xform = lambda x: min_display + slope * (x-min_val)
        # transform display coordinates back to val coordinates
        self._reverse = lambda x: min_val + (x - min_display) / slope
    
    def tick_values(self, vmin, vmax):
        return self._reverse(super().tick_values(self._xform(vmin), self._xform(vmax)))
    
    def view_limits(self, dmin, dmax):
        return self._reverse(super().view_limits(self._xform(dmin), self._xform(dmax)))

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
            results.append(i+0.5)
    return results

def add_colorbar(figure, boundaries, names):
    """
    Creates a colorbar that corresponds to the Baltimore categories in the input,
    adds it to the figure.
    """
    augmented_bounds = [figure.get_axes()[0].get_ylim()[1],] + boundaries + [figure.get_axes()[0].get_ylim()[0],]
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
