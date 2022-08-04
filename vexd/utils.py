from pathlib import Path

from flask import g, current_app
import numpy as np
import pandas as pd
from scipy.stats import norm

def file_info(path):
    if not path.exists():
        return None
    return {
        'name': path.name,
        'size': path.stat().st_size
    }

def collect_files(dir, basename):
    return {
        'TSV': file_info(dir / f'{basename}.tsv'),
        'PARQUET': file_info(dir / f'{basename}.parquet')
    }

def parse_dl_files(directory):
    result = {}
    dir = Path(directory)
    result['annot'] = collect_files(dir, 'deg_analyses')
    result['all'] = collect_files(dir, 'all_viruses')
    result['virus'] = {}
    for tsv in sorted(dir.glob('*.tsv')):
        if tsv.stem in set(['all_viruses', 'deg_analyses']):
            continue
        virus = tsv.stem.replace('_', ' ')
        result['virus'][virus] = collect_files(dir, tsv.stem)
    return result

def get_background_distribution():
    if 'bg' not in g:
        with current_app.open_instance_resource('background.csv') as f:
            g.bg = pd.read_csv(f)
    return g.bg

def get_background_size():
    if 'bg_size' not in g:
        with current_app.open_instance_resource('bg_count.txt') as f:
            g.bg_size = int(f.readline())
    return g.bg_size

def get_rank(val):
    bg = get_background_distribution()
    return np.rint(np.interp(val, bg['logfc'].to_numpy(), bg['megarank'].to_numpy()*1000))

def mannwhitney(vals):
    ranks = get_rank(vals)
    n1 = len(vals)
    n2 = get_background_size() - n1
    u1 = np.sum(ranks) - (n1 * (n1+1)) / 2
    u_mean = n1 * n2 / 2
    u = u1
    rho = u1 / (n1*n2)
    if u1 > u_mean:
        u1 = n1 * n2 - u1
    p = norm(loc=u_mean, scale=np.sqrt(u_mean * (n1 + n2 + 1)/6.0)).cdf(u1)*2
    return p, u, rho

def brunner_munzel(ranks, total_size):
    """
    Implements the Brunner-Munzel test to identify differences between ranked groups.  
    
    Brunner-Munzel is a non-parametric form of the Mann-Whitney U test.  Unlike
    Mann-Whitney, Brunner-Munzel does not assume that the variances in ranks are the
    same in each group, and it apparently handles unbalanced groups better.

    Differs from scipy.stats.brunnermunzel by optimizing for the case where one group
    is much smaller than the other.  Inputs are the ranks of items in the smaller group,
    and the size of both groups together.

    Returns a dictionary containing:
        
        * *statistic* - the test statistic, which has a standard normal distribution (for large numbers)
        * *pvalue* - the two-sided p-value corresponding to *statistic*
        * *support* - sometimes called 'stochastic superiority statistic' or
                    'common language effect size', this is the probability that
                    a randomly chosen member of group 1 (the input) will have a
                    higher rank than a randomly chosen member of group 2.
    """
    rank1 = np.sort(ranks)
    n1 = len(rank1)
    n2 = total_size - n1
    mean_rank1 = np.mean(rank1)
    exp_mean1 = (n1+1)/2.0
    mean_rank2 = exp_mean2 = (n2+1)/2.0
    # Now remove all of the values in rank1 from the mean of rank2
    # https://math.stackexchange.com/a/22351
    for i in range(n1):
        mean_rank2 += (mean_rank2 - rank1[i]) / (total_size-i-1)
    var1 = var2 = 0
    prev_rank = 0
    # Calculate the two variances
    # 10.1016/j.csda.2006.05.024
    for i in range(n1):
        var1 += (rank1[i] - i - mean_rank1 + exp_mean1) ** 2
        var2 += (rank1[i] - prev_rank) * (i - mean_rank2 + exp_mean2) ** 2
        prev_rank = rank1[i]
    var2 += (total_size - prev_rank) * (n1 - mean_rank2 + exp_mean2) ** 2
    var1 /= n1 - 1
    var2 /= n2 - 1
    # Finally, estimate the pooled variance (from the same Neubert&Brunner pub)
    Vn = np.sqrt(total_size*(var1/n2 + var2/n1))
    Tn = ((mean_rank2 - mean_rank1)/Vn)*np.sqrt(n1*n2/total_size)
    p = norm.cdf(-abs(Tn))*2
    return {
        'statistic': Tn, 
        'pvalue': p, 
        'support': (mean_rank1 - exp_mean1)/n2, 
        'n1': n1,
        'n2': n2,
    }

def enrichment(logfc_list):
    """
    Given a set of log fold changes, computes the enrichment of them compared to
    the background of all computed logfc's in the database
    """
    return brunner_munzel(get_rank(logfc_list), get_background_size())
