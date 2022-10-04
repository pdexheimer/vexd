import numpy as np
import pandas as pd
import click
from flask import current_app

# A Mongo aggregation pipeline to get all differential expression results
_deg_aggregation = [
    {
        '$match': {
            'ensembl_id': { '$ne': '' }
        }
    }, {
        '$group': {
            '_id': {
                'study': '$study', 
                'virus': '$virus', 
                'bto_id': '$bto_id', 
                'bto_name': '$bto_name', 
                'platform': '$platform'
            }, 
            'results': {
                '$push': {
                    'ensembl_id': '$$ROOT.ensembl_id', 
                    'logfc': '$$ROOT.logfc', 
                    'adj_p': '$$ROOT.adj_p'
                }
            }
        }
    }, {
        '$project': {
            'study': '$_id.study', 
            'virus': '$_id.virus', 
            'bto_id': '$_id.bto_id', 
            'bto_name': '$_id.bto_name', 
            'platform': '$_id.platform', 
            '_id': 0, 
            'results': 1
        }
    }
]

# A Mongo aggregation to get all differential expression analyses
# similar to geosearch.search_studies()
_analysis_aggregation = [
    { 
        '$match': {'samples.valid_experiment': True}
    }, {
        '$unwind': {'path': '$samples'}
    }, {
        '$match': {
            'samples.valid_experiment': True, 
            'samples.normalized_virus': {'$ne': 'Uninfected'}
        }
    }, {
        '$group': {
            '_id': {
                'study': '$id', 
                'bto_name': '$samples.bto_name', 
                'virus': '$samples.normalized_virus', 
                'platform': '$samples.cdf', 
                'bto_id': '$samples.bto_id'
            }, 
            'total_samples': {'$sum': 1}, 
            'usable_samples': {'$sum': {'$cond': [{'$eq': [{'$type': '$samples.missing_reason'}, 'missing']}, 1, 0]}}, 
            'cell_type': {'$first': '$samples.cell_type'}
        }
    }, {
        '$lookup': {
            'from': 'geo', 
            'as': 'controls', 
            'let': {
                'study': '$_id.study', 
                'cell_type': '$_id.bto_name', 
                'platform': '$_id.platform'
            }, 
            'pipeline': [
                { 
                    '$match': {'$expr': {'$eq': ['$id', '$$study']}}
                }, {
                    '$unwind': {'path': '$samples'}
                }, {
                    '$match': {'$expr': {
                            '$and': [
                                {'$eq': ['$samples.normalized_virus', 'Uninfected']}, 
                                {'$eq': ['$samples.valid_experiment', True]}, 
                                {'$eq': ['$id', '$$study']},
                                {'$eq': ['$samples.cdf', '$$platform']},
                                {'$eq': ['$samples.bto_name', '$$cell_type']}
                            ]
                    } }
                }, {
                    '$replaceRoot': {'newRoot': '$samples'}
                }, {
                    '$group': {
                        '_id': None, 
                        'total_controls': {'$sum': 1}, 
                        'usable_controls': {'$sum': {'$cond': [{'$eq': [{'$type': '$missing_reason'}, 'missing']}, 1, 0]}}
                    }
                }
            ]
        }
    }, {
        '$lookup': {
            'from': 'results', 
            'as': 'gene_count', 
            'let': {
                'study': '$_id.study', 
                'cell_type': '$_id.bto_name', 
                'platform': '$_id.platform', 
                'virus': '$_id.virus'
            }, 
            'pipeline': [
                {
                    '$match': {'$expr': {
                            '$and': [
                                {'$eq': ['$$study', '$study']},
                                {'$eq': ['$$cell_type', '$bto_name']},
                                {'$eq': ['$$platform', '$platform']},
                                {'$eq': ['$$virus', '$virus']},
                                {'$lte': ['$adj_p', 0.05]},
                                {'$or': [
                                    {'$gte': ['$logfc', 1]},
                                    {'$lte': ['$logfc', -1]}

                                ]}
                            ]
                    } }
                }, {
                    '$group': {
                        '_id': None, 
                        'signif': {'$sum': 1}, 
                        'down': {'$sum': {'$cond': [{'$lt': ['$logfc', 0]}, 1, 0]}}, 
                        'up': {'$sum': {'$cond': [{'$gt': ['$logfc', 0]}, 1, 0]}}
                    }
                }, {
                    '$project': {'_id': 0}
                }
            ]
        }
    }, {
        '$lookup': {
            'from': 'geo', 
            'as': 'geo', 
            'let': {'studyid': '$_id.study'}, 
            'pipeline': [
                {
                    '$match': {'$expr': {'$eq': ['$id', '$$studyid']}}
                }, {
                    '$project': {'_id': 0, 'title': 1}
                }
            ]
        }
    }, {
        '$replaceRoot': {
            'newRoot': {
                '$mergeObjects': [
                    '$_id', 
                    {'$arrayElemAt': ['$gene_count', 0]}, 
                    {'$arrayElemAt': ['$geo', 0]}, 
                    {'$arrayElemAt': ['$controls', 0]}, 
                    {
                        'usable_samples': '$usable_samples', 
                        'total_samples': '$total_samples', 
                        'cell_type': '$cell_type'
                    }
                ]
            }
        }
    }, {
        '$set': {
            'enough_samples': {'$and': [
                    {'$gte': ['$usable_samples', 2]}, 
                    {'$gte': ['$usable_controls', 2]}
            ]}, 
            'study_order': {'$toInt': {'$ltrim': {'input': '$study', 'chars': 'GSE'}}}
        }
    }, {
        '$sort': {
            'signif': -1, 
            'enough_samples': -1, 
            'study_order': -1
        }
    }, {
        '$project': {
            '_id': 0, 
            'enough_samples': 0, 
            'study_order': 0
        }
    }
]

def _write_dataframe(df, basename, directory=None):
    path = basename if directory is None else f'{directory}/{basename}'
    df = df.dropna(how='all')
    df.to_csv(f'{path}.tsv', sep='\t')
    df.to_parquet(f'{path}.parquet')


def prepare_db_dumps(geo, directory):
    click.echo(f'Preparing to create download files in {directory}')
    db = geo.client.vexd
    click.echo('Collecting DEG results...')
    result = db.results.aggregate(_deg_aggregation, allowDiskUse=True)
    click.echo('Creating result matrix...')

    summary = None
    for comp in result:
        info = pd.DataFrame(
            [ (x['logfc'], x['adj_p']) for x in comp['results'] ],
            index=[ x['ensembl_id'] for x in comp['results'] ],
            columns=pd.MultiIndex.from_frame(pd.DataFrame({
                'study': [comp['study'],]*2,
                'virus': [comp['virus'],]*2,
                'bto_id': [comp['bto_id'],]*2,
                'bto_name': [comp['bto_name'],]*2,
                'platform': [comp.get('platform', 'RNA-seq'),]*2,
                'data': ['log_fc', 'adj_p'],
            }))
        )
        if summary is None:
            summary = info
        else:
            summary = summary.merge(info, how='outer', left_index=True, right_index=True)

    summary.rename_axis(index='ensembl_id', inplace=True)
    # Add the gene symbol as a level on the index
    result = db.genes.aggregate([{'$project': {'_id': 0, 'ensembl_id': 1, 'symbol': 1}}])
    gene_lookup = pd.DataFrame(list(result)).set_index('ensembl_id')
    summary, gene_lookup = summary.align(gene_lookup, join='left', axis=0, copy=False)
    summary.set_index(gene_lookup['symbol'], append=True, inplace=True)

    click.echo('Writing DEG files...')
    _write_dataframe(summary, 'all_viruses', directory)
    for virus, group in summary.groupby('virus', axis=1):
        _write_dataframe(group, virus.replace(' ', '_'), directory)

    click.echo('Collecting metadata...')
    result = db.geo.aggregate(_analysis_aggregation)
    click.echo('Writing files...')
    _write_dataframe(pd.DataFrame(list(result)), 'deg_analyses', directory)
    click.echo('Done creating downloads')

def extrapolate(max_rank, max_val, max_reduced_rank, max_reduced_val, stride):
    return max_reduced_val + (stride/(max_rank-max_reduced_rank))*(max_val-max_reduced_val)

def deduplicate(df, col_with_dupes='logfc', col_to_interp='megarank'):
    dedup = pd.DataFrame(df)
    dedup.insert(0, 'duplicated', df[col_with_dupes].duplicated(keep=False))
    dedup.drop_duplicates(subset=col_with_dupes, keep='first', inplace=True)
    dedup.loc[dedup['duplicated'], col_to_interp] = np.nan
    dedup[col_to_interp].interpolate(inplace=True)
    return dedup.drop(columns='duplicated')

def reduce_distribution(values, stride=1000):
    """
    Given a sorted Series of values, computes a reduced distribution of the same
    values.  The intention is that the reduced distribution is much smaller and
    easier to work with, and still provides a reasonable base for imputing the
    missing data later on.  The final (len(values) % stride) values are extrapolated
    to the next appropriate data point

    Returns a DataFrame with (rank/stride) as the index and the appropriate values
    in a column.  The column's name will be values.name from the input Series
    (defaulting to 'logfc')
    """
    reduced_values = values.iloc[::stride]
    extrapolated_val = extrapolate(len(values)-1, values.iloc[-1], len(reduced_values)-1, reduced_values.iloc[-1], stride)
    # Add the extrapolated value and simultaneously reset the reduced index
    reduced_values = pd.concat([reduced_values, pd.Series([extrapolated_val,])], ignore_index=True)
    df = deduplicate(pd.DataFrame({'logfc': reduced_values, 'megarank': reduced_values.index.to_series()}))
    if values.name is not None:
        df.rename(columns={'logfc': values.name}, inplace=True)
    return df.set_index('megarank')

def save_background_distribution(geo):
    click.echo("Retrieving background distribution...")
    all_results = pd.DataFrame(geo.get_all_results())\
        .drop('ensembl_id', axis='columns')\
        .set_index(['virus', 'bto_id'])
    click.echo("Sorting distribution...")
    all_results = all_results.sort_values(by='logfc', ascending=True)
    click.echo("Subsetting and saving...")
    background = reduce_distribution(all_results['logfc'])
    subsets = [('N/A', 'N/A')]
    subset_sizes = [all_results.shape[0],]
    click.echo("  Subsetting by virus...")
    for g, subdf in all_results.groupby('virus', sort=False):
        background = background.join(reduce_distribution(subdf['logfc']), rsuffix='a')
        subsets.append((g,'N/A'))
        subset_sizes.append(subdf.shape[0])
    click.echo("  Subsetting by cell type...")
    for g, subdf in all_results.groupby('bto_id', sort=False):
        background = background.join(reduce_distribution(subdf['logfc']), rsuffix='a')
        subsets.append(('N/A',g))
        subset_sizes.append(subdf.shape[0])
    click.echo("  Subsetting by virus x cell type...")
    for g, subdf in all_results.groupby(['virus', 'bto_id'], sort=False):
        background = background.join(reduce_distribution(subdf['logfc']), rsuffix='a')
        subsets.append(g)
        subset_sizes.append(subdf.shape[0])
    click.echo("  Saving...")
    subset_index = pd.MultiIndex.from_tuples(subsets, names=['virus', 'bto_id'])
    background.columns = subset_index
    with current_app.open_instance_resource('background.csv', 'w') as f:
        background.to_csv(f, index=True)
    with current_app.open_instance_resource('counts.csv', 'w') as f:
        pd.Series(subset_sizes, index=subset_index, name="count").to_csv(f)
    click.echo("Background distribution saved")

