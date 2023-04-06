import base64
from collections import Counter
import re

from flask import Blueprint, current_app, redirect, render_template, send_from_directory, url_for
from flask.helpers import make_response
import pandas as pd
from webargs import fields, validate
from webargs.flaskparser import use_kwargs

from . import plot
from .geo import geo
from .utils import enrichment, parse_dl_files

bp = Blueprint('ui', __name__, url_prefix='/')

#########################################
## Top level pages
#########################################

# Home page
@bp.route('/')
def vexd_home():
    return redirect(url_for('ui.gene_search'))

# Genes
@bp.route('/gene_lookup')
def gene_search():
    return render_template(
        'gene_search.html',
        study_count=geo().num_studies(),
        virus_count=geo().num_viruses()
    )

# Enrichment
@bp.route('/enrich')
def multiple_genes():
    return render_template(
        'multiple_genes.html',
        virus_tissue=geo().virus_tissue_combos()
    )

# Studies
@bp.route('/study')
def study_search():
    return render_template('study_search.html')

# API
@bp.route('/api')
def api_docs():
    return render_template('api_docs.html')

# Downloads
@bp.route('/downloads')
def download_page():
    return render_template('downloads.html', 
        file_info=parse_dl_files(current_app.config['DOWNLOAD_DIRECTORY'])
    )

# Methods
@bp.route('/methods')
def methods():
    return render_template('methods.html')

# About VExD
@bp.route('/about')
def about():
    return render_template('about.html', virus_counts=geo().count_by_virus())

# Help
@bp.route('/help')
def help():
    return render_template('help.html')

#########################################
## (Single) Gene-related pages
#########################################

# Results page
@bp.route('/gene')
@use_kwargs({
    'q': fields.Str(),
    'tissue': fields.Str(missing=""),
    'descendants': fields.Bool(missing=True)
}, location='query')
def display_gene_info(q, tissue, descendants):
    return render_template(
        'gene.html', 
        gene=geo().get_gene_info(q), 
        gene_results=geo().get_results_by_gene(q, tissue, descendants),
        tissue=tissue,
        use_descendants=descendants,
        tissue_name=geo().get_name_of_bto_term(tissue)
    )

# Boxplot figure
@bp.route('/gene_plot')
@use_kwargs({
    'q': fields.Str(),
    'tissue': fields.Str(missing=""),
    'descendants': fields.Bool(missing=True)
}, location='query')
def gene_boxplot(q, tissue, descendants):
    response = make_response(
        plot.gene_boxplot(
            geo().get_gene_info(q),
            geo().get_results_by_gene(q, tissue, descendants)
        )
    )
    response.content_type = 'image/png'
    return response

# Gene Results as text (for download)
@bp.route('/gene/txt')
@use_kwargs({
    'q': fields.Str(),
    'tissue': fields.Str(missing=""),
    'descendants': fields.Bool(missing=True)
}, location='query')
def gene_info_text(q, tissue, descendants):
    gene_info = geo().get_gene_info(q)
    response = make_response(render_template(
        'gene_results.txt', 
        gene=gene_info, 
        gene_results=geo().get_results_by_gene(q, tissue, descendants)
    ))
    response.headers.set('Content-disposition', 'attachment', filename=f"{gene_info['symbol']}_results.txt")
    response.content_type = 'text/plain'
    return response

#########################################
## (Multiple Gene) Enrichment-related pages
#########################################

# Heatmap Generation/Statistics Return
@bp.route('/gene/multiple', methods=['POST'])
@use_kwargs({
    'raw_genes': fields.Str(), 
    'remove_unknown': fields.Boolean(missing=False),
    'sort_genes': fields.Boolean(missing=False),
    'virus': fields.Str(),
    'bto_id': fields.Str()
}, location='form')
def gene_heatmap(raw_genes, remove_unknown, sort_genes, virus, bto_id):
    gene_list = []
    escaped = raw_genes.translate(str.maketrans({
        "-": r"\-", "*": r"\*", "+": r"\+", "^": r"\^", "$": r"\$",
        "(": r"\(", ")": r"\)", "[": r"\[", "]": r"\]", ".": r"\.", "\\": r"\\"
    }))
    # Get the unique input symbols, but maintain their order (so a set won't work)
    input_genes = list(Counter(re.split(r'[\r\n ,;]+', escaped)))
    for gene_id in input_genes:
        if len(gene_list) >= 150:
            break
        if gene_id == '':
            continue
        gene_info = geo().find_gene(gene_id)
        if not (remove_unknown and gene_info is None):
            gene_list.extend(gene_info if gene_info is not None else [{'ensembl_id': 'Unknown', 'symbol': gene_id}])
    gene_results = pd.DataFrame(geo().get_multiple_gene_results([g['ensembl_id'] for g in gene_list if g['ensembl_id'] != 'Unknown'], virus, bto_id))
    if gene_results.empty:
        return make_response({
            'num_genes': 0,
            'num_measurements': 0,
            'num_background': 'NA',
            'support': 'NA',
            'pval': 'NA',
            'heatmap': base64.b64encode(b'').decode(),
        })
    enrich = enrichment(gene_results['logfc'], virus, bto_id)
    return make_response({ 
        'num_genes': gene_results['ensembl_id'].nunique(),
        'num_measurements': enrich['n1'],
        'num_background': enrich['n2'],
        'support': enrich['support'],
        'pval': enrich['pvalue'],
        'heatmap': base64.b64encode(plot.gene_heatmap(gene_list, gene_results, sort_genes)).decode(),
    })

#########################################
## Study-related pages
#########################################

# Search Results
@bp.route('/virus')
@use_kwargs({'q': fields.Str()}, location='query')
@use_kwargs({'tissue': fields.Str()}, location='query')
@use_kwargs({'descendants': fields.Bool(missing=False)}, location='query')
def virus_results(q, tissue, descendants):
    return render_template(
        'combined_results.html', 
        virus=q, 
        tissue=tissue,
        search_kids=descendants,
        search_results=geo().search_studies(q, tissue, descendants)
    )

# Download single-study results
result_args = {
    'virus': fields.Str(),
    'study': fields.Str(),
    'platform': fields.Str(),
    'cell_type': fields.Str()
}
@bp.route('/results/<geneSet>')
@use_kwargs({'geneSet': fields.Str(validate=validate.OneOf(['all','sig','up','down']))}, location='view_args')
@use_kwargs(result_args, location='query')
def download_results(virus, study, platform, cell_type, geneSet):
    response = make_response(render_template('results.txt', 
        results=geo().get_analysis_results(study, virus, cell_type, platform, geneSet)))
    response.headers.set('Content-disposition', 'attachment', filename="results.txt")
    response.content_type = 'text/plain'
    return response

# Search Results as text (for download)
@bp.route('/virus/txt')
@use_kwargs({'q': fields.Str()}, location='query')
@use_kwargs({'tissue': fields.Str()}, location='query')
@use_kwargs({'descendants': fields.Bool(missing=False)}, location='query')
def virus_results_as_text(q, tissue, descendants):
    response = make_response(render_template(
        'combined_results.txt', 
        virus=q, 
        tissue=tissue,
        search_kids=descendants,
        search_results=geo().search_studies(q, tissue, descendants)
    ))
    response.headers.set('Content-disposition', 'attachment', filename="studies.txt")
    response.content_type = 'text/plain'
    return response

# Single GEO Entry page
@bp.route('/geo/', defaults={'gse_id': None})
@bp.route('/geo/<gse_id>')
@use_kwargs({'gse_id': fields.Str()}, location='query')
def display_geo_entry(gse_id):
    return render_template('gse.html', gse=geo().get_gse_info(gse_id))

#########################################
## Download-related pages
#########################################

# Individual file download
@bp.route('/download/<filename>')
@use_kwargs({'filename': fields.Str()}, location='view_args')
def download_file(filename):
    return send_from_directory(
        current_app.config['DOWNLOAD_DIRECTORY'],
        filename,
        as_attachment=True
    )
