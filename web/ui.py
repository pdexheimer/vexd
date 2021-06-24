from flask import Blueprint, render_template
from flask.helpers import make_response
from webargs import fields
from webargs.flaskparser import use_kwargs
from .geo import geo

bp = Blueprint('ui', __name__, url_prefix='/')

@bp.route('/')
def vexd_home():
    return render_template('homepage.html')

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

@bp.route('/geo/', defaults={'gse_id': None})
@bp.route('/geo/<gse_id>')
@use_kwargs({'gse_id': fields.Str()}, location='query')
def display_geo_entry(gse_id):
    return render_template('gse.html', gse=geo().get_gse_info(gse_id))

@bp.route('/gene')
@use_kwargs({'q': fields.Str()}, location='query')
def display_gene_info(q):
    return render_template('gene.html', gene=geo().get_gene_info(q), gene_results=geo().get_results_by_gene(q))

@bp.route('/geo_lookup')
def geo_lookup():
    return render_template('geo_search.html')

@bp.route('/gene_lookup')
def gene_search():
    return render_template('gene_search.html')

@bp.route('/geo_search')
def geo_search():
    return render_template('uber_search.html')

@bp.route('/search')
@use_kwargs({'q': fields.Str()}, location='query')
def search_results(q):
    return render_template('search_results.html', gse_list=geo().search_gse(virus_name=q))

@bp.route('/statistics')
def statistics():
    return render_template('stats.html', virus_counts=geo().count_by_virus())

@bp.route('/heatmap')
@use_kwargs({'q': fields.Str()}, location='query')
def heatmap(q):
    return render_template('morpheus.html', 
        virus=q,
        genes=geo().get_significant_genes_by_virus(virus_name=q),
        studies=geo().get_studies_by_virus(q)
    )

@bp.route('/methods')
def methods():
    return render_template('methods.html')

@bp.route('/about')
def about():
    return render_template('homepage.html')

@bp.route('/help')
def help():
    return render_template('homepage.html')
