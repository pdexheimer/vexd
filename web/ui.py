from flask import Blueprint, render_template
from webargs import fields
from webargs.flaskparser import use_kwargs
from .geo import geo

bp = Blueprint('ui', __name__, url_prefix='/')

@bp.route('/geo/', defaults={'gse_id': None})
@bp.route('/geo/<gse_id>')
@use_kwargs({'gse_id': fields.Str()}, location='query')
def display_geo_entry(gse_id):
    return render_template('gse.html', gse=geo().get_gse_info(gse_id))

@bp.route('/gene')
@use_kwargs({'q': fields.Str()}, location='query')
def display_gene_info(q):
    return render_template('gene.html', gene=geo().get_gene_info(q), gene_results=geo().get_results_by_gene(q))

@bp.route('geo_lookup/')
def geo_search():
    return render_template('geo_search.html')

@bp.route('gene_lookup/')
def gene_search():
    return render_template('gene_search.html')

@bp.route('/')
def app_home():
    return render_template('uber_search.html')

@bp.route('/search')
@use_kwargs({'q': fields.Str()}, location='query')
def search_results(q):
    return render_template('search_results.html', gse_list=geo().search_gse(virus_name=q))

@bp.route('/statistics')
def statistics():
    return render_template('stats.html', virus_counts=geo().count_by_virus())