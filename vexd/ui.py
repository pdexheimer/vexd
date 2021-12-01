from flask import Blueprint, current_app, render_template, send_from_directory
from flask.helpers import make_response
from webargs import fields, validate
from webargs.flaskparser import use_kwargs
from .geo import geo
from .utils import parse_dl_files

bp = Blueprint('ui', __name__, url_prefix='/')

@bp.route('/')
def vexd_home():
    return render_template(
        'homepage.html',
        study_count=geo().num_studies(),
        virus_count=geo().num_viruses()
    )

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

@bp.route('/geo/', defaults={'gse_id': None})
@bp.route('/geo/<gse_id>')
@use_kwargs({'gse_id': fields.Str()}, location='query')
def display_geo_entry(gse_id):
    return render_template('gse.html', gse=geo().get_gse_info(gse_id))

@bp.route('/gene')
@use_kwargs({'q': fields.Str()}, location='query')
def display_gene_info(q):
    return render_template('gene.html', gene=geo().get_gene_info(q), gene_results=geo().get_results_by_gene(q))

@bp.route('/gene/txt')
@use_kwargs({'q': fields.Str()}, location='query')
def gene_info_text(q):
    gene_info = geo().get_gene_info(q)
    response = make_response(render_template(
        'gene_results.txt', 
        gene=gene_info, 
        gene_results=geo().get_results_by_gene(q)
    ))
    response.headers.set('Content-disposition', 'attachment', filename=f"{gene_info['symbol']}_results.txt")
    response.content_type = 'text/plain'
    return response

@bp.route('/gene_lookup')
def gene_search():
    return render_template('gene_search.html')

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

@bp.route('/methods')
def methods():
    return render_template('methods.html')

@bp.route('/downloads')
def download_page():
    return render_template('downloads.html', 
        file_info=parse_dl_files(current_app.config['DOWNLOAD_DIRECTORY'])
    )

@bp.route('/download/<filename>')
@use_kwargs({'filename': fields.Str()}, location='view_args')
def download_file(filename):
    return send_from_directory(
        current_app.config['DOWNLOAD_DIRECTORY'],
        filename,
        as_attachment=True
    )

@bp.route('/api')
def api_docs():
    return render_template('api_docs.html')

@bp.route('/about')
def about():
    return render_template('about.html', virus_counts=geo().count_by_virus())

@bp.route('/help')
def help():
    return render_template('help.html')
