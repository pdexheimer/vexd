from flask import Blueprint, abort, jsonify
from webargs import fields, validate
from webargs.flaskparser import use_args

from .geo import geo

bp = Blueprint('api', __name__, url_prefix='/api')

@bp.route('/virus_info/<name>')
def virus_info(name):
    result = geo().get_virus(name)
    if result is None:
        abort(404, f'{name} is not a recognized virus species')
    return result

@bp.route('/virus_typeahead/<text>')
def virus_typeahead(text):
    return jsonify(geo().find_virus_by_prefix(text))

@bp.route('/gene_typeahead/<text>')
def gene_typeahead(text):
    return jsonify(geo().find_gene_by_prefix(text))

@bp.route('/public_data/<id>')
def geo_lookup(id):
    id = id.upper()
    result = geo().get_gse_info(id)
    if result is None:
        abort(404, f'{id} has not been curated')
    return result

geo_search_args = {
    'virus': fields.Str(missing=None),
    'platform': fields.Str(missing=None, validate=validate.OneOf(['microarray', 'ngs'])),
    'with_pubmed': fields.Bool(missing=False),
    'only_valid': fields.Bool(missing=True),
    'human_only_viruses': fields.Bool(missing=True),
}

@bp.route('/public_data/')
@use_args(geo_search_args, location='query')
def geo_search(args):
    result = geo().search_gse(
        virus_name=args['virus'],
        platform=args['platform'],
        with_pubmed=args['with_pubmed'],
        only_valid=args['only_valid'],
        human_only_virus=args['human_only_viruses']
    )
    if result is None:
        return jsonify([])
    return jsonify(result)