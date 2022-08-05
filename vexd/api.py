import re

from flask import Blueprint, abort, jsonify
from webargs import fields, validate
from webargs.flaskparser import use_args, use_kwargs

from .geo import geo

bp = Blueprint('api', __name__, url_prefix='/api')

@bp.errorhandler(400)
@bp.errorhandler(404)
@bp.errorhandler(422)
def json_error(e):
    return jsonify(error=str(e)), e.code

# Typeahead endpoints for autocomplete boxes
# Not part of the public API

@bp.route('/virus_typeahead/<text>')
def virus_typeahead(text):
    return jsonify(geo().find_virus_by_prefix(text))

@bp.route('/gene_typeahead/<text>')
def gene_typeahead(text):
    return jsonify(geo().find_gene_by_prefix(text))

# Other endpoints used for dynamical web-ness
# Not part of the public API

result_args = {
    'virus': fields.Str(),
    'study': fields.Str(),
    'platform': fields.Str(),
    'cell_type': fields.Str()
}
@bp.route('/preview/<geneSet>')
@use_kwargs({'geneSet': fields.Str(validate=validate.OneOf(['sig','up','down']))}, location='view_args')
@use_kwargs(result_args, location='query')
def analysis_preview(virus, study, platform, cell_type, geneSet):
    return jsonify(geo().get_analysis_results(study, virus, cell_type, platform, geneSet, 10))

@bp.route('/random/gene')
@use_kwargs({'num': fields.Int()}, location='query')
def random_geneids(num=10):
    return jsonify([ g['_id'] for g in geo().random_genes(num) ])

# Virus Endpoints

@bp.route('/v1/virus/name/<name>')
def virus_info(name):
    result = geo().get_virus(name)
    if result is None:
        abort(404, description=f'{name} is not a recognized virus species')
    return result

@bp.route('/v1/virus/alias/<alias>')
def virus_lookup(alias):
    return jsonify({
        'alias': alias,
        'results': geo().resolve_virus_alias(alias)
    })

@bp.route('/v1/virus', methods=['POST'])
@use_kwargs({'name': fields.Str(), 'alias': fields.Str()}, location='json')
def meta_virus_lookup(**kwargs):
    if 'name' not in kwargs and 'alias' not in kwargs:
        abort(400, description="Must specify either 'name' or 'alias'")
    if 'name' in kwargs and 'alias' in kwargs:
        abort(400, description="Must specify one of 'name' or 'alias', not both")
    if 'name' in kwargs:
        return virus_info(kwargs['name'])
    return virus_lookup(kwargs['alias'])

# GEO Endpoints

@bp.route('/v1/geo/<id>')
def geo_lookup(id):
    id = id.upper()
    result = geo().get_gse_info(id)
    if result is None:
        abort(404, description=f'{id} has not been curated')
    return result

# Search Endpoints

study_search_args = {
    'virus': fields.Str(missing=''),
    'bto_id': fields.Str(
        missing='', 
        validate=validate.Regexp('BTO:\d+', flags=re.IGNORECASE)#, error="'{input}' is not a valid BTO ID (should be BTO:xxxxxx)")
    ),
    'include_bto_children': fields.Bool(missing=True)
}

@bp.route('/v1/search', methods=["POST"])
@use_args(study_search_args)
def search_studies(args):
    return jsonify(geo().search_studies(args['virus'], args['bto_id'].upper(), args['include_bto_children']))

result_args = {
    'study': fields.Str(
        missing='',
        validate=validate.Regexp('GSE\d+', flags=re.IGNORECASE)
    ),
    'virus': fields.Str(missing=''),
    'bto_id': fields.Str(
        missing='', 
        validate=validate.Regexp('BTO:\d+', flags=re.IGNORECASE)#, error="'{input}' is not a valid BTO ID (should be BTO:xxxxxx)")
    ),
    'platform': fields.Str(missing=''),
    'gene': fields.Str(
        missing='',
        validate=validate.Regexp('ENSG\d+', flags=re.IGNORECASE)
    ),
    'set': fields.Str(
        missing='sig',
        validate=validate.OneOf(['all', 'sig', 'up', 'down'])
    )
}

@bp.route('/v1/results', methods=['POST'])
@use_args(result_args)
def search_results(args):
    if args['gene'] and (args['study'] or args['virus'] or args['bto_id'] or args['platform']):
        abort(400, description='Either search by gene or analysis (study/virus/bto), not both')
    if args['gene']:
        return jsonify(geo().get_results_by_gene(args['gene']))
    return jsonify(geo().get_analysis_results(
        args['study'],
        args['virus'],
        args['bto_id'],
        args['platform'],
        args['set']
    ))