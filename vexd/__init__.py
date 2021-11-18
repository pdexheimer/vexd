import os

from flask import Flask, request
from markupsafe import Markup

from . import api, geo, ui


def create_app(test_config=None):
    # create and configure the app
    app = Flask(__name__, instance_relative_config=True)
    app.config.from_mapping(
        SECRET_KEY='dev',
        MONGODB_HOST='localhost',
        MONGODB_PORT=27017,
        MONGODB_USER=None,
        MONGODB_PASS=None,
        ALWAYS_REPORT_HTTPS=False,
    )

    if test_config is None:
        # load the instance config, if it exists, when not testing
        app.config.from_pyfile('config.py', silent=True)
    else:
        # load the test config if passed in
        app.config.from_mapping(test_config)

    geo.init(app.config['MONGODB_HOST'],
        app.config['MONGODB_PORT'],
        app.config['MONGODB_USER'],
        app.config['MONGODB_PASS'])

    # ensure the instance folder exists
    try:
        os.makedirs(app.instance_path)
    except OSError:
        pass
    
    app.register_blueprint(api.bp)
    app.register_blueprint(ui.bp)

    @app.before_request
    def fake_https_if_asked():
        if app.config['ALWAYS_REPORT_HTTPS']:
            request.scheme = 'https'

    @app.template_filter('pluralize')
    def pluralize(word, count, plural_word=None):
        if count == 1:
            return word
        if plural_word is None:
            return word+'s'
        return plural_word

    @app.template_filter('virus_class')
    def get_sample_category(sample):
        key = 'normalized_virus' if sample['valid_experiment'] else 'invalid_reason'
        return sample[key]
    
    @app.template_filter()
    def highlight(value, is_signif):
        template = Markup("<b>{}</b>") if is_signif else "{}"
        return template.format(value)
    
    @app.template_test('virus_class_is')
    def sample_category_equals(sample, value):
        return get_sample_category(sample) == value

    return app
