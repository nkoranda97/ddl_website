import os
from flask import Flask


def create_app(test_config = None):
    app = Flask(__name__, instance_relative_config = True)
    
    app.config.from_mapping(
        SECRET_KEY = '87b2beec9695cc7c2028e5ff6ce0127855741d4c0bc2ce5506dbfb4033302b91',
        DATABASE = os.path.join(app.instance_path, 'flask_app.sqlite'),
        UPLOAD_FOLDER = os.path.join(app.instance_path, 'uploads/')
    )

    
    if not os.path.exists(app.config['UPLOAD_FOLDER']):
        os.makedirs(app.config['UPLOAD_FOLDER'])

    if test_config is None:
        app.config.from_pyfile('config.py', silent=True),
    else:
        app.config.from_mapping(test_config)

    try:
        os.makedirs(app.instance_path)
    except OSError:
        pass


    
    from . import db
    db.init_app(app)
    
    from . import index
    app.register_blueprint(index.bp)
    
    from . import select
    app.register_blueprint(select.bp)
    
    from . import analyze
    app.register_blueprint(analyze.bp)
    

    
    @app.route('/static/<path:filename>')
    def serve_static(filename):
        return app.send_static_file(filename)
    
    
    return app
