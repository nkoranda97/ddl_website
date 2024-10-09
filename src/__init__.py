import os
from flask import Flask, send_from_directory
from dotenv import load_dotenv
from flask_wtf.csrf import CSRFProtect


def create_app(test_config = None):
    app = Flask(__name__, instance_relative_config = True)
    load_dotenv(override=True)
    app.config.from_mapping(
        SECRET_KEY=os.getenv('SECRET_KEY'),
        DATABASE=os.path.join(app.instance_path, 'flask_app.sqlite'),
        UPLOAD_FOLDER=os.path.join(app.instance_path, 'uploads/')
    )

    csrf = CSRFProtect(app)

    
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
    
    from . import select_project
    app.register_blueprint(select_project.bp)
    
    from . import analyze
    app.register_blueprint(analyze.bp)
    
    @app.route('/favicon.ico')
    def favicon():
        return send_from_directory(os.path.join(app.root_path, 'static/images'),
                                   'favicon.ico', mimetype='image/vnd.microsoft.icon')

    
    @app.route('/static/<path:filename>')
    def serve_static(filename):
        return app.send_static_file(filename)
    
    
    return app
