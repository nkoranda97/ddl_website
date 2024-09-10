import os

from flask import Flask, render_template

def create_app(test_config = None):
    app = Flask(__name__, instance_relative_config = True)
    
    app.config.from_mapping(
        SECRET_KEY = 'dev',
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

    @app.route('/')
    def index():
        return render_template('index/index.html')
    
    from . import db
    db.init_app(app)
    
    
    from . import select
    app.register_blueprint(select.bp)
    
    return app
