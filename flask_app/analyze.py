from datetime import date
import os
import zipfile
import shutil

from flask import (
    Blueprint, flash, redirect, render_template, request, url_for, current_app
)
   
from werkzeug.utils import secure_filename

from flask_app.db import get_db
from .preprocess import preprocess

bp = Blueprint('analyze', __name__, url_prefix='/analyze')

@bp.route('/workspace/<int:project_id>')
def workspace(project_id):
    db = get_db()
    project = db.execute(
        'SELECT project_name, project_author, creation_date, directory_path, vdj_path, adata_path FROM projects WHERE project_id = ?',
        (project_id,)
    ).fetchone()
    

    
    return render_template('analyze/workspace.html', project=project)
