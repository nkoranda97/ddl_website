from datetime import date
import os
import zipfile
import shutil

from flask import (
    Blueprint, flash, redirect, render_template, request, url_for, current_app, session
)
   
from werkzeug.utils import secure_filename

from flask_app.db import get_db
from .preprocess import preprocess

bp = Blueprint('analyze', __name__, url_prefix='/analyze')

@bp.route('/workspace/<int:project_id>')
def workspace(project_id):
    db = get_db()
    project = db.execute(
        'SELECT * FROM projects WHERE project_id = ?',
        (project_id,)
    ).fetchone()
    
    session['project_data'] = {
        'project_id': project['project_id'],
        'project_name': project['project_name'],
        'project_author': project['project_author'],
        'creation_date': project['creation_date'],
        'vdj_path': project['vdj_path'],
        'adata_path': project['adata_path']
    }

    
    return redirect('/dashapp')
