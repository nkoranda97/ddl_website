from datetime import date
import os
import shutil

from flask import (
    Blueprint, flash, redirect, render_template, request, url_for, current_app
)
   
from werkzeug.utils import secure_filename

from src.db import get_db
from .utils.ddl import preprocess
from .index import login_required
from .utils.forms import DataUpload
from .utils.file_handle import file_extraction

bp = Blueprint('select', __name__, url_prefix = '/select')

@bp.route('/upload', methods = ('GET', 'POST'))
@login_required
def upload():
    form = DataUpload()
    
    if form.validate_on_submit():
        project_name = form.project_name.data
        author_name = form.author_name.data
        data_uploaded = form.data_uploaded.data
        species = form.species.data
        current_date = date.today()
        files = form.zip_folder.data
        db = get_db()
        project_folder = os.path.join(current_app.config['UPLOAD_FOLDER'], project_name)
        if not os.path.exists(project_folder):
            os.makedirs(project_folder)
        
        sample_folders = file_extraction(files=files, project_folder=project_folder,landing_page=request.url)
        
        error = None
        
        if data_uploaded == 'Both':
            
            adata_path, vdj_path = preprocess(project_folder, sample_folders, data_uploaded, species)

            if error is None:
                try:
                    db.execute(
                        'INSERT INTO projects (project_name, project_author, creation_date, directory_path, vdj_path, adata_path) VALUES (?, ?, ?, ?, ?, ?)',
                        (project_name, author_name, current_date, project_folder, vdj_path, adata_path)
                    )
                    db.commit()
                except db.IntegrityError:
                    error = f'Project Name {project_name} already exists.'
                else:
                    return redirect(url_for('select.project_list'))
        elif data_uploaded == 'VDJ':
            vdj_path = preprocess(project_folder, sample_folders, data_uploaded, species)

            if error is None:
                try:
                    db.execute(
                        'INSERT INTO projects (project_name, project_author, creation_date, directory_path, vdj_path, adata_path) VALUES (?, ?, ?, ?, ?, ?)',
                        (project_name, author_name, current_date, project_folder, vdj_path, 'NULL')
                    )
                    db.commit()
                except db.IntegrityError:
                    error = f'Project Name {project_name} already exists.'
                else:
                    return redirect(url_for('select.project_list'))

        
    return render_template('select/upload.html', form=form)

@bp.route('/project_list')
@login_required
def project_list():
    db = get_db()
    projects = db.execute(
        'SELECT project_id, project_name, project_author, creation_date FROM projects'
    ).fetchall()
    return render_template('select/project_list.html', projects = projects)

@bp.route('/delete_project/<int:project_id>', methods=['POST'])
@login_required
def delete_project(project_id):
    db = get_db()
    project = db.execute('SELECT directory_path FROM projects WHERE project_id = ?', (project_id,)).fetchone()
    
    if project is None:
        flash('Project not found.')
        return redirect(url_for('select.project_list'))
    
    directory_path = project['directory_path']
    
    if os.path.exists(directory_path):
        shutil.rmtree(directory_path)
    
    db.execute('DELETE FROM projects WHERE project_id = ?', (project_id,))
    db.commit()
    
    flash('Project deleted successfully.')
    return redirect(url_for('select.project_list'))

