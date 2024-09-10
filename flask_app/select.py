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

bp = Blueprint('select', __name__, url_prefix = '/')

@bp.route('/upload', methods = ('GET', 'POST'))
def upload():
    if request.method == 'POST':
        project_name = request.form['project_name']
        author_name = request.form['author_name']
        project_folder = os.path.join(current_app.config['UPLOAD_FOLDER'], project_name)
        db = get_db()
        error = None 
        current_date = date.today()
        
        if not project_name:
            error = 'Project Name Required'
        
        if not author_name:
            error = 'Name Required'
            
        if not os.path.exists(project_folder):
            os.makedirs(project_folder)
            
        files = request.files.getlist('folders')
        sample_folders = []
        
        for file in files:
            if file and file.filename.endswith('.zip'):
                filename = secure_filename(file.filename)
                zip_path = os.path.join(project_folder, filename)
                file.save(zip_path)

                try:
                    with zipfile.ZipFile(zip_path, 'r') as zip_ref:
                        zip_ref.extractall(project_folder)
                    os.remove(zip_path) 

                    extracted_folder = filename.rsplit('.', 1)[0]
                    sample_folder = os.path.join(project_folder, extracted_folder)
                    sample_folders.append(sample_folder)

                    files_in_folder = os.listdir(sample_folder)
                    if not any(f.endswith('.csv') for f in files_in_folder) or not any(f.endswith('.fasta') for f in files_in_folder):
                        flash(f'Folder {extracted_folder} does not contain both .csv and .fasta files.')
                        return redirect(request.url)
                except zipfile.BadZipFile:
                    flash(f'Failed to unzip file {filename}.')
                    return redirect(request.url)
                except Exception as e:
                    flash(f'An error occurred while processing file {filename}.')
                    return redirect(request.url)
                
        adata_path, vdj_path = preprocess(project_folder, sample_folders)

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
        
    return render_template('select/upload.html')

@bp.route('/project_list')
def project_list():
    db = get_db()
    projects = db.execute(
        'SELECT project_id, project_name, project_author, creation_date FROM projects'
    ).fetchall()
    return render_template('select/project_list.html', projects = projects)

@bp.route('/delete_project/<int:project_id>', methods=['POST'])
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

