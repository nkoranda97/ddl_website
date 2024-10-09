import os
import zipfile
import shutil

from flask import flash, redirect
   
from werkzeug.utils import secure_filename  
        
        
def file_extraction(files, project_folder, landing_page):  
    sample_folders = []
    for file in files:
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
                return redirect(landing_page)
        except zipfile.BadZipFile:
            flash(f'Failed to unzip file {filename}.')
            return redirect(landing_page)
        except Exception as e:
            flash(f'An error occurred while processing file {filename}.')
            return redirect(landing_page)
        
    return sample_folders
        