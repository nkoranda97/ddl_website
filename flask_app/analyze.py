from flask import Blueprint, render_template, session
from flask_app.db import get_db
from bokeh.embed import components
from bokeh.plotting import figure
from bokeh.resources import CDN
import dandelion as ddl
import scanpy as sc
from .plotting import plot

bp = Blueprint('analyze', __name__, url_prefix='/analyze')

@bp.route('/workspace/<int:project_id>')
def workspace(project_id):
    db = get_db()
    project = db.execute(
        'SELECT * FROM projects WHERE project_id = ?',
        (project_id,)
    ).fetchone()
    
    vdj, adata = load_project(project)
    
    df = vdj.metadata
    
    
    # Create a simple scatter plot
    p1 = plot.bar(df, x = 'v_call_VDJ')
    

    # Create a simple bar chart
    p2 = plot.pie(df, x = 'v_call_VDJ')

    # Create a simple line chart
    p3 = figure(title="Line Chart")
    p3.line([1, 2, 3], [4, 5, 6])
    

    # Get the script and div components
    script, div = components([p1, p2, p3])

    return render_template('analyze/workspace.html',
                           script=script, 
                           div=div, 
                           project=project,
                           resources=CDN.render())

def load_project(project):
    vdj_path = project['vdj_path']
    adata_path = project['adata_path']
    vdj = ddl.read_h5ddl(vdj_path)
    if adata_path != 'NULL':
        adata = sc.read(adata_path)
        return vdj, adata
    else:
        return vdj, None