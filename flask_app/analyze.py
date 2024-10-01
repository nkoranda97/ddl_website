from flask import Blueprint, render_template, session
from flask_app.db import get_db
from bokeh.embed import components
from bokeh.plotting import figure
from bokeh.models import Tabs, TabPanel
from bokeh.resources import CDN

from .plotting import plot, alignment_viewer
from .bokeh_logo.logo_generator import generate_logo
from .ddl import load_project

bp = Blueprint('analyze', __name__, url_prefix='/analyze')

@bp.route('/graphs/<int:project_id>')
def graphs(project_id):
    db = get_db()
    project = db.execute(
        'SELECT * FROM projects WHERE project_id = ?',
        (project_id,)
    ).fetchone()
    
    vdj, adata = load_project(project)
    
    df = vdj.metadata
    
    
    p1_1 = plot.bar(df, x='v_call_VDJ')
    p1_2 = plot.bar(df, x='c_call_VDJ')
    p1_3 = plot.bar(df, x='j_call_VDJ')
    p1_4 = plot.bar(df, x='isotype')

    p1 = Tabs(tabs=[
        TabPanel(child=p1_1, title='v_call_VDJ'),
        TabPanel(child=p1_2, title='c_call_VDJ'),
        TabPanel(child=p1_3, title='j_call_VDJ'),
        TabPanel(child=p1_4, title='isotype'),
    ])
    
    # Create a simple bar chart
    p2_1 = plot.pie(df, x = 'v_call_VDJ')
    p2_2 = plot.pie(df, x = 'c_call_VDJ')
    p2_3 = plot.pie(df, x = 'j_call_VDJ')
    p2_4 = plot.pie(df, x = 'isotype')
    
    p2 = Tabs(tabs = [
        TabPanel(child = p2_1, title = 'v_call_VDJ'),
        TabPanel(child = p2_2, title = 'c_call_VDJ'),
        TabPanel(child = p2_3, title = 'j_call_VDJ'),
        TabPanel(child = p2_4, title = 'isotype'),
    ])
    
    p3 = plot.table(df)
    
    p4 = alignment_viewer.view_alignment(vdj.data, 'IGKV3-4*01', title = 'IGKV3-4*01')
    
    p5 = generate_logo(vdj.data, 'seqlogo', color='proteinClustal', width=16, gene = 'IGKV3-4*01')
    
    p1p2p3 = Tabs(tabs = [
        TabPanel(child = p1, title = 'Bar Graph'),
        TabPanel(child = p2, title = 'Pie Graph'),
        TabPanel(child = p3, title = 'Date Table'),
        TabPanel(child = p4, title = 'Alignment Viewer' )
    ])

    # Create a simple line chart
    p3 = figure(title="Line Chart")
    p3.line([1, 2, 3], [4, 5, 6])
    

    # Get the script and div components
    script, div = components([p1p2p3, p5, p3])

    return render_template('analyze/graphs.html',
                           script=script, 
                           div=div, 
                           project=project,
                           resources=CDN.render())

