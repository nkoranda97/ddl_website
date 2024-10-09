from flask import Blueprint, render_template, session
from src.db import get_db
from bokeh.embed import components
from bokeh.plotting import figure
from bokeh.models import Tabs, TabPanel
from bokeh.resources import CDN

from .plotting import plot, alignment_viewer
from .bokeh_logo.logo_generator import generate_logo
from .utils.ddl import load_project
from .index import login_required
from .utils.forms import GeneSelect

bp = Blueprint('analyze', __name__, url_prefix='/analyze')

@bp.route('/graphs/<int:project_id>')
@login_required
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
    
    p2_1 = plot.pie(df, x='v_call_VDJ')
    p2_2 = plot.pie(df, x='c_call_VDJ')
    p2_3 = plot.pie(df, x='j_call_VDJ')
    p2_4 = plot.pie(df, x='isotype')
    
    p2 = Tabs(tabs=[
        TabPanel(child=p2_1, title='v_call_VDJ'),
        TabPanel(child=p2_2, title='c_call_VDJ'),
        TabPanel(child=p2_3, title='j_call_VDJ'),
        TabPanel(child=p2_4, title='isotype'),
    ])
    
    p3 = plot.table(df)
    
    
    
    p1p2p3 = Tabs(tabs=[
        TabPanel(child=p1, title='Bar Graph'),
        TabPanel(child=p2, title='Pie Graph'),
        TabPanel(child=p3, title='Data Table'),
    ])
    
    # Get the script and div components
    script, div = components([p1p2p3])

    return render_template('analyze/graphs.html',
                           script=script, 
                           div=div, 
                           project=project,
                           project_id=project_id,
                           resources=CDN.render(),
                           active_tab = 'graphs')

@bp.route('/alignment/<int:project_id>')
@login_required
def alignment(project_id):
    db = get_db()
    project = db.execute(
        'SELECT * FROM projects WHERE project_id = ?',
        (project_id,)
    ).fetchone()
    
    vdj, adata = load_project(project)
    
    df = vdj.metadata
    
    p = alignment_viewer.view_alignment(vdj.data, 'IGKV3-4*01', title='IGKV3-4*01')
    
    script, div = components([p])
    
    return render_template('analyze/alignment.html',
                        script=script, 
                        div=div, 
                        project=project,
                        project_id=project_id,
                        resources=CDN.render(),
                        active_tab = 'alignments')
    
@bp.route('/logo/<int:project_id>', methods=['GET', 'POST'])
@login_required
def logo(project_id):
    db = get_db()
    project = db.execute(
        'SELECT * FROM projects WHERE project_id = ?',
        (project_id,)
    ).fetchone()
    
    script = None; div = None
    
    vdj, adata = load_project(project)
    genes = vdj.data['v_call'].unique()
    
    form = GeneSelect(genes)
    if form.validate_on_submit():
        gene = form.gene.data
        p = generate_logo(vdj.data, 'seqlogo', color='proteinClustal', width=16, gene=gene)

        script, div = components([p])

    
    return render_template('analyze/logo.html',
                    script=script, 
                    div=div, 
                    project=project,
                    project_id=project_id,
                    resources=CDN.render(),
                    active_tab = 'logo',
                    form = form)