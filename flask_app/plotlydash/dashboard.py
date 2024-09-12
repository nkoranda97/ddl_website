import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output, State
import plotly.express as px
import pandas as pd
from flask import session
import dandelion as ddl
import scanpy as sc

# Global dictionary to store preprocessed data
preprocessed_data_cache = {}

def init_dashboard(server):
    """Create a Plotly Dash dashboard."""
    dash_app = dash.Dash(
        server=server,
        routes_pathname_prefix="/dashapp/",
        external_stylesheets=[
            "/static/styles.css"
        ],
    )

    def load_project_data():
        """Load project data from session."""
        if 'project_data' in session:
            project_data = session['project_data']
            vdj_path = project_data['vdj_path']
            adata_path = project_data['adata_path']
            vdj = ddl.read_h5ddl(vdj_path)
            adata = sc.read(adata_path)
            return project_data['project_name'], vdj, adata
        return None, None, None

    def preprocess_data(vdj, adata):
        """Preprocess data."""
        vdj, adata = ddl.pp.check_contigs(vdj, adata)
        ddl.tl.find_clones(vdj)
        ddl.tl.generate_network(vdj)
        ddl.tl.clone_size(vdj)
        df = vdj.metadata
        df = df[df.isotype_status != 'Multi']
        return df

    def create_figure(df, chart_type):
        """Create figure based on chart type."""
        if chart_type == 'Bar':
            fig = bar_graph(df, x='v_call_VDJ', color='v_call_VDJ', sort=True, title='V Call VDJ Usage')
            fig.update_layout(showlegend=False)
            fig.update_traces(dict(marker_line_width=0))
        elif chart_type == 'Pie':
            fig = px.pie(df, names='v_call_VDJ', title='V Call VDJ Usage')
            fig.update_traces(textposition='inside', textinfo='percent+label')
        return fig

    @dash_app.callback(
        Output('project-info', 'children'),
        [Input('url', 'pathname')]
    )
    def display_project_info(pathname):
        project_name, vdj, adata = load_project_data()
        if vdj and adata:
            if project_name not in preprocessed_data_cache:
                df = preprocess_data(vdj, adata)
                preprocessed_data_cache[project_name] = df
            else:
                df = preprocessed_data_cache[project_name]

            fig_1 = create_figure(df, 'Bar')
            fig_2 = bar_graph(df, x='isotype_status', color='locus_status', title='Isotype Usage')
            fig_2.update_traces(dict(marker_line_width=0))

            return html.Div([
                html.A("Back to Projects", href="/project_list"),
                html.H3(f"Project Name: {project_name}"),
                dcc.Dropdown(
                    id='chart-type',
                    options=[
                        {'label': 'Bar Chart', 'value': 'Bar'},
                        {'label': 'Pie Chart', 'value': 'Pie'}
                    ],
                    value='Bar',
                    clearable=False
                ),
                dcc.Graph(id='figure-1', figure=fig_1),
                dcc.Graph(figure=fig_2),
                html.Div(id='chart-type-store', style={'display': 'none'})
            ])
        return html.Div("No project selected.")

    @dash_app.callback(
        Output('figure-1', 'figure'),
        [Input('chart-type', 'value')],
        [State('url', 'pathname')]
    )
    def update_chart(chart_type, pathname):
        project_name, _, _ = load_project_data()
        if project_name and project_name in preprocessed_data_cache:
            df = preprocessed_data_cache[project_name]
            return create_figure(df, chart_type)
        return {}

    dash_app.layout = html.Div([
        dcc.Location(id='url', refresh=False),
        html.Div(id='project-info')
    ])

    return dash_app.server

def bar_graph(df, x, color, title, sort=True):
    if sort:
        counts = df[x].value_counts().reset_index()
        counts.columns = [x, 'count']
        sorted_categories = counts.sort_values(by='count', ascending=False)[x].tolist()
        return px.bar(df, x=x, color=color, title=title, category_orders={x: sorted_categories})
    return px.bar(df, x=x, color=color, title=title)
