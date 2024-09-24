from bokeh.models import ColumnDataSource, HoverTool, CustomJS, Dropdown
from bokeh.plotting import figure, show
from bokeh.layouts import column
import pandas as pd

def bar(df, x, title=''):
    df[x] = df[x].astype(str)
    group = df.groupby(x).size().reset_index(name='counts')
    group = group.sort_values(by='counts', ascending=False)
    
    source = ColumnDataSource(group)
    
    # Convert DataFrame to dictionary
    df_dict = df.to_dict(orient='list')
    
    p = figure(y_range=group[x].astype(str), title=title,
               toolbar_location="above", tools="pan,wheel_zoom,box_zoom,reset")
    p.hbar(y=x, right='counts', height=0.9, source=source)
    
    hover = HoverTool()
    hover.tooltips = [("Category", f"@{x}"), ("Count", "@counts")]
    p.add_tools(hover)
    
    p.ygrid.grid_line_color = None
    p.x_range.start = 0
    
    menu = [(col, col) for col in df.columns]
    dropdown = Dropdown(label=x, button_type="warning", menu=menu)
    
    with open('flask_app/static/JavaScript/callback.js', 'r') as f:
        js_code = f.read()
    
    callback = CustomJS(args=dict(source=source, p=p, df_dict=df_dict, dropdown=dropdown, hover=hover), code=js_code)
    
    dropdown.js_on_event('menu_item_click', callback)
    
    return column(dropdown, p)

if __name__ == "__main__":
    df = pd.read_csv('/Users/nick/dandelion_tutorial/ddl_test.csv')
    plot = bar(df, 'v_call_VDJ')
    show(plot)