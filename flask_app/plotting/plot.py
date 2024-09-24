from bokeh.models import ColumnDataSource, HoverTool
from bokeh.palettes import Turbo256
from bokeh.plotting import figure, show
from bokeh.transform import factor_cmap
import pandas as pd

def bar(df, x, title=''):
    
    df[x] = df[x].astype(str)
    group = df.groupby(x).size().reset_index(name='counts')
    
    group = group.sort_values(by='counts', ascending=False)
    
    x_cmap = factor_cmap(x, palette=Turbo256, factors=group[x].tolist())
    
    source = ColumnDataSource(group)
    
    p = figure(y_range=group[x].astype(str), title=title,
               toolbar_location="above", tools="pan,wheel_zoom,box_zoom,reset")
    p.hbar(y=x, right='counts', height=0.9, source=source, line_color=x_cmap, fill_color=x_cmap)
    
    hover = HoverTool()
    hover.tooltips = [("Category", f"@{x}"), ("Count", "@counts")]
    p.add_tools(hover)
    
    p.ygrid.grid_line_color = None
    p.x_range.start = 0

    return p

if __name__ == "__main__":
    
    df = pd.read_csv('/Users/nick/dandelion_tutorial/ddl_test.csv')
    
    plot = bar(df,'v_call_VDJ')
    
    from bokeh.io import show
    show(plot)