from bokeh.models import ColumnDataSource
from bokeh.palettes import Turbo256
from bokeh.plotting import figure
from bokeh.transform import factor_cmap
import pandas as pd

def bar(df, x, title=''):
    
    df[x] = df[x].astype(str)
    group = df.groupby(x).size().reset_index(name='counts')
    group = group.sort_values(by='counts', ascending=False)

    x_cmap = factor_cmap(x, palette=Turbo256, factors=sorted(df[x].unique()))
    
    # Create a ColumnDataSource
    source = ColumnDataSource(group)
    
    # Create a horizontal bar plot
    p = figure(y_range=group[x].astype(str), title=title,
               toolbar_location=None, tools="")
    p.hbar(y=x, right='counts', height=0.9, source=source, line_color=x_cmap, fill_color=x_cmap)
    
    p.ygrid.grid_line_color = None
    p.x_range.start = 0

    return p

# Example usage
if __name__ == "__main__":
    # Create a sample DataFrame
    
    df = pd.read_csv('/Users/nick/dandelion_tutorial/ddl_test.csv')
    
    # Generate the bar plot
    plot = bar(df,'v_call_VDJ')
    
    # To show the plot in a standalone script, you can still call show(plot)
    from bokeh.io import show
    show(plot)