from bokeh.models import ColumnDataSource, HoverTool, CustomJS, Dropdown, Slider, ColorPicker, TextInput, DataTable, TableColumn, Div
from bokeh.plotting import figure, show
from bokeh.layouts import column, row
from bokeh.palettes import Category20
from bokeh.transform import cumsum
import pandas as pd
from math import pi
from collections import Counter
import json

def bar(df, x, url):
    df[x] = df[x].astype(str)

    # Assign the same URL to all entries

    # Group data and sort
    group = df.groupby(x).size().reset_index(name='counts')
    group = group.sort_values(by='counts', ascending=False)
    
    x_agg = df.groupby(x).agg({col: list for col in df.columns if col != x})
    def format_dict(lst):
        counter = Counter(lst)
        sorted_items = sorted(counter.items(), key=lambda item: item[1], reverse=True)
        return '\n'.join([f"{k}: {v}" for k, v in sorted_items])
    
    x_agg = x_agg.applymap(format_dict)
    x_agg.reset_index(inplace = True)
    group = group.merge(x_agg, on=x)

    # Add URLs to the data source
    source = ColumnDataSource(group)

    # Create the plot
    p = figure(y_range=group[x].astype(str), title='Graph Title',
               toolbar_location="right", tools="tap,pan,wheel_zoom,box_zoom,reset,save",
               height=600, width=1200)
    p.hbar(y=x, right='counts', height=0.9, source=source)

    # Add hover tool
    hover = HoverTool()
    hover.tooltips = [("Category", f"@{x}"), ("Count", "@counts")]
    p.add_tools(hover)

    p.ygrid.grid_line_color = None
    p.x_range.start = 0


    '''
    # Add JavaScript callback for hyperlinking
    callback = CustomJS(args=dict(source=source), code="""
        var selected = source.selected.indices;
        var index = selected[0];
        var url = source.data['urls'][index];
        window.open(url, "_blank");
 
    """)
    '''
    callback = CustomJS(args=dict(source=source), code="""
        var selected = source.selected.indices;
        if (selected.length > 0) {
            var index = selected[0];
            var gene = source.data['""" + x + """'][index];
            var count = source.data['counts'][index];
            var content = "<h3>Gene: " + gene + "</h3>" +
                        "<p>Count: " + count + "</p>";
                        
            var options = ['j_call_VDJ', 'v_call_VDJ', 'c_call_VDJ', 'v_call_VJ', 'j_call_VJ', 'c_call_VJ'];

            options.forEach(function(option) {
                if (option !== '""" + x + """') {
                    var data = source.data[option][index];
                    content += "<p>" + option + ":</p><pre>" + data + "</pre><br>";
                }
            });

            openPanel(content);
        }
    """)

    
    p.js_on_event('tap', callback)

    toolbox = row()

    # Title input
    title_input = TextInput(value='Graph Title', title="Graph Title:")
    title_input_callback = CustomJS(args=dict(p=p), code="""
        p.title.text = cb_obj.value;
    """)
    title_input.js_on_change('value', title_input_callback)
    toolbox.children.append(title_input)

    # Color picker for bars
    color_picker = ColorPicker(title="Bar Color", color="#31AADE")
    color_picker_callback = CustomJS(args=dict(renderer=p.renderers[0]), code="""
        renderer.glyph.fill_color = cb_obj.color;
        renderer.glyph.line_color = cb_obj.color;
    """)
    color_picker.js_on_change('color', color_picker_callback)
    toolbox.children.append(color_picker)

    # Dropdown for bar width
    bar_width_dropdown = Dropdown(label="Bar Width", menu=[("0.5", "0.5"), ("0.7", "0.7"), ("0.9", "0.9")])
    bar_width_callback = CustomJS(args=dict(renderer=p.renderers[0]), code="""
        renderer.glyph.height = parseFloat(cb_obj.item);
    """)
    bar_width_dropdown.js_on_event('menu_item_click', bar_width_callback)
    toolbox.children.append(bar_width_dropdown)

    layout = column(p, toolbox, sizing_mode='stretch_both')

    return layout


def pie(df, x, title=''):
    df_copy = df.copy()
    df_copy = df_copy.dropna(subset=[x])
    df_copy.loc[:, x] = df_copy[x].astype(str)
    
    group = df_copy.groupby(x).size().reset_index(name='counts')
    group = group.sort_values(by='counts', ascending=False)
    group.reset_index(drop=True, inplace=True)
    
    group['angle'] = group['counts'] / group['counts'].sum() * 2 * pi
    group['percentage'] = group['counts'] / group['counts'].sum() * 100
    colors = (Category20[20] * (len(group) // 20 + 1))[:len(group)]
    group['color'] = colors
    
    source = ColumnDataSource(group)
    
    p = figure(height=600, width = 1000, title=title,
               toolbar_location="above", tools="pan,wheel_zoom,box_zoom,reset",
               x_range=(-0.5, 1.0))

    p.wedge(x=0, y=1, radius=0.4,
            start_angle=cumsum('angle', include_zero=True), end_angle=cumsum('angle'),
            line_color="white", fill_color='color', legend_field=x, source=source)
    
    p.axis.axis_label = None
    p.axis.visible = False
    p.grid.grid_line_color = None
    
    hover = HoverTool()
    hover.tooltips = [
        ("Category", f"@{x}"),
        ("Count", "@counts"),
        ("Percentage", "@percentage{0.2f}%")
    ]
    p.add_tools(hover)
    
    slider = Slider(start=0, end=50, value=0, step=0.01, title="Percent Cutoff")
    
    with open('src/static/JavaScript/update_pie.js', 'r') as js_file:
        js_code = js_file.read()
    
    callback = CustomJS(args=dict(source=source, slider=slider, original_data=source.data, x=x), code=js_code)
    
    slider.js_on_change('value', callback)
    
    
    
    layout = column(slider, p, sizing_mode='stretch_both')
    
    return layout

from bokeh.models import ColumnDataSource, DataTable, TableColumn, HTMLTemplateFormatter
from bokeh.layouts import column
from bokeh.io import show
import pandas as pd

def table(df):
    Columns = [
        TableColumn(
            field=Ci, 
            title=Ci, 
            formatter=HTMLTemplateFormatter(template=f'<div title="<%= value %>"><%= value %></div>')
        ) 
        for Ci in df.columns
    ]  # bokeh columns with tooltips
    return DataTable(columns=Columns, source=ColumnDataSource(df), sizing_mode='stretch_both')  # bokeh table


 
if __name__ == "__main__":
    df = pd.read_csv('/Users/nick/dandelion_tutorial/ddl_test.csv')
    plot = bar(df, 'v_call_VDJ')
    
    show(plot)
    