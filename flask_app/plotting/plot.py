from bokeh.models import ColumnDataSource, HoverTool, CustomJS, Dropdown, Slider, ColorPicker, TextInput
from bokeh.plotting import figure, show
from bokeh.layouts import column
from bokeh.palettes import Category20
from bokeh.transform import cumsum
import pandas as pd
from math import pi

def bar(df, x):
    df[x] = df[x].astype(str)
    group = df.groupby(x).size().reset_index(name='counts')
    group = group.sort_values(by='counts', ascending=False)
    
    source = ColumnDataSource(group)
    
    # Convert DataFrame to dictionary
    df_dict = df.to_dict(orient='list')
    
    p = figure(y_range=group[x].astype(str), title = 'Graph Title',
               toolbar_location="right", tools="pan,wheel_zoom,box_zoom,reset,save",
               height = 600, width = 1200)
    hbar = p.hbar(y=x, right='counts', height=0.9, source=source)
    
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
    
    picker = ColorPicker(color = 'red')
    picker.js_link('color', hbar.glyph, 'fill_color')
    picker.js_link('color', hbar.glyph, 'line_color')
    
    
    color_picker_callback = CustomJS(args=dict(picker=picker), code="""
        picker.title = 'Color: ' + picker.color;
        picker.style.backgroundColor = picker.color;
    """)
    picker.js_on_change('color', color_picker_callback)
    
    text_input = TextInput(value='Graph Title', title="Graph Title:")
    text_input_callback = CustomJS(args=dict(p=p), code="""
        p.title.text = cb_obj.value;
    """)
    text_input.js_on_change('value', text_input_callback)
    
    return column(text_input,dropdown, p, picker)


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
    
    with open('flask_app/static/JavaScript/update_pie.js', 'r') as js_file:
        js_code = js_file.read()
    
    callback = CustomJS(args=dict(source=source, slider=slider, original_data=source.data, x=x), code=js_code)
    
    slider.js_on_change('value', callback)
    
    
    
    layout = column(slider, p)
    
    return layout
    
if __name__ == "__main__":
    df = pd.read_csv('/Users/nick/dandelion_tutorial/ddl_test.csv')
    plot = bar(df, 'v_call_VDJ')
    
    show(plot)
    