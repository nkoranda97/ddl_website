import math
from . import colorMaps, utils
from bokeh.io import show
from bokeh.models import ColumnDataSource, Plot, LinearAxis, Grid, Range1d
from bokeh.models.glyphs import ImageURL
from bokeh.layouts import column
from datetime import datetime
from flask import url_for

base_path = utils.get_base_path()
print(f"Base path: {base_path}")  # Debugging: Print the base path

class SeqLogo(object):
    
    def __init__(self, plot_width, plot_height, steps):
        self.plots = []
        self.plot_width = plot_width
        self.plot_height = plot_height
        # Steps of the ticker on the x axis
        self.steps = steps

    def draw(self, parsed_sequences, color_scheme, web=True):
        ratio = utils.calculate_ratio(parsed_sequences)
        seq_length = len(ratio)
        # data numbering starts at 1 -> so should the ratio
        ratio.insert(0, {})

        # Each subplot should be plot_width letter 'long' at max.
        subplot_count = math.ceil(seq_length / self.plot_width)
        # Create subplots.
        for k in range(0, subplot_count):
            x_start = 1 + self.plot_width * k

            if k == subplot_count - 1:
                x_end = seq_length + 1
                self.plot_width = (x_end - x_start)
            else:
                x_end = self.plot_width + 1 + self.plot_width * k

            # x and y have the values of the axes of the plot - same for every sequence
            x = list(range(x_start, x_end))  # Convert range to list
            y = [0] * (x_end - x_start)
            source_seq = ColumnDataSource(dict(x=x, y=y))

            subplot = Plot(title=None, width=int(80 * self.plot_width), height=self.plot_height,
                           x_range=Range1d(start=x_start - 1, end=x_end),
                           y_range=Range1d(start=0, end=1),
                           min_border_top=10, toolbar_location=None,
                           )

            _sum = 0
            for position in range(x_start, x_end):
                for letter, rat in ratio[position].items():                      
                    remainder = 0
                    if not rat == 1.0:
                        remainder = _sum % 1
                        
                    _sum += rat
                    if not rat == 0.0:
                        
                        image_url = url_for('static', filename=f'/images/{color_scheme}_{letter}.svg', _external=True)
                        image = ImageURL(name=letter, url=dict(value=image_url), x=position, y=remainder,
                                            w=1, h=rat, anchor="bottom_center")
                        subplot.add_glyph(source_seq, image)

                _sum = 0

            xaxis = LinearAxis(axis_label="Position")
            xaxis.bounds = (x_start, x_end)
            xaxis.ticker = [x_start] + list(range(x_start + self.steps - 1, x_end, self.steps))  # Convert range to list
            subplot.add_layout(xaxis, 'below')

            subplot.add_layout(Grid(dimension=0, ticker=xaxis.ticker))
            self.plots.append(subplot)      
        return self.plots

    def show(self):
        return column(self.plots)