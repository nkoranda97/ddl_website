from . import colorMaps
from bokeh.models import ColumnDataSource, Plot, LinearAxis, Grid, Range1d, HoverTool, ZoomInTool, ZoomOutTool, PanTool
from bokeh.models.glyphs import Text
from bokeh.io import show, output_file, save
from bokeh.io.export import export_png, export_svg
from bokeh.transform import factor_cmap
from bokeh.layouts import column
import math
from datetime import datetime

class Alignment(object):
    def __init__(self, plot_width, plot_height):
        self.plots = []
        self.plot_width = plot_width
        self.plot_height = plot_height

    def draw(self, parsed_sequences, color_scheme):
        color_map = colorMaps.get_colormap(color_scheme)
        sequence_count = len(parsed_sequences)
        seq_lengths = [parsed_sequences[i].get('seq_length') for i in range(0, sequence_count)]
        seq_names = [parsed_sequences[i].get('id') for i in range(0, sequence_count)]
        max_seq_length = max(seq_lengths)

        # Each subplot should be plot_width letter 'long' at max.
        subplot_count = math.ceil(max_seq_length / self.plot_width)

        # Create subplots.
        for k in range(0, subplot_count):
            x_start = 1 + self.plot_width * k

            if k == subplot_count - 1:
                x_end = max_seq_length + 1
                self.plot_width = (x_end - x_start) * 1.25
            else:
                x_end = self.plot_width + 1 + self.plot_width * k

            # X has the values of the X axis of the plot - same for every sequence
            # Y has the values of the Y axis of the plot - different for every sequence
            x = list(range(x_start, x_end))

            subplot = Plot(title=None, width=int(10 * self.plot_width), height=30 * sequence_count,
                           x_range=Range1d(start=x_start - 1, end=x_end),
                           y_range=Range1d(start=0, end=sequence_count),
                           min_border_top=10, min_border_left=50,
                           toolbar_location='below',
                           tools=[HoverTool(tooltips=[("Position", "@x"), ("Sequence", "@y")]), ZoomInTool(), ZoomOutTool(), PanTool()]
                           )

            # Add sequences to the plot.
            for i in range(0, sequence_count):
                y_seq = []
                # The y value for the i. sequence will be i for all letters.
                for j in range(x_start, x_end):
                    y_seq.append(i + 1)

                seq = list(parsed_sequences[i].get('seq'))[x_start - 1:x_end - 1]
                source_seq = ColumnDataSource(dict(x=x, y=y_seq, text=seq))

                glyph_seq = Text(x="x", y="y", text="text",
                                 text_color=factor_cmap('text', palette=list(color_map.values()),
                                                        factors=list(color_map.keys())),
                                 text_font_size="9pt",
                                 x_offset=-3.3,
                                 text_line_height=0.8,
                                 text_baseline="top")
                subplot.add_glyph(source_seq, glyph_seq)

            yaxis = LinearAxis(axis_label="Sequence name")
            yaxis.bounds = (1, sequence_count)
            yaxis.ticker = list(range(1, sequence_count + 1))
            label_dict = {i + 1: seq_names[i] for i in range(0, sequence_count)}
            yaxis.major_label_overrides = label_dict
            subplot.add_layout(yaxis, 'left')

            xaxis = LinearAxis(axis_label="Position")
            xaxis.bounds = (x_start, x_end)
            xaxis.ticker = [x_start] + list(range(x_start + 19, x_end, 20))
            subplot.add_layout(xaxis, 'below')

            subplot.add_layout(Grid(dimension=0, ticker=xaxis.ticker))
            self.plots.append(subplot)

        return self.plots

    def show(self):
        return column(self.plots)
