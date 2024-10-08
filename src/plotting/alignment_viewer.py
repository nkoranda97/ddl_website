from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO
from Bio.Seq import Seq

import os

import numpy as np

from bokeh.plotting import figure
from bokeh.models import ColumnDataSource, Plot, Grid, Range1d, Title
from bokeh.models.glyphs import Text, Rect
from bokeh.layouts import gridplot

def create_alignmnet(df, gene: str):
    sequences = df[df['v_call'] == gene][['sequence_id', 'sequence_alignment']].to_dict()['sequence_alignment']
    alignment = [SeqRecord(Seq(value.replace('.', '')).translate(), id=os.path.basename(key)) for key, value in sequences.items()]
    return MultipleSeqAlignment(alignment)

def get_colors(seqs):
    """make colors for bases in sequence"""
    text = [i for s in list(seqs) for i in s]
    clrs = {
        'A': 'red', 'C': 'blue', 'D': 'purple', 'E': 'pink', 'F': 'green',
        'G': 'orange', 'H': 'cyan', 'I': 'yellow', 'K': 'magenta', 'L': 'lime',
        'M': 'teal', 'N': 'brown', 'P': 'olive', 'Q': 'navy', 'R': 'maroon',
        'S': 'gold', 'T': 'silver', 'V': 'coral', 'W': 'indigo', 'Y': 'violet',
        '.': 'white', '-': 'white'
    }
    colors = [clrs[i] for i in text]
    return colors

def make_alignment(aln, fontsize="9pt", plot_width=800, title = None):
    """Bokeh sequence alignment view"""

    #make sequence and id lists from the aln object
    seqs = [rec.seq for rec in (aln)]
    ids = [rec.id for rec in aln]    
    text = [i for s in list(seqs) for i in s]
    colors = get_colors(seqs)    
    N = len(seqs[0])
    S = len(seqs)    
    width = .4

    x = np.arange(1,N+1)
    y = np.arange(0,S,1)
    #creates a 2D grid of coords from the 1D arrays
    xx, yy = np.meshgrid(x, y)
    #flattens the arrays
    gx = xx.ravel()
    gy = yy.flatten()
    #use recty for rect coords with an offset
    recty = gy+.5
    h= 1/S
    #now we can create the ColumnDataSource with all the arrays
    source = ColumnDataSource(dict(x=gx, y=gy, recty=recty, text=text, colors=colors))
    plot_height = len(seqs)*15+50
    x_range = Range1d(0,N+1, bounds='auto')
    if N>100:
        viewlen=100
    else:
        viewlen=N
    #view_range is for the close up view
    view_range = (0,viewlen)
    tools="xpan, xwheel_zoom, reset, save"

    #entire sequence view (no text, with zoom)
    p = figure(title=None, width= plot_width, height=50,
               x_range=x_range, y_range=(0,S), tools=tools,
               min_border=0, toolbar_location='below')
    rects = Rect(x="x", y="recty",  width=1, height=1, fill_color="colors",
                 line_color=None, fill_alpha=0.6)
    p.add_glyph(source, rects)
    p.yaxis.visible = False
    p.grid.visible = False  

    #sequence text view with ability to scroll along x axis
    p1 = figure(title=None, width=plot_width, height=plot_height,
                x_range=view_range, y_range=ids, tools="pan,reset,wheel_zoom",
                min_border=0, toolbar_location='below')#, lod_factor=1)          
    glyph = Text(x="x", y="y", text="text", text_align='center',text_color="black",
                text_font="monospace",text_font_size=fontsize)
    rects = Rect(x="x", y="recty",  width=1, height=1, fill_color="colors",
                line_color=None, fill_alpha=0.4)
    p1.add_glyph(source, glyph,)
    p1.add_glyph(source, rects)

    p1.grid.visible = False
    p1.xaxis.major_label_text_font_style = "bold"
    p1.yaxis.minor_tick_line_width = 0
    p1.yaxis.major_tick_line_width = 0
    p.sizing_mode = 'stretch_both'

    p = gridplot([[p1]], toolbar_location='below', sizing_mode='stretch_both')
    return p

def view_alignment(df, gene:str, title):
    alignment = create_alignmnet(df, gene)
    return make_alignment(alignment, plot_width=2000, title = title)

if __name__ == "__main__":
    from bokeh.io import show

    aln = AlignIO.read('/Users/nick/ddl_website/flask_app/plotting/output.aln','fasta')
    p = make_alignment(aln, plot_width=900)
    show(p)