from .bioviz.msa import draw_alignment_from_file, draw_seqlogo_from_file
from bokeh.layouts import column

def generate_logo(df, logo_type, color, width=None, gene = 'any'):
    if logo_type == 'alignment':
        logo = draw_alignment_from_file(df, color, plot_width=width if width is not None else 100, gene = gene)
    elif logo_type == 'seqlogo':
        logo = draw_seqlogo_from_file(df, color, plot_width=width if width is not None else 100, gene = gene)

    return logo.show()
