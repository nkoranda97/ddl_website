from . import seqLogo, parser, colorMaps, alignment
from bokeh.io import show



color_maps = colorMaps.get_all_color_map_names()

def draw_seqlogo_from_file(df, color_scheme, web=False, plot_width=20, plot_height=160, steps=5, gene = 'all'):
    """ Gets a multiple sequence alignment file and draws the plot.

    :param file: The path to the clustal / clustal_num / msf / fasta file.
    :type file: str
    :param color_scheme: The name of the color scheme.
    :type color_scheme: str
    :param plot_width: An integer number for the width of the plot. Default is 20px.
    :type plot_width: int
    :param plot_height: An integer number for the height of the plot. Default is 160px.
    :type plot_height: int
    :param steps: An integer number to use as the step size of the 'x' axis. Default is 5.
    :type steps: int
    :return: A SeqLogo object.
    :rtype: SeqLogo
    """

    sl = seqLogo.SeqLogo(plot_width, plot_height, steps)
    parsed_sequences = parser.parse_df(df, gene = gene)
    sl.draw(parsed_sequences, color_scheme, web = False)
    return sl


def draw_alignment_from_file(df, color_scheme, web=False, plot_width=20, plot_height=160):
    """ Gets a multiple sequence alignment file and draws the plot.

    :param file: The path to the clustal / clustal_num / msf / fasta file.
    :type file: str
    :param color_scheme: The name of the color scheme.
    :type color_scheme: str
    :param plot_width: An integer number for the width of the plot. Default is 20px.
    :type plot_width: int
    :param plot_height: An integer number for the height of the plot. Default is 160px.
    :type plot_height: int
    :return: An Alignment object.
    :rtype: Alignment
    """
    aln = alignment.Alignment(plot_width, plot_height)
    parsed_sequences = parser.parse_df(df, 'all')
    if isinstance(parsed_sequences, Exception):
        return f'Drawing alignment failed with the following error: {parsed_sequences}'
    aln.draw(parsed_sequences, color_scheme)
    return aln

'''
def draw_dendrogram(file, web=False, plot_width=500, plot_height=500):
    """ Gets a dendrogram file and draws the plot.

    :param file: The path to the dendrogram file.
    :type file: str
    :param plot_width: An integer number for the width of the plot. Default is 500px.
    :type plot_width: int
    :param plot_height: An integer number for the height of the plot. Default is 500px.
    :type plot_height: int
    :return: A Dendrogram object.
    :rtype: Dendrogram
    """
    if file.split('.')[-1] not in file_formats:
        logging.error("File format must be clustal / clustal_num / msf / fasta!")
        return InvalidFileFormatException(f'{file} is not valid for this diagram type. Please use files with clustal msf or fasta extension.' )
    dend = dendrogram.Dendrogram(plot_width, plot_height)
    parsed_sequences = parser.parse_file(file)
    if isinstance(parsed_sequences, Exception):
        return f'Drawing dendrogram failed with the following error: {parsed_sequences}'
    dend.draw(parsed_sequences, web)
    return dend
    
'''

def get_all_color_map_names():
    """
    :return: List of all available color maps.
    :rtype: list
    """
    return colorMaps.get_all_color_map_names()


def get_nucleotide_color_map_names():
    """
    :return: List of all available color maps for nucleotides.
    :rtype: list
    """
    return colorMaps.get_nucleotide_color_map_names()


def get_protein_color_map_names():
    """
    :return: List of all available color maps for proteins.
    :rtype: list
    """
    return colorMaps.get_protein_color_map_names()

def show(logo):
    logo.show()