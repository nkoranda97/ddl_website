from ..logo_generator import generate_logo
from bokeh.io import show

if __name__ == "__main__":
    import pandas as pd
    
    df = pd.read_csv('test_file.csv')
    
    p = generate_logo(df, 'seqlogo', chain = 'H', color='proteinClustal', width=80)
    show(p)