from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment


# This parser requires a file and returns a Seq.
def parse_df(df, gene):

    sequences = []

    
    # Extract sequences based on the condition
    if gene and gene != 'all':
        sequences_dict = df[df['v_call'] == gene][['sequence_id', 'sequence_alignment']].to_dict()['sequence_alignment']
    else:
        sequences_dict = df[['sequence_id', 'sequence_alignment']].to_dict()['sequence_alignment']    
        
    
   
    alignment = [SeqRecord(Seq(value.replace('.', '')).translate(), id=str(key)) for key, value in sequences_dict.items()]

 
    max_length = max(len(record.seq) for record in alignment)
    
    alignment = [SeqRecord(Seq(str(record.seq).ljust(max_length, '-')), id=record.id) for record in alignment]


    for seq_record in alignment:
        sequences.append({
            'id': seq_record.id,
            'seq': str(seq_record.seq),
            'seq_length': len(seq_record)
        })
    return sequences

if __name__ == "__main__":
    import pandas as pd
    df = pd.read_csv('test_file.csv')
    
    print(parse_df(df, gene = 'all'))