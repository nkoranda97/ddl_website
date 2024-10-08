from flask_wtf import FlaskForm
from flask_wtf.file import FileRequired, FileAllowed, MultipleFileField
from wtforms import StringField, RadioField, SelectField
from wtforms.validators import DataRequired


class LoginForm(FlaskForm):
    username = StringField('Name', validators=[DataRequired()])
    password = StringField('Password', validators=[DataRequired()])
    
class DataUpload(FlaskForm):
    project_name = StringField('Project Name', validators=[DataRequired()])
    author_name = StringField('Author', validators=[DataRequired()])
    data_uploaded = RadioField('Data Uploaded', choices=[('VDJ', 'VDJ'), ('GEX', 'GEX'), ('Both', 'Both')], validators=[DataRequired()])
    species = RadioField('Species', choices=[('human', 'Human'), ('mouse', 'Mouse')], validators=[DataRequired()])
    zip_folder = MultipleFileField('Zip Folder', validators=[FileRequired(),
                                                     FileAllowed(['zip'], 'Zip File Upload Only')])

class GeneSelect(FlaskForm):
    gene = SelectField('Gene', validators=[DataRequired()])
    
    def __init__(self, genes, *args, **kwargs):
        super(GeneSelect, self).__init__(*args, **kwargs)
        self.genes = genes
        self.gene.choices = [(gene, gene) for gene in genes]