from flask_wtf import FlaskForm
from wtforms import StringField, SubmitField, IntegerField, Form, FloatField, FieldList, FieldList, FormField, HiddenField
from wtforms.validators import InputRequired, Optional

class ElementForm(Form):
    element = StringField("Element Symbol", validators=[InputRequired()])
    charge = IntegerField("Charge Constraint",validators=[Optional()])

class KnownForm(Form):
    known = StringField("Element Symbol", validators=[Optional()])

class PrecursorForm(Form):
    precursor = StringField("Element Symbol", validators=[Optional()])

class SearchForm(FlaskForm):
    message_field = HiddenField()
    output_data = HiddenField()
    n_points = IntegerField(validators=[Optional()])
    elements = FieldList(FormField(ElementForm), min_entries=3)
    knowns=FieldList(FormField(KnownForm),min_entries=1)
    precursors=FieldList(FormField(PrecursorForm),min_entries=1)
    submit = SubmitField("Calculate")
    submit
