from flask_wtf import FlaskForm
from wtforms import StringField, SubmitField, IntegerField, Form, FloatField, FieldList, FieldList, FormField, HiddenField
from wtforms.validators import InputRequired, Optional

class ElementForm(Form):
    element = StringField("Element Symbol", validators=[InputRequired()])
    charge = IntegerField("Charge Constraint",validators=[Optional()])

class SearchForm(FlaskForm):

    message_field = HiddenField()
    output_data = HiddenField()
    n_points = IntegerField(validators=[Optional()])
    elements = FieldList(FormField(ElementForm), min_entries=3)
    submit = SubmitField("Calculate")
    submit