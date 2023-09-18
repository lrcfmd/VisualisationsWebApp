from flask_wtf import FlaskForm
from wtforms import StringField, SubmitField, IntegerField, Form, FloatField, FieldList, FieldList, FormField, HiddenField
from wtforms.validators import InputRequired

class ElementForm(Form):
    element = FloatField(0, [InputRequired()])
    charge = IntegerField()


class SearchForm(FlaskForm):

    message_field = HiddenField()
    output_data = HiddenField()
    n_points = IntegerField()
    elements = FieldList(FormField(ElementForm), min_entries=2)
    submit = SubmitField("Calculate")
    submit