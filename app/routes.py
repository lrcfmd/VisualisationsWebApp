import os
import re
from flask_restful import Resource, Api, reqparse
import numpy as np
import joblib
import pandas as pd
from logging.config import dictConfig
from app.model import Model
from app.visualise_cube_web import Cube
from app.visualise_square_web import Square
from jinja2 import BaseLoader, TemplateNotFound,ChoiceLoader, FileSystemLoader
from urllib import request, parse
from flask import render_template, jsonify, make_response
from plotly import io as pio
from app import app, api
from app.forms import SearchForm

#We split loading of templates between local and remote as common templates
#For LMDS are hosted seperately to reduce code duplication
#Thus we define a quick class to use to fetch these remotely
class UrlLoader(BaseLoader):
    def __init__(self, url_prefix):
        self.url_prefix = url_prefix

    def get_source(self, environment, template):
        url = parse.urljoin(self.url_prefix, template)
        app.logger.debug(url)
        try:
            t = request.urlopen(url)
            if t.getcode() is None or t.getcode() == 200:
                return t.read().decode('utf-8'), None, None
        except IOError:
            pass
        raise TemplateNotFound(template)


#Configure logging
dictConfig({"version": 1,
            "disable_existing_loggers": False,
            "formatters": {"default": {
                        "format": '[%(asctime)s] %(levelname)s in %(module)s: %(message)s',
                }},

            "handlers": {
                "wsgi": {
                    "class": "logging.StreamHandler",
                    "stream": "ext://flask.logging.wsgi_errors_stream",
                    "formatter": "default",
                    }
                },

            "root": {"level": "DEBUG", "handlers": ["wsgi"]},
            })


app.jinja_loader = ChoiceLoader([app.jinja_loader, UrlLoader("https://lmds.liverpool.ac.uk/static")])

@app.route("/API_info", methods=['GET', 'POST'])
def api_info():
    return render_template("api_info.html")


def render_3d(elements, n_points=None):
    model=Model()
    if isinstance(elements, dict):
        model.setup_charged(elements)
    else:
        model.setup_uncharged(elements)
    model.find_corners_edges(model.normal_vectors)
    plotter=Cube(model)
    if n_points is not None:
        model.find_spreadout_points(n_points)
    else:
        model.find_spreadout_points(5)
    plotter.plot_plotting_df(columns='Composition')
    fig = plotter.show(return_fig=True)
    return  pio.to_html(fig, full_html=False)

def render_2d(elements, n_points=None):
    model=Model()
    if isinstance(elements, dict):
        model.setup_charged(elements)
    else:
        model.setup_uncharged(elements)
    model.find_corners_edges(model.normal_vectors)
    plotter=Square(model)
    if n_points is not None:
        model.find_spreadout_points(n_points)
    else:
        model.find_spreadout_points(5)
    plotter.plot_plotting_df(columns='Composition')
    fig = plotter.show(return_fig=True)
    return  pio.to_html(fig, full_html=False)


#Define route
@app.route("/", methods=['GET', 'POST'])
@app.route("/predict", methods=['GET', 'POST'])
def predict():
    form = SearchForm()
    if form.validate_on_submit():
        try:
            #Parse strings
            if form.elements[0].charge.data:
                elements = {x.element.data: x.charge.data for x in form.elements if x.charge.data and x.element.data}
                app.logger.debug(elements)
                if len (elements) > 5 or len (elements) <3:
                    return render_template("MOF_ml.html", form=form, message="Please enter between 3 and 5 elements") 
                if len(elements) > 4:
                    return render_template("MOF_ml.html", form=form, results=render_3d(elements, n_points=form.n_points.data), message="")
                if len(elements) > 3:
                    return render_template("MOF_ml.html", form=form, results=render_2d(elements, n_points=form.n_points.data), message="")
                return render_template("MOF_ml.html", form=form, results=render_2d(list(elements.keys()), n_points=form.n_points.data), message="")
            
            elements = [x.element.data for x in form.elements if x.element.data]
            app.logger.debug(elements)
            if len (elements) > 4 or len (elements) < 3:
                return render_template("MOF_ml.html", form=form, message="Without charge constraints between 3 and 4 elements are supported")
            if len(elements) > 3:
                return render_template("MOF_ml.html", form=form, results=render_3d(elements, n_points=form.n_points.data), message="")
            return render_template("MOF_ml.html", form=form, results=render_2d(elements, n_points=form.n_points.data), message="")
        
        except Exception as e:
            app.logger.debug(e)
            return render_template("MOF_ml.html", form=form, message="Failed to process input, check it is properly formatted.")
            

    return render_template("MOF_ml.html", form=form, message="")

