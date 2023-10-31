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
from app.validator import add_knowns, add_precursors
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


def render_3d(elements, n_points=None,knowns=None,precursors=None):

    model=Model()

    if isinstance(elements, dict):
        model.setup_charged(elements)
    else:
        model.setup_uncharged(elements)

    model.find_corners_edges(model.normal_vectors,add_corners_to_knowns=True)
    plotter=Cube(model)

    if precursors:
        model.add_precursors(precursors)
        plotter.plot_mesh(
            model.precursors_constrained,poly=True,name='Precursor accessible')

    if knowns:
        model.add_knowns(knowns,make_plotting_df=True)
        plotter.plot_plotting_df(['Label'],name='Known phases',c='red')

    if n_points is not None:
        num_steps=2500
        T=50
        if hasattr(model,'knowns_constrained'):
            model.find_spreadout_points(num_steps,n_points,
                                        make_plotting_df=True,
                                        use_cut_omega=bool(precursors),
                                        T=T)
            plotter.plot_plotting_df(
            plotting_columns='Composition',name='Suggested points',c='green')
    fig = plotter.show(return_fig=True)
    return  pio.to_html(fig,
                        full_html=False,default_width='100vw',
                        default_height="100vh")

def render_2d(elements, n_points=None,knowns=None,precursors=None):
    model=Model()
    if isinstance(elements, dict):
        model.setup_charged(elements)
    else:
        model.setup_uncharged(elements)
    model.find_corners_edges(model.normal_vectors,add_corners_to_knowns=True)
    plotter=Square(model)

    if precursors:
        model.add_precursors(precursors)
        plotter.plot_mesh(
            model.precursors_constrained,poly=True,name='Precursor accessible')

    if knowns:
        model.add_knowns(knowns,make_plotting_df=True)
        plotter.plot_plotting_df(['Label'],name='Known phases',c='red')

    if n_points is not None:
        num_steps=2500
        T=50
        model.find_spreadout_points(
            num_steps,n_points,make_plotting_df=True,
            use_cut_omega=bool(precursors),T=T)
    plotter.plot_plotting_df(
        plotting_columns=['Composition'],c='green',name='Suggested_points')
    fig = plotter.show(return_fig=True)
    return  pio.to_html(fig, full_html=False)


#Define route
@app.route("/", methods=['GET', 'POST'])
@app.route("/predict", methods=['GET', 'POST'])
def predict():
    form = SearchForm()
    if form.validate_on_submit():
        #add any known phases
        try:
            knowns=[]
            if form.knowns[0].known.data:
                for i in form.knowns:
                    if i.known.data:
                        try:
                            knowns.append(add_knowns(i.known.data))
                        except ValueError as e:
                            return render_template(
                                "MOF_ml.html", form=form,
                                message=(f"Error:\n Known phase {i.known.data} is not"
                                +" recognised")) 
            precursors=[]
            if form.precursors[0].precursor.data:
                for i in form.precursors:
                    if i.precursor.data:
                        try:
                            precursors.append(add_precursors(i.precursor.data))
                        except ValueError as e:
                            return render_template(
                                "MOF_ml.html", form=form,
                                message=(f"Error:\n Precursor {i.precursor.data} is not"
                                +" recognised")) 
            #Parse strings
            if form.elements[0].charge.data:
                elements = {x.element.data: x.charge.data for x in
                            form.elements if x.charge.data and x.element.data}
                app.logger.debug(elements)
                if len (elements) > 5 or len (elements) <3:
                    return render_template(
                        "MOF_ml.html", form=form,
                        message="Please enter between 3 and 5 elements")
                if len(elements) > 4:
                    return render_template(
                        "MOF_ml.html", form=form,results=render_3d(
                            elements,
                            n_points=form.n_points.data,
                            knowns=knowns,
                            precursors=precursors),
                        message="")
                if len(elements) > 3:
                    return render_template(
                        "MOF_ml.html", form=form,results=render_2d(
                            elements,n_points=form.n_points.data,
                            knowns=knowns,
                            precursors=precursors),
                        message="")
                return render_template(
                    "MOF_ml.html", form=form, results=render_2d(
                        list(elements.keys()),
                        n_points=form.n_points.data,
                        knowns=knowns,
                        precursors=precursors),
                    message="")
            elements = [x.element.data for x in form.elements if x.element.data]
            app.logger.debug(elements)
            if len (elements) > 4 or len (elements) < 3:
                return render_template(
                    "MOF_ml.html", form=form, 
                    message=("Without charge constraints between 3 and 4 "
                             +"elements are supported"))
            if len(elements) > 3:
                return render_template(
                    "MOF_ml.html", form=form, results=render_3d(
                        elements, n_points=form.n_points.data,knowns=knowns,
                        precursors=precursors ),
                    message="")
            return render_template(
                "MOF_ml.html", form=form, results=render_2d(
                    elements, n_points=form.n_points.data,knowns=knowns,
                    precursors=precursors),
                message="")

        except Exception as e:
            app.logger.debug(e)
            return render_template(
                "MOF_ml.html", form=form,
                message=("Failed to process input, check it is properly "
                         +"formatted.")
            )

    return render_template("MOF_ml.html", form=form, message="")

