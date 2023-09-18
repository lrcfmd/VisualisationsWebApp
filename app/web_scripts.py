from model import *
from visualise_square_web import *
from visualise_cube_web import *

#2d charged
model=Model()
model.setup_charged({'Li':1,'Al':3,'O':-2,'F':-1})
model.find_corners_edges(model.normal_vectors)
plotter=Square(model)
model.find_spreadout_points(5)
plotter.plot_plotting_df(columns='Composition')
plotter.show()
#2d uncharged
'''
model=Model()
model.setup_uncharged(['Li','O','F'])
model.find_corners_edges(model.normal_vectors)
plotter=Square(model)
model.find_spreadout_points(5)
plotter.plot_plotting_df(columns='Composition')
plotter.show()
'''
#3d charged
'''
model=Model()
model.setup_charged({'Li':1,'Al':3,'B':3,'O':-2,'F':-1})
model.find_corners_edges(model.normal_vectors)
plotter=Cube(model)
model.find_spreadout_points(5)
plotter.plot_plotting_df(columns='Composition')
plotter.show()
'''
#3d uncharged
'''
model=Model()
#find affine transformation to visualisation space
model.setup_uncharged(['Fe','Ni','Cu','Au'])
#find corners and edges of phase field
model.find_corners_edges(model.normal_vectors)
#create phase field visualisation
plotter=Cube(model)
#find spread out points
model.find_spreadout_points(5)
#plot the points the model found
plotter.plot_plotting_df(columns='Composition')
#show and not save fig
plotter.show()
#save fig as string_provided.html and dont show
plotter.show(save='test',show=False)
'''
