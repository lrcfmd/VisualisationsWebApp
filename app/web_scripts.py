from model import *
from visualise_square_web import *
from visualise_cube_web import *
from pymatgen.core import Composition

#2d charged
'''
model=Model()
model.setup_charged({'Li':1,'Al':3,'O':-2,'F':-1})
model.find_corners_edges(model.normal_vectors)
plotter=Square(model)
model.find_spreadout_points(5)
plotter.plot_plotting_df(columns='Composition')
plotter.show()
'''
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
model=Model()
model.setup_charged({'Li':1,'Al':3,'B':3,'O':-2,'F':-1})
precursors=["Li2O","LiAlO2","BF3","LiAlB2O4F1"]
precursors=[Composition(x) for x in precursors]
model.find_corners_edges(model.normal_vectors,add_corners_to_knowns=True)
plotter=Cube(model)
model.add_precursors(precursors)
plotter.plot_mesh(
    model.precursors_constrained,poly=True,name='Precursor accessible')
#plotter.plot_points(model.omega_cut,s=0.5,c='green')
#plotter.plot_points(model.omega,s=0.5,name='all')


#knowns=[Composition("LiAlBO3F"),Composition("Li4AlB2O5F5")]
#model.add_knowns(knowns,make_plotting_df=True)
#plotter.plot_plotting_df(
#    plotting_columns=['Label'],c='red',name="Known phases")
model.find_spreadout_points(
    250,15,T=50,make_plotting_df=True,use_cut_omega=True)
plotter.plot_plotting_df(
    plotting_columns=['Composition'],c='green',name="Suggested points")
plotter.show()

model.set_precursor_amounts_for_suggested()

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
#2d charged updated
'''
model=Model()
model.setup_charged({'Li':1,'Al':3,'O':-2,'F':-1})
precursors=["Li2O","LiAlO2","LiAlO1F2"]
precursors=[Composition(x) for x in precursors]
model.find_corners_edges(model.normal_vectors,add_corners_to_knowns=True)
plotter=Square(model)
model.add_precursors(precursors)
plotter.plot_mesh(
    model.precursors_constrained,poly=True,name='Precursor accessible')
#plotter.plot_points(model.omega_cut,s=0.5,c='green')
#plotter.plot_points(model.omega,s=0.5,name='all')


knowns=[Composition("LiAlO1F2"),Composition("Li4AlO3F3")]
model.add_knowns(knowns,make_plotting_df=True)
plotter.plot_plotting_df(
    plotting_columns=['Label'],c='red',name="Known phases")
model.find_spreadout_points(
    25,15,T=50,make_plotting_df=True,use_cut_omega=True)
plotter.plot_plotting_df(
    plotting_columns=['Composition'],c='green',name="Suggested points")
plotter.show()
'''
