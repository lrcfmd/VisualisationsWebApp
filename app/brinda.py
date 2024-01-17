from model import *
from visualise_square_web import *
from visualise_cube_web import *
from pymatgen.core import Composition
import plotly.express as px

#Fulle dimension
'''
model=Model()
model.setup_charged({'Cu':1,'Sn':2,'S':-2,'Cl':-1,'I':-1})
precursors=["Cu2S","CuCl","CuI","SnI2","SnS","SnCl2","SnI2"]
precursors=[Composition(x) for x in precursors]
model.find_corners_edges(model.normal_vectors,add_corners_to_knowns=True)
plotter=Cube(model)
model.add_precursors(precursors)
#plotter.plot_mesh(
#    model.precursors_constrained,poly=True,name='Precursor accessible')
#plotter.plot_points(model.omega_cut,s=0.5,c='green')
#plotter.plot_points(model.omega,s=0.5,name='all')


knowns=['Cu4S4Sn','Cu3S3.6Sn','Cu3S4Sn','Cu2.7S4Sn','Cu2S3Sn','Cu2S7Sn3',
        'Cu4S16Sn7','Cu2S9Sn4','CuS8Sn3.75','Cu5S7Sn2','Cu4S6Sn','CuI3Sn',
        'Cu3IS','Cu3IS','I6SSn4','I2SSn2','I2SSn2','I4S16Sn','Cl3ISn2',
        'ClISn']
knowns=[Composition(x) for x in knowns]
model.add_knowns(knowns,make_plotting_df=True)
plotter.plot_plotting_df(
    plotting_columns=['Label'],c='red',name="Known phases")

ratioa=[1,0,0,0,0]
ratiob=[0,1,0,0,0]
min_amount=0.1
model.add_constraints(ratioa,min_amount)
model.add_constraints(ratiob,min_amount)

model.find_spreadout_points(
    100,num_points=6,T=90,make_plotting_df=True,use_cut_omega=False,
    show_chain=False)
plotter.plot_plotting_df(
    plotting_columns=['Composition'],c='green',name="Suggested points")

plotter.plot_points(model.omega,name='points',s=1,composition=True)

plotter.show(save='../../brinda/initial_scan')

model.set_precursor_amounts_for_suggested()
model.get_precursor_df_as_csv('../../brinda/initial_scan.csv')
'''
#Iodine only
model=Model()
model.setup_charged({'Cu':1,'Sn':2,'S':-2,'I':-1})
precursors=["Cu2S","CuI","SnS","SnI2"]
precursors=[Composition(x) for x in precursors]
model.find_corners_edges(model.normal_vectors,add_corners_to_knowns=True)
plotter=Square(model)
model.add_precursors(precursors)
#plotter.plot_mesh(
#    model.precursors_constrained,poly=True,name='Precursor accessible')
#plotter.plot_points(model.omega_cut,s=0.5,c='green')
#plotter.plot_points(model.omega,s=0.5,name='all')


knowns=['Cu4S4Sn','Cu3S3.6Sn','Cu3S4Sn','Cu2.7S4Sn','Cu2S3Sn','Cu2S7Sn3',
        'Cu4S16Sn7','Cu2S9Sn4','CuS8Sn3.75','Cu5S7Sn2','Cu4S6Sn','CuI3Sn',
        'Cu3IS','Cu3IS','I6SSn4','I2SSn2','I2SSn2','I4S16Sn','Cu2SnI4']
knowns=[Composition(x) for x in knowns]
model.add_knowns(knowns,make_plotting_df=True)
plotter.plot_plotting_df(
    plotting_columns=['Label'],c='red',name="Known phases")

ratioa=[1,0,0,0]
ratiob=[0,1,0,0]
min_amount=0.05
model.add_constraints(ratioa,min_amount)
model.add_constraints(ratiob,min_amount)

model.find_spreadout_points(
    1000,num_points=5,T=90,make_plotting_df=False,use_cut_omega=False,
    show_chain=True)


#plotter.plot_points(model.omega,name='points',s=1,composition=True,prec=3)


model.set_precursor_amounts_for_suggested(lcm=6)

model.make_df_from_points_standard(model.suggested_points_standard)
plotter.plot_plotting_df(
    plotting_columns=['Composition'],c='green',name="Suggested points")
plotter.show(save='../../brinda/initial_scan')

model.get_precursor_df_as_csv('../../brinda/initial_scan.csv',comp_prec=3)
