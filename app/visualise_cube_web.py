import matplotlib.pyplot as plt
import numpy as np
from algorithm import *
from matplotlib import cm
from scipy.spatial import ConvexHull, convex_hull_plot_2d
import plotly.express as px
import plotly.graph_objects as go


class Cube:

    def __init__(self,model):
        self.model=model
        self.cube_size=model.cube_size
        self.phase_field=model.phase_field
        self.create_fig_corners(
            model.corners,model.edges,model.corner_compositions)

    def create_fig_corners(self,corners,edges,labels,s=10):
        '''
        creates figs with corners and edges as provided
        corners: list of corners in standard basis
        edges: connectivity of corners
        labels: labels for corners
        s: size of corners
        '''
        corners=self.model.convert_to_constrained_basis(corners)
        corners_df=pd.DataFrame(columns=['x','y','z'])
        corners_df['x']=corners[:,0]
        corners_df['y']=corners[:,1]
        corners_df['z']=corners[:,2]
        corners_df['Label']=labels
        self.xmin=corners_df['x'].min()
        self.xmax=corners_df['x'].max()
        self.ymin=corners_df['y'].min()
        self.ymax=corners_df['y'].max()
        self.zmin=corners_df['z'].min()
        self.zmax=corners_df['z'].max()
        self.fig=px.scatter_3d(
            corners_df,x='x',y='y',z='z',text='Label',
            hover_data={'x':False,'y':False,'z':False,'Label':False})
        self.fig.update_traces(marker={'size': s})
        for edge in edges:
            x=[corners[edge[0]][0],corners[edge[1]][0]]
            y=[corners[edge[0]][1],corners[edge[1]][1]]
            z=[corners[edge[0]][2],corners[edge[1]][2]]
            fig1 = go.Figure(data=go.Scatter3d(
                x=x, y=y,z=z, mode='lines',line_color='black',
                hoverinfo='skip',showlegend=False))
            self.fig=go.Figure(data=self.fig.data+fig1.data)

    def show(self,title=None,show=True,save=None,csv=None, return_fig=False):
        self.fig.update_scenes(xaxis_visible=False)
        self.fig.update_scenes(yaxis_visible=False)
        self.fig.update_scenes(zaxis_visible=False)
        self.fig.update_layout(
            font_size=15,
            hoverlabel_font_size=20,
            coloraxis_colorbar=dict(yanchor="top", y=1, x=0,ticks="outside"),
        )
        if title is not None:
            self.fig.update_layout(
                title=dict(text=title),
                font_size=30)
        if return_fig:
            return self.fig
        if show:
            self.fig.show()
        if save is not None:
            self.fig.write_html(save+".html")
        if csv is not None:
            df=self.model.points_df
            df=df.drop('Label',axis=1)
            df.to_csv(csv)

    def plot_plotting_df(self,columns=[],s=5,c='black'):
        '''
        plots model.plotting df
        columns: columns of df to display on hover
        s: size of popints
        c: colour of points
        '''
        if not hasattr(self.model,'plotting_df'):
            raise Exception('Cant plot it if it doesnt exist')
        df=self.model.plotting_df
        if not {'x','y','z'}.issubset(df.columns):
            raise Exception('Dataframe must have x,y,z columns')

        hover_data={'x':False,
                    'y':False,
                    'z':False,}
        for i in df.columns:
            if i in columns:
                hover_data[i]=True
            else:
                hover_data[i]=False
        fig1 = px.scatter_3d(
            df,x='x',y='y',z='z',hover_data=hover_data)
        fig1.update_traces(marker_size=s)
        fig1.update_traces(marker_color=c)
        self.fig=go.Figure(data=self.fig.data+fig1.data)
