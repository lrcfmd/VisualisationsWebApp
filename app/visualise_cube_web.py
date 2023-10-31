import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
import pandas as pd
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

    def plot_plotting_df(
            self,plotting_columns=[],s=None,c=None,name=None,df=None):
        '''
        plots model.plotting df
        columns: columns of df to display on hover
        s: size of popints
        c: colour of points
        '''
        if df is None:
            if not hasattr(self.model,'plotting_df'):
                raise Exception('Cant plot it if it doesnt exist')
            df=self.model.plotting_df
        if not {'x','y','z'}.issubset(df.columns):
            raise Exception('Dataframe must have x,y,z columns')

        hover_data={'x':False,
                    'y':False,
                    'z':False,}
        for i in df.columns:
            if i in plotting_columns:
                hover_data[i]=True
            else:
                hover_data[i]=False
        fig1 = px.scatter_3d(
            df,x='x',y='y',z='z',hover_data=hover_data)
        if s is not None:
            fig1.update_traces(marker_size=s)
        if c is not None:
            fig1.update_traces(marker_color=c)
        if name is not None:
            fig1.update_traces(name=name,showlegend=True)
        self.fig=go.Figure(data=self.fig.data+fig1.data)

    def plot_mesh(self,points,name=None,poly=False):

        x=points[:,0]
        y=points[:,1]
        z=points[:,2]
        if poly:
            hull=ConvexHull(points)
            i=hull.simplices[:,0]
            j=hull.simplices[:,1]
            k=hull.simplices[:,2]
            fig = go.Figure(
                data=[go.Mesh3d(
                    x=x, y=y, z=z, i=i,j=j,k=k,color='lightpink',
                    opacity=0.50)])
        if name is not None:
            fig.update_traces(showlegend=True,name=name)
        self.fig=go.Figure(data=self.fig.data+fig.data)

    def plot_points(self,points,name=None,s=None,c=None):
        df=pd.DataFrame()
        df['x']=points[:,0]
        df['y']=points[:,1]
        df['z']=points[:,2]
        self.plot_plotting_df(df=df,name=name,s=s,c=c)

