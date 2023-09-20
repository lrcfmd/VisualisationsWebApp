import numpy as np
import plotly.express as px
import plotly.graph_objects as go
import pandas as pd
import math

class Square():
    #plotting fun ction for when constrained space is two dimensional
    
    def __init__(self,model):
        '''
        takes Model class and creates basic figure (canvas corners and edges)
        corners,edges and corner compositions of model must have already been
        found and stored
        model: Model class
        '''
        self.model=model
        self.cube_size=model.cube_size
        self.phase_field=model.phase_field
        self.create_fig_corners(
            model.corners,model.edges,model.corner_compositions)

    def set_corner_text_position(self,corners):
        '''
        Tries to find a good place to put corner text
        corners: list of corners in standard basis
        '''
        text_position=[]
        for x in corners:
            x_c=self.model.convert_to_constrained_basis(x)
            x_x=None
            if abs(x_c[0]-self.xmin)<abs(x_c[0]-self.xmax):
                x_x='left'
            else:
                x_x='right'
            x_y=None
            if abs(x_c[1]-self.ymin)<abs(x_c[1]-self.ymax):
                x_y='bottom'
            else:
                x_y='top'
            if x_c[0]==self.xmin or x_c[0]==self.xmax:
                if x_c[1]!=self.ymin and x_c[1]!=self.ymax:
                    x_y='middle'
            if x_c[1]==self.ymin or x_c[1]==self.ymax:
                if x_c[0]!=self.xmin and x_c[0]!=self.xmax:
                    x_x='center'
            text_position.append(x_y+" "+x_x)
        self.corner_text_position=text_position

    def create_fig_corners(self,corners,edges,labels,s=5):
        '''
        creates figs with corners and edges as provided
        corners: list of corners in standard basis
        edges: connectivity of corners
        labels: labels for corners
        s: size of corners
        '''
        corners_c=self.model.convert_to_constrained_basis(corners)
        columns=['x','y']
        df=pd.DataFrame(data=corners_c,columns=columns)
        self.xmin=df['x'].min()
        self.xmax=df['x'].max()
        self.ymin=df['y'].min()
        self.ymax=df['y'].max()
        self.set_corner_text_position(corners)
        df['Label']=labels
        hover_data={}
        for i in df.columns:
            hover_data[i]=False
        self.fig=px.scatter(df,x='x',y='y',text='Label',hover_data=hover_data)
        self.fig.update_traces(marker={'size':s},
                              textposition=self.corner_text_position)
        for edge in edges:
            x=[corners_c[edge[0]][0],corners_c[edge[1]][0]]
            y=[corners_c[edge[0]][1],corners_c[edge[1]][1]]
            fig1 = go.Figure(data=go.Scatter(
                x=x,y=y,mode='lines',line_color='black',hoverinfo='skip',
                showlegend=False))
            self.fig=go.Figure(data=self.fig.data+fig1.data)

    def show(self,title=None,show=True,save=None,legend_title=None,s=20, return_fig=False):
        '''
        shows and optiuonally saves figure
        '''
        self.fig.update_xaxes(showticklabels=False,autorange=False)
        self.fig.update_yaxes(showticklabels=False,autorange=True)
        self.fig.update_layout(
            font_size=s,
            hoverlabel_font_size=20,
            legend_title=legend_title,
            height=1000,
            coloraxis_colorbar=dict(yanchor="top", y=1, x=0,ticks="outside"),
            xaxis_range=[self.xmin-self.cube_size/3,self.xmax+self.cube_size/3]
        )
        self.fig.update_yaxes(
            scaleanchor = "x",
            scaleratio = 1,
        )
        if return_fig:
            return self.fig
        if title is not None:
            self.fig.update_layout(title=dict(text=title))
        if show:
            self.fig.show()
        if save is not None:
            self.fig.write_html(save+".html")
        return

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
        if not {'x','y'}.issubset(df.columns):
            raise Exception('Dataframe must have x,y columns')
        if 'z' in df.columns:
            raise Exception('Dataframe contains z but this is a 2d plotter?')

        hover_data={'x':False,
                    'y':False,
                   }
        for i in df.columns:
            if i in columns:
                hover_data[i]=True
            else:
                hover_data[i]=False
        fig1 = px.scatter(
            df,x='x',y='y',hover_data=hover_data)
        fig1.update_traces(marker_size=s)
        fig1.update_traces(marker_color=c)
        self.fig=go.Figure(data=self.fig.data+fig1.data)
