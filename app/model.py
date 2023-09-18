import numpy as np
import math
import cdd as pcdd
from scipy.optimize import linprog
import pandas as pd

class Model():
    #overall class that contains most functions for processing information

    def setup(self,charge_vector=None,cube_size=100,basis=None):
        '''
        charge vector: list(n), expected formal charges for each element.
        cube size: integer: controlls resolution of discretisation.
        converts params to right format and calls grid creator
        '''
        normalb=np.ones(len(self.phase_field))
        if charge_vector is not None:
            normala=np.array(charge_vector)
            normal_vectors=np.vstack((normala,normalb))
        else:
            normal_vectors=np.array([normalb])
        contained_point=self.get_contained_point(normal_vectors)
        contained_point=cube_size*np.array(contained_point)
        self.create_omega_constrained(
            normal_vectors,cube_size,contained_point,basis=basis)

    def setup_charged(
            self,species_dict,cube_size=30,basis=None):
        '''
        species_dict: key of element string and value of formal charge
        cube_size: resolution of grid
        basis: optional basis to use (to avoid random rotations)
        '''
        charge_vector=[]
        phase_field=[]
        for key,value in species_dict.items():
            charge_vector.append(value)
            phase_field.append(key)
        self.phase_field=phase_field
        self.setup(charge_vector,cube_size,basis)

    def setup_uncharged(
            self,species_list,cube_size=30,basis=None):
        '''
        species_list: list of element string
        cube_size: resolution of grid
        basis: optional basis to use (to avoid random rotations)
        '''
        self.phase_field=species_list
        self.setup(cube_size=cube_size,basis=basis)

    def get_contained_point(self,A):
        '''
        function to get the translation vector required for affine
        transformation to visualisation space
        A: np matrix with columns as vectors orthogonal to basis of
        visualisation space
        '''
        c=np.zeros(A.shape[1])
        c[0]=1
        b=[0]*A.shape[0]
        b[-1]=1
        return linprog(c,A_eq=A,b_eq=b).x

    def cartesian(self,arrays, out=None):
        #function for building grid (cartesian product)
        arrays = [np.asarray(x) for x in arrays]
        dtype = arrays[0].dtype

        n = np.prod([x.size for x in arrays])
        if out is None:
            out = np.zeros([n, len(arrays)], dtype=dtype)

        #m = n / arrays[0].size
        m = int(n / arrays[0].size) 
        out[:,0] = np.repeat(arrays[0], m)
        if arrays[1:]:
            self.cartesian(arrays[1:], out=out[0:m, 1:])
            for j in range(1, arrays[0].size):
            #for j in xrange(1, arrays[0].size):
                out[j*m:(j+1)*m, 1:] = out[0:m, 1:]
        return out

    def find_orthonormal(self,A):
        #finds a vector orthonormal to the column space of A
        rand_vec=np.random.rand(A.shape[0],1)
        A = np.hstack((A,rand_vec))
        b = np.zeros(A.shape[1])
        b[-1] = 1
        x = np.linalg.lstsq(A.T,b,rcond=None)[0]
        return x/np.linalg.norm(x)

    def create_omega_constrained(
            self,normal_vectors,cube_size,contained_point,basis=None):
        '''
        takes setup params and contructs finite grid spanning phase field.
        uses linear equations to reduce
        dimensionality of grid to num elements-2.
        '''

        dim=normal_vectors.shape[1]
        plane_dim = dim-len(normal_vectors)

        #basis which are equivalent up to rotation are 
        if basis is None:
            x=np.empty((plane_dim,dim))
            A = normal_vectors.T
            for i in range(plane_dim):
                x[i] = self.find_orthonormal(A)
                A = np.hstack((A,np.array([x[i]]).T))
            self.basis=x
            print('Basis used:')
            print("np."+repr(self.basis))
        else:
            self.basis=basis

        self.contained_point=contained_point

        #generate big grid of points (centred on contained p and bounded
        #by cube of length the diagonal of n-dim cube
        max_length = np.sqrt(dim*cube_size**2)
        max_co = math.floor(max_length)
        ii = np.array(range(-max_co,max_co+1,1),dtype='float64')
        if plane_dim == 3:
            omega=self.cartesian((ii,ii,ii))
        elif plane_dim ==2:
            omega=self.cartesian((ii,ii))
        elif plane_dim==1:
            omega=self.cartesian((ii))
        elif plane_dim==4:
            omega=self.cartesian((ii,ii,ii,ii))
        elif plane_dim==5:
            omega=self.cartesian((ii,ii,ii,ii,ii))
        elif plane_dim==6:
            omega=self.cartesian((ii,ii,ii,ii,ii,ii))
        else:
            raise Exception(
                'Error, plane dim: '+str(plane_dim)+', is not implemented')

        #get constraint equations specifying (n-1) dim hyperplanes
        #coresponding to 0 of each element
        line_coefficients=self.get_constraint_lines()
        #delete points the dont satisfy equations
        for line in line_coefficients:
            omega=np.delete(
                omega,np.where(
                    np.einsum(
                        '...i,i->...',omega,line[:-1])<-1*line[-1]),
                axis=0)

        self.omega=omega #the grid of points
        self.constrained_dim=plane_dim #dimension of constrained space
        self.normal_vectors=normal_vectors #normal vectors
        self.cube_size=cube_size #resolution of discretisation

    def convert_to_standard_basis(self,points,norm=None):
        '''
        converts list of points in constrained basis to standard
        basis
        points: list of points
        norm: optional L1 norm of points returned (defaults to cube size)
        '''
        if norm is None:
            norm=self.cube_size
        A=self.basis
        p_standard=self.contained_point+np.einsum('ji,...j->...i',A,points)
        if p_standard.ndim==1:
            p_standard=norm*p_standard/np.sum(p_standard)
        else:
            norma=np.sum(p_standard,axis=1)
            p_standard=(p_standard.T/norma).T
            p_standard*=norm
        return p_standard

    def convert_to_constrained_basis(self,points):
        '''
        converts list of points in standard representation to constrained
        representation.
        points: list of points
        '''
        point=np.array(points)
        if point.ndim==1:
            point=np.array(point)*self.cube_size/sum(point)
            point_s=point-self.contained_point
            point=np.einsum('ij,j',self.basis,point_s)
            return point
        else:
            norma=np.sum(point,axis=1)/self.cube_size
            point=(point.T/norma).T
            point=point-self.contained_point
            point=np.einsum('ij,...j->...i',self.basis,point,dtype=float)
            return point

    def find_spreadout_points(self,n):
        '''
        finds n spread out points in space and stores them in plot ready df
        n: number of points
        '''
        if self.constrained_dim!=2 and self.constrained_dim!=3:
            raise Exception('only implemented for dim=2 or 3')
        points=self.omega[
            np.random.choice(self.omega.shape[0],n,replace=False)]
        x=points[:,0]
        y=points[:,1]
        df=pd.DataFrame()
        df['x']=x
        df['y']=y
        if self.constrained_dim==3:
            df['z']=points[:,2]
        points_stan=self.convert_to_standard_basis(points)
        for n,el in enumerate(self.phase_field):
            df[el]=points_stan[:,n]
        df['Composition']=self.get_composition_strings(points_stan)
        self.plotting_df=df


    def get_composition_strings(self,points_stan,l1_norm=1,out='html'):
        '''
        function to get a string representation for a composition
        should probably replace this with pymatgen
        points_stan: list of points in standard basis
        l1_norm: L1 norm
        out: output format(for the subscripts basically)
        '''
        points_stan=np.array(points_stan)
        if points_stan.ndim==1:
            points_stan=np.array([points_stan])
        comps=[]
        for point in points_stan:
            point=l1_norm*point/sum(point)
            comp=''
            if out=='html':
                for x,el in zip(point,self.phase_field):
                    x=round(x,2)
                    if x!=0:
                        comp+=el+'<sub>'
                        comp+=str(x)+'</sub>'
                comps.append(comp)
            else:
                raise Exception('Unknwon output format')
        if len(comps)==1:
            return comps[0]
        else:
            return comps

    def find_corners_edges(self,A):
        '''
        gets corners and edges of polyhedra
        A:2d np matrix with columns as the constraining linear equations
        returns in standard rep
        '''
        A = A.T
        cp=self.get_contained_point(A.T)
        dim=A.shape[0]
        plane_dim=dim-A.shape[1]
        x=np.empty((plane_dim,dim))
        for i in range(plane_dim):
            x[i] = self.find_orthonormal(A)
            A = np.hstack((A,np.array([x[i]]).T))
        charge_basis=x
        cons=[]
        for i in range(dim):
            con=[0]*dim
            con[i]=1
            cons.append(np.array(con))
        cons_c=[]
        if plane_dim==2:
            for i in cons:
                a=np.dot(i,charge_basis[0])
                b=np.dot(i,charge_basis[1])
                c=np.dot(i,cp)
                cons_c.append([c,a,b])
        elif plane_dim==3:
            for i in cons:
                a=np.dot(i,charge_basis[0])
                b=np.dot(i,charge_basis[1])
                c=np.dot(i,charge_basis[2])
                d=np.dot(i,cp)
                cons_c.append([d,a,b,c])
        else:
            raise Exception('Not implemented for polytope above 3d')
        #create matric for polyhedra construction
        mat=pcdd.Matrix(cons_c,linear=False)
        mat.rep_type = pcdd.RepType.INEQUALITY
        #create polyhedron
        poly = pcdd.Polyhedron(mat)
        #get corners
        ext = poly.get_generators()
        corners=[]
        for i in ext:
            corners.append(i[1:])
        corners=np.array(corners)
        #get edges
        edges=[]
        adjacencies = list(poly.get_adjacency())
        for n,i in enumerate(adjacencies):
            for j in list(i):
                if j>n:
                    edges.append([n,j])
        corners=cp+np.einsum('ji,...j->...i',charge_basis,corners)
        corner_compositions=self.get_composition_strings(corners)
        self.corners=corners
        self.edges=edges
        self.corner_compositions=corner_compositions

    def get_constraint_lines(self):
        #function that uses the assigned basis vectors to find the planes
        #corresponding to x_i>0 in constrained space
        line_coefficients=[]
        for i in range(self.basis.shape[1]):
            coefficients=[]
            for j in range(self.basis.shape[0]):
                coefficients.append(self.basis[j][i])
            coefficients.append(1*self.contained_point[i])
            line_coefficients.append(np.array(coefficients))
        return line_coefficients

