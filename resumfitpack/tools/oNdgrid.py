#!/usr/bin/env python
import numpy as np
import pylab as py
from scipy.spatial         import ConvexHull
from scipy.spatial         import Delaunay
from sklearn.preprocessing import MinMaxScaler
from itertools import combinations

class ONDGRID():

    def __init__(self,points,alpha=1,rad=None,umin=0,umax=1,gridgoal=1000):

        self.points  = points
        self.alpha   = alpha
        self.length  = len(points)
        self.dim     = points.shape[1]
        self.rad     = rad

        self.umin   = umin
        self.umax   = umax

        self.scaler  = MinMaxScaler(feature_range=(0.1, 0.9))
        self.upoints = self.scaler.fit_transform(points)
        self.upmin=[]
        self.upmax=[]
        for i in range(self.dim):
            self.upmin.append(min(self.upoints.T[i]))
            self.upmax.append(max(self.upoints.T[i]))
        print self.upmin,self.upmax
        #self.hull1   = ConvexHull(self.upoints)
        self.hull1   = ConvexHull(self.points)
        print 'hull is done'
        self.gridgoal=gridgoal
        #self.gen_hull2()
        self.gen_grid()

    def get_min_rad(self):
        r2=1e100
        comb=list(combinations(range(self.dim),2))
        for simplex in self.hull1.simplices:
            coord=[self.upoints[simplex,i] for i in range(self.dim)]
            #coord=[self.points[simplex,i] for i in range(self.dim)]
            #--for each simplex, there are self.dim choose 2 radii to calculate
            for idx in comb:
                i0,i1=idx[0],idx[1]
                _r2=0.0
                for j in range(len(coord)):
                    _r2+=(coord[j][i0]-coord[j][i1])**2
                if r2>_r2:  r2=_r2
        rad=np.sqrt(r2)*self.alpha
        if self.rad!=None: rad=self.rad
        return rad

    def gen_hull2(self):

        print 'generating hull2'
        r=self.get_min_rad()
        deltacoord=[np.zeros(2*self.dim) for _ in range(self.dim)]
        for i in range(self.dim):
            deltacoord[i][i]=r
            deltacoord[i][i+self.dim]=-r

        epoints=[]
        for simplex in self.hull1.simplices:
            coord=[self.upoints[simplex,i] for i in range(self.dim)]
            #coord=[self.points[simplex,i] for i in range(self.dim)]
            for i in range(self.dim):
                ex=[]
                for ii in range(2*self.dim):
                    coordchange=[]
                    for j in range(self.dim):
                        coordchange.append(coord[j][i]+deltacoord[j][ii])
                    ex.append(coordchange)
                epoints.extend(ex)
        epoints=np.array(epoints)

        print 'epoints gathered',epoints.shape
        print np.unique(epoints).shape
        self.hull2   = ConvexHull(epoints)
        print 'hull2 drawn'
        self.epoints = epoints

    def gen_grid(self,grid=None):
        print 'setting up grid'

        #dhull = Delaunay(self.hull2.points)
        dhull = Delaunay(self.hull1.points)
        ogrid=[]
        goal=self.gridgoal
        if grid==None:
            rand=np.random.rand(goal,self.dim)
            for i in range(self.dim):
                for j in range(len(rand)):
                    #rand[j][i]=self.pstd[i]*rand[j][i]+self.pmean[i]
                    rand[j][i]=self.upmin[i]+(self.upmax[i]-self.upmin[i])*rand[j][i]
                    #if self.upmax[i]-self.upmin[i]<1e-3: rand[j][i]=self.upmin[i]
            #grid=0.1+0.8*np.random.rand(goal,self.dim)
            grid=rand
            print grid
            length=self.length
            cnt=0
            while length<length+1:#goal:
                if cnt==100: print 'cnt at 100'
                if cnt==1000: print 'cnt at 1000'
                if cnt==10000: print 'cnt at 10000'
                #print 'doing the grid, length is %d'%length
                if length==16: print 'length gained 1'
                if length==100: print 'length at 100'
                if length==200: print 'length at 200'
                for point in grid:
                    if dhull.find_simplex(point)>=0:
                        ogrid.append([point[i] for i in range(len(point))])
                        length+=1
                #grid=np.random.rand(goal,self.dim)
                rand=np.random.rand(goal,self.dim)
                for i in range(self.dim):
                    for j in range(len(rand)):
                        #rand[j][i]=self.pstd[i]*rand[j][i]+self.pmean[i]
                        rand[j][i]=self.upmin[i]+(self.upmax[i]-self.upmin[i])*rand[j][i]
                        #if self.pmax[i]-self.pmin[i]<1e-3: rand[j][i]=self.pmin[i]*1.00000001
                #grid=0.1+0.8*np.random.rand(goal,self.dim)
                grid=rand
                cnt+=1
            print 'got out of the loop; length gained 1'

        #ogrid=[]
        for point in grid:
            if dhull.find_simplex(point)>=0 :
                #ogrid.append([point[0],point[1],point[2]])
                ogrid.append([point[i] for i in range(len(point))])

        self.grid  = grid
        self.ogrid = ogrid

    def get_grid(self):
        return self.scaler.inverse_transform(self.ogrid)

    def plot(self):

        #--Needs improvement

        ax=py.subplot(111)

        #--plot original scaled points
        X,Y=np.transpose(self.upoints)
        ax.plot(X,Y,'ro',label='original points')

        #--plot convex hull1
        for simplex in self.hull1.simplices:
            ax.plot(self.upoints[simplex, 0],self.upoints[simplex, 1], 'b-')

        #--plot convex hull2
        for simplex in self.hull2.simplices:
            ax.plot(self.epoints[simplex, 0],self.epoints[simplex, 1], 'b-')

        #--plot optimal grid
        X,Y=np.transpose(self.ogrid)
        ax.plot(X,Y,'k.',label='optimal grid points')

        ax.legend()
        #ax.set_ylim(0,1)
        #ax.set_xlim(0,1)
        py.savefig('o3grid.pdf')
        py.close()
 
def example00():


    points=np.random.randn(10,3)
    og=O2DGRID(points,1)
    #og.plot()

if __name__=='__main__':

    example00()




