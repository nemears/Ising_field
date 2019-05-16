import numpy as np
import pyglet
import imageio
from random import choice, random
import sys
import matplotlib.pyplot as plt
import csv
import Tools as tp

def noField(i,j,time):
    return 0

def cField(i,j,time):
    return 1

def negField(i,j,time):
    return -1

def sino(i,j,time):
    return np.sin(0.018*time)

def linear(i,j,time):
    return (i-8)/16

def travellingWave(i,j,time):
    A = 1
    w = 0.018
    return A*np.sin(16*i-w*time)

class Lattice(object):

    def __init__(self,xdim,ydim,beta,field,animate=False,save=False):
        self.xdim = xdim #dimensions of lattice
        self.ydim = ydim
        self.beta = beta # inverse temperature
        self.field = field # function of space and time representing the field H
        self.animate = animate # boolean whether to animate
        self.save = save # boolean whether to save
        self.lattice = np.random.choice([1,-1],(xdim,ydim))
        self.save_array = []
        self.final_array = []

        # Setting up animation
    
        if self.animate:
            window = pyglet.window.Window()
            width = 500
            height = 500
            window.set_size(width,height)
            xsize = window.width//xdim
            ysize = window.height//ydim
            xgap = window.width - (xsize*xdim)
            ygap =  window.height - (ysize*ydim)

            self.upImage = pyglet.image.load('bs.png')
            self.upImage.width = xsize
            self.upImage.height = ysize
            self.downImage = pyglet.image.load('ws.png')
            self.downImage.width = xsize
            self.downImage.height = ysize

            self.sprites = []
            self.main_batch = pyglet.graphics.Batch()

            for i in range(self.xdim):
                for j in range(self.ydim):
                    if self.lattice[i][j] == 1:
                        self.sprites.append(pyglet.sprite.Sprite(self.upImage,i*xsize+xgap//2
                                                ,j*ysize+ygap//2,batch=self.main_batch))
                    else:
                        self.sprites.append(pyglet.sprite.Sprite(self.downImage,i*xsize+xgap//2
                                                ,j*ysize+ygap//2,batch=self.main_batch))

    def getAdjacent(self,i,j): # gets all adjacent coords to (i,j)
        coords = []
        if i != 0:
            coords.append(self.lattice[i-1][j])
        else:
            coords.append(self.lattice[self.xdim-1][j])
        if i != self.xdim-1:
            coords.append(self.lattice[i+1][j])
        else:
            coords.append(self.lattice[0][j])
        if j != 0:
            coords.append(self.lattice[i][j-1])
        else:
            coords.append(self.lattice[i][self.ydim-1])
        if j != self.ydim-1:
            coords.append(self.lattice[i][j+1])
        else:
            coords.append(self.lattice[i][0])
        return coords

    def getH(self,i,j,t): # gets total energy of spin site
        spin = self.lattice[i][j]
        adjacent = self.getAdjacent(i,j)
        H = 0
        for a in adjacent:
            H+= -a*spin
        return H - self.field(i,j,t)*spin

    def tick(self,time): # checks certain number of spins to flip
        self.save_array.append(self.lattice.copy())
        for e in range(1):#(self.xdim*self.ydim)//64):
            i,j = choice(range(self.xdim)),choice(range(self.ydim))
            H = self.getH(i,j,time)
            if H > 0 or random() < np.exp(-self.beta*-2*H):
                self.lattice[i][j] = self.lattice[i][j]*-1
            if self.animate:
                if self.lattice[i][j] == 1:
                    self.sprites[i*self.ydim+j].image = self.upImage
                else:
                    self.sprites[i*self.ydim+j].image = self.downImage
        if self.animate:
            self.main_batch.draw()

    def findEq(self,tolerance,thresh,dom,time): # reduces runtime
        if time<tolerance: #tolerance run time
            return True
        elif time<4*tolerance:
            vv = np.std([np.average(s) for s in self.save_array[time-dom::10]])
            return vv >= thresh*(time/tolerance) # False if at equilibrium
        else:
            return False
    
    def run(self,nframe,tillEq=False,message=True): # runs the lattice simulation
        if self.animate:
            pyglet.clock.schedule_interval(self.tick, 1/60)
            pyglet.app.run()
        else:
            if tillEq:
                i = 0
                dom = 0
                while np.abs(np.average(self.lattice)) < 1 and self.findEq(nframe,0.001,dom,i):
                    self.tick(i)
                    i+=1
                    if dom < nframe**0.5:
                        dom+=1
                    if i%3600 == 0 and message:
                        print('{} mins'.format(i/3600))
            else:
                for i in range(nframe):
                    self.tick(i)
                    if i%3600==0 and message:
                        print('{} mins'.format(i/3600))

    def getM(self):
        M=[]
        for s in self.save_array:
            M.append(np.abs(np.average(s)))
        return M

    def getStd(self, H=False):
        S =[]
        if H:
            t = 0
            h_array = np.zeros((self.xdim,self.ydim))
            for s in self.save_array:
                for i in range(self.xdim):
                    for j in range(self.ydim):
                        h_array[i][j] = self.getH(i,j,t)
                t+=1
                S.append(np.std(h_array))
        else:
            for s in self.save_array:
                S.append(np.std(s))
        return S

    def saveFrame(self,lattice):
        lwidth,lheight = np.shape(self.lattice)
        temp = np.zeros((512,512),dtype=np.uint8)
        cpixel = 0
        for i in range(lwidth):
            for j in range(lheight):
                if lattice[i][j] == -1:
                    cpixel = 0
                else:
                    cpixel = 1
                for k in range(512//lwidth):
                    for l in range(512//lheight):
                        temp[i*(resolution[0]//lwidth)+k][j*(resolution[1]//lheight)+l] = (cpixel*255)
        self.final_array.append(temp)

    def saveMov(self,name):
        for s in self.save_array:
            self.saveFrame(s)
        imagio.mimwrite(name+'.mp4',self.final_array,fps=60)

def magvstd(): # runs sim for single lattice and gets m and stdS
    M = []
    S = []
    nFrame = 4000
    sim_num = 100
    for i in range(sim_num):
        lat = Lattice(16,16,10,sino,False,False)
        lat.run(nFrame)
        M.append(lat.getM())
        S.append(lat.getStd())
    K = [np.abs(np.average([M[m][i] for m in range(sim_num)])) for i in range(nFrame)]
    T = [np.average([S[m][i] for m in range(sim_num)]) for i in range(nFrame)]
    return K,T

def tCurie():
    sim_num = 30
    bmax = 1
    lT = []
    x = np.linspace(0,bmax,100)
    for b in x:
        print('b =', b)
        temp =[]
        for i in range(sim_num):
            lat = Lattice(16,16,b,sino)
            lat.run(5000,tillEq=True,message=True)
            temp.append(np.abs(np.average(lat.lattice)))
        lT.append(np.average(temp))
        print('average for {} is {}'.format(b,np.average(temp)))
    with open('sinLongMvsB.csv','w') as file:
        writer = csv.writer(file)
        writer.writerow(x)
        writer.writerow(lT)
    plt.plot(x,lT)
    plt.plot(x,[(1-(np.sinh(2*b))**(-4))**(0.125) for b in x])
    plt.show()
        
if __name__ == '__main__':
    
    def tWave(A,w,k,i,j,time):
        return A*np.sin(k*i-w*time)

    def tDec(A,w,k,func):
        def inner(i,j,time):
            return func(A,w,k,i,j,time)
        return inner

    funcs = []

    A = np.linspace(0,5,10)
    W = np.linspace(0,0.002,10)
    K = np.linspace(0,100,10)
    
    for i in W:
        funcs.append(tDec(5,i,32,tWave))

    means = []
    std = []
    lats = []
    acfs = []
    peaks = []
    i = 0
    
    for f in funcs:
        print(i)
        lat = Lattice(16,16,10,f)
        lat.run(100000)
        means.append(np.abs(np.mean(lat.lattice)))
        std.append(np.std(lat.lattice))
        lats.append(lat.lattice)
        acfs.append(tp.acf(lat.lattice))
        peaks.append(tp.getAcfPeaks(acfs[i]))
        i+=1

    for i in range(len(A)):
        plt.plot(peaks[i],[K[i] for p in peaks[i]],'bo')
    plt.show()
















                
