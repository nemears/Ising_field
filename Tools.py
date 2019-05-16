import imageio
import numpy as np

def saveFrame(lattice,sMax,sMin, resolution, color = False): # Takes lattice and gets it ready for export to mp4
    t_lat = lattice.copy()
    for f in np.nditer(t_lat):
        f = f-sMin
    lwidth,lheight = np.shape(t_lat)
    temp = np.zeros(resolution,dtype=np.uint8)
    cpixel = 0
    for i in range(lwidth):
        for j in range(lheight):
            cpixel = t_lat[i][j]/(sMax-sMin)
            for k in range(resolution[0]//lwidth):
                for l in range(resolution[1]//lheight):
                    temp[i*(resolution[0]//lwidth)+k][j*(resolution[1]//lheight)+l] = (cpixel*255)
    if color:
        red = 127 + 0.5*temp
        green = temp
        blue = (255-temp)
        return red,green,blue
    else:
        return temp

def checkBounds(coords,x,y,b1,b2): # Makes sure coords are in the domain of lattice
    if x >= b1 and x <= b2 and y >= b1 and y <= b2:
        coords.append((x,y))
    elif x<b1: # Periodic boundaries only on x
        coords.append((b2+x, y))
    elif x > b2:
        coords.append((x-b2,y))
        
def getCircle(cx,cy,r,bounds=(0,0)):
    coords = []
    x = r
    y = 0
    while x>=y:
        coords.append((x,y))
        if x**2 + (y+1)**2 > r**2:
            x-=1
        y+=1
    tcoords = []
    if bounds == (0,0):
        for c in coords:
            tcoords.append((-c[0]+cx,-c[1]+cy))
            tcoords.append((-c[0]+cx,c[1]+cy))
            tcoords.append((c[0]+cx,-c[1]+cy))
            tcoords.append((c[0]+cx,c[1]+cy))
            tcoords.append((-c[1]+cx,-c[0]+cy))
            tcoords.append((-c[1]+cx,c[0]+cy))
            tcoords.append((c[1]+cx,-c[0]+cy))
            tcoords.append((c[1]+cx,c[0]+cy))
        return tcoords
    else:
        for c in coords:
            checkBounds(tcoords,-c[0]+cx,-c[1]+cy,b1,b2)
            checkBounds(tcoords,-c[0]+cx,c[1]+cy,b1,b2)
            checkBounds(tcoords,c[0]+cx,-c[1]+cy,b1,b2)
            checkBounds(tcoords,c[0]+cx,c[1]+cy,b1,b2)
            checkBounds(tcoords,-c[1]+cx,-c[0]+cy,b1,b2)
            checkBounds(tcoords,-c[1]+cx,c[0]+cy,b1,b2)
            checkBounds(tcoords,c[1]+cx,-c[0]+cy,b1,b2)
            checkBounds(tcoords,c[1]+cx,c[0]+cy,b1,b2)
        return tcoords

def acf(array,bounds=(0,0)):
    hwidth, hheight = np.shape(array)
    fftim = np.fft.fft2(array)
    conjim = np.conj(fftim)
    ans1 = fftim*conjim
    ans2 = np.fft.ifft2(ans1)/(hwidth*hheight)
    ans3 = ans2.real
    global ans4
    ans4 = np.zeros(np.shape(ans3))

    for i in range(hwidth):
        for j in range(hheight):
            ans4[i][j]=ans3[i][j]-np.min(ans3)
    
    ans4 = ans4/np.max(ans4)
    ans5 = np.fft.fftshift(ans4)
    cor = []
    for r in range(min(hwidth//2,hheight//2)):
        coords = getCircle(hwidth//2,hheight//2,r,bounds)
        tot = 0
        for c in coords:
            tot+=ans5[c[0]][c[1]]
        cor.append(tot/len(coords))
    return cor

def getAcfPeaks(acf):
    peaks = [0]
    for i in range(1,len(acf)):
        if acf[i] > acf[i-1]:
            if peaks[-1] == i-1:
                peaks.pop(-1)
            peaks.append(i)
    return peaks

def getAdjacent(i,j,periodicX = False, periodicY = False,bounds = (0,0)): # gets all adjacent coords to (i,j)
        coords = []
        if periodicX and bounds != (0,0):
            if i != bounds[0]:
                coords.append((i-1,j))
            else:
                coords.append((bounds[1],j))
            if i != bounds[1]:
                coords.append((i+1,j))
            else:
                coords.append((bounds[0],j))
            if j != bounds[0]:
                coords.append((i,j-1))
            if j != bounds[1]:
                coords.append((i,j+1))
            return coords
        elif periodicY and bounds != (0,0):
            if i != bounds[0]:
                coords.append((i-1,j))
            if i != bounds[1]:
                coords.append((i+1,j))
            if j != bounds[0]:
                coords.append((i,j-1))
            else:
                coords.append((i,bounds[1]))
            if j != bounds[1]:
                coords.append((i,j+1))
            else:
                coords.append((i,bounds[0]))
            return coords
        elif bounds == (0,0):
            coords.append(self.lattice[i-1][j])
            coords.append(self.lattice[i+1][j])
            coords.append(self.lattice[i][j-1])
            coords.append(self.lattice[i][j+1])
            return coords
        else:
            if i != bounds[0]:
                coords.append((i-1,j))
            else:
                coords.append((bounds[1],j))
            if i != bounds[1]:
                coords.append((i+1,j))
            else:
                coords.append((bounds[0],j))
            if j != 0:
                coords.append((i,j-1))
            else:
                coords.append((i,M-1))
            if j != M-1:
                coords.append((i,j+1))
            else:
                coords.append((i,0))
            return coords

def getPD(vals,spread):
    vals.sort()
    mn = min(vals)
    mx = max(vals)
    x = np.linspace(mn+((mx-mn)/(spread-1)),mx,spread-1)
    y = np.zeros(spread-1)
    i = 0
    for j in x:
        while vals[0] < j:
            y[i]+=1
            vals.pop(0)
        i+=1
    return x,y
