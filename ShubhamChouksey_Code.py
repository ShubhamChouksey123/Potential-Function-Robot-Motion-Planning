

#importing the necessary libraries
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg 
import math # for using floor function
from queue import Queue  #using queue in brushfire Algorithm


# Converting a RGB image into grayscale in Python
def rgb2gray(rgb):
    return np.dot(rgb[...,:3], [0.2989, 0.5870, 0.1140])

img = mpimg.imread('cs.png')     
gray = rgb2gray(img)    
plt.imshow(gray, cmap=plt.get_cmap('gray'), vmin=0, vmax=1)
I = gray


N = 25;         # Selecting 25 as the number of grids for BrushFire Algorithm
Inorm = (I);    # Inorm is 100 X 100 Matrix
world = np.zeros((N+1, N+1)) # world is 25 X 25 Matrix Grid with each entry as 0  #1 based indexing 



for i in range(0, 100, int(100/N) ):
    for j in range(0, 100, int(100/N)):
        minval = 1000000
        for k in range(i, i + int(100/N) , 1):
            for l in range(j, j + int(100/N) , 1):
                minval = min(minval, Inorm[k][l])
        if minval == 0:
            # Making Polygon as a box of red color with four vertices cordinate which are appended in the obs1 array
            obs1 = [(j,i), (j+(100/N)-1,i), (j+(100/N)-1,i+(100/N)-1), (j, i+(100/N)-1)]
            obs1.append(obs1[0])   #repeating the first point to create a 'closed loop'
            xs1, ys1 = zip(*obs1)  #creating lists of x and y values
            plt.fill(xs1,ys1, 'r', alpha=0.7) 
            world[1+math.floor(i/(100/N))][1+ math.floor(j/(100/N))] =1     
            
        else :
            obs1 = [(j,i), (j+(100/N)-1,i), (j+(100/N)-1,i+(100/N)-1), (j, i+(100/N)-1)]
            obs1.append(obs1[0])   #repeating the first point to create a 'closed loop'
            xs1, ys1 = zip(*obs1)  #creating lists of x and y values
            plt.fill(xs1,ys1, 'g', alpha=0.7) 
            world[1+math.floor(i/(100/N))][1+ math.floor(j/(100/N))] = 0     




#### BrushFire Algorithm on the World Matrix ##############


def CheckZeroInMat(world):
    for i in range(1, 26, 1):
        for j in range(1, 26, 1):
            if world[i][j] == 0:
                return 1
    return 0

def LiesInRange(x, y):
    if x >= 1 and x <= N and y >= 1 and y <= N: 
        return 1
    return 0

cur = 1
count = 1000
# neighbour = [[1, 0],[0 ,1],[-1, 0], [0, -1]]

neighbour = [[1,0], [1, 1], [0, 1], [-1, 1], [-1, 0], [-1, -1], [0, -1],[1, -1]]



que = Queue()

for i in range(1, 26, 1): 
    for j in range(1, 26, 1): 
        if world[i][j] == 1: 
            que.put([i,j, 2]); 


while que.empty() == 0 : 
    val = que.get()
    x = val[0]; y = val[1]; dist = val[2];
    if LiesInRange(x + 1, y) and world[x + 1][y] == 0: 
        world[x + 1][y] = dist
        que.put([x+1, y, dist+1])
    if LiesInRange(x , y+1) and world[x][y+1] == 0: 
        world[x][y+1] = dist;
        que.put([x, y+1, dist+1]);
    if LiesInRange(x-1 , y) and world[x-1][y] == 0: 
        world[x-1][y] = dist;
        que.put([x-1, y, dist+1]);
    if LiesInRange(x , y-1) and world[x][y-1] == 0: 
        world[x][y-1] = dist;
        que.put([x, y-1, dist+1]);
    if LiesInRange(x + 1, y+1) and world[x + 1][y+1] == 0: 
        world[x + 1][y+1] = dist
        que.put([x+1, y+1, dist+1])
    if LiesInRange(x-1 , y+1) and world[x-1][y+1] == 0: 
        world[x-1][y+1] = dist;
        que.put([x-1, y+1, dist+1]);
    if LiesInRange(x-1 , y-1) and world[x-1][y-1] == 0: 
        world[x-1][y-1] = dist;
        que.put([x-1, y-1, dist+1]);
    if LiesInRange(x+1 , y-1) and world[x+1][y-1] == 0: 
        world[x+1][y-1] = dist;
        que.put([x+1, y-1, dist+1]);  

for i in range(1, 26, 1): 
    for j in range(1, 26, 1):
        x = (i)*int(100/N) -2
        y = (j)*int(100/N) -2
        plt.annotate('%d'%(world[i][j]), xy=(y-1.8 , x+0.5))
        
plt.show()



###########################################################################################
# Repulive Function 


# print("\n\n New World = ")
for i in range(0, 25, 1):
    for j in range(0, 25, 1):
        world[i][j] = world[i+1][j+1]


a = np.delete(world, 25, 0)
a = np.delete(a, 25, 1)
world = a

eta = 5;                     
[dworldx, dworldy] = np.gradient(world);

#swapping dworldx and dworldy as they are coming swaped
temp = dworldx
dworldx = dworldy
dworldy = temp

  
# Creating MeshGrid Method 3
X , Y = np.meshgrid(np.linspace(0, 99, 100), np.linspace(0, 99, 100))

Urep = np.zeros((100, 100))  # Urep an array of size 100 x 100
Qstar = 4;
gradUrep = np.zeros((100, 100, 2))


for i in range(0,100, 1):
    for j in range(0,100, 1):
        Dq = world[math.floor(i/(100/N))][math.floor(j/(100/N))]
        gradDq = [dworldx[math.floor(i/(100/N))][math.floor(j/(100/N))], dworldy[math.floor(i/(100/N))][math.floor(j/(100/N))] ]
        if Dq <= Qstar:
            Urep[i][j]  = 0.5*eta*math.pow(((1/Dq) - (1/Qstar)), 2)
            gradUrep[i][j][0] = eta*((1/Qstar) - (1/Dq))*(1/pow(Dq, 2))*gradDq[0]; 
            gradUrep[i][j][1] = eta*((1/Qstar) - (1/Dq))*(1/pow(Dq, 2))*gradDq[1];
        else :
            Urep[i][j] = 0; 
            gradUrep[i][j][0] = 0
            gradUrep[i][j][1] = 0;

fig = plt.figure() # creates space for a figure to be drawn 
# Uses a 3d prjection as model is supposed to be 3D
axes = fig.gca(projection ='3d')
# Plots the three dimensional data consisting of x, y and z 
axes.plot_surface(X, Y, Urep)   # plotting repulsive potential function
# show command is used to visualize data plot   
plt.show()            


gradUrep[:, :, 0], gradUrep[: ,:, 1] = np.gradient(Urep)

dx , dy = np.gradient(Urep)
img = mpimg.imread('cs.png')     
gray = rgb2gray(img)    
plt.imshow(gray, cmap=plt.get_cmap('gray'), vmin=0, vmax=1)

# for coloring the quiver plot 
# n = -2
# color = np.sqrt(((dx-n)/2)*2 + ((dy-n)/2)*2)

plt.quiver(X, Y, gradUrep[:, :, 1], -gradUrep[:, :, 0], 0, alpha = 0.5) 
  
plt.show()

########################################################################################
# quadratic cum conic potential
# Attractive Potential

zeta = 0.1;
dstar = 6;
q_goal = [20, 60];
X , Y = np.meshgrid(np.linspace(0, 99, 100), np.linspace(0, 99, 100))
Uatt = np.zeros((100, 100))  # Urep an array of size 100 x 100
gradUatt = np.zeros((100, 100, 2))


for i in range(0,100, 1):
    for j in range(0,100, 1):
        q = [i, j]
        if np.linalg.norm(np.array(q_goal) - np.array(q)) <= dstar:
            Uatt[i][j]  = 0.5*zeta*math.pow(np.linalg.norm(np.array(q_goal) - np.array(q)), 2)
            gradUatt[i][j][0] = zeta*(q_goal[0] - q[0]); 
            gradUatt[i][j][1] = zeta*(q_goal[1] - q[1]);
        else :
            Uatt[i][j] = ( dstar*zeta*np.linalg.norm(np.array(q_goal) - np.array(q)) ) - (0.5*zeta*math.pow(dstar, 2)) 
            gradUatt[i][j][:] = dstar*zeta*(np.array(q_goal) - np.array(q))/np.linalg.norm(np.array(q_goal) - np.array(q))
        
fig = plt.figure() # creates space for a figure to be drawn 
# Uses a 3d prjection as model is supposed to be 3D
axes = fig.gca(projection ='3d')
# Plots the three dimensional data consisting of x, y and z 
axes.plot_surface(X, Y, Uatt)   # plotting repulsive potential function
# show command is used to visualize data plot   
plt.show()            



dx , dy = np.gradient(Uatt)
img = mpimg.imread('cs.png')     
gray = rgb2gray(img)    
plt.imshow(gray, cmap=plt.get_cmap('gray'), vmin=0, vmax=1)

# for coloring the quiver plot 
# n = -2
# color = np.sqrt(((dx-n)/2)*2 + ((dy-n)/2)*2)


plt.quiver(X, Y, gradUatt[:, :, 0], gradUatt[:, :, 1], 0,scale=0.7, alpha = 0.5,  units="xy") 
  
plt.show()


################################## Total Potential ##############################

U = np.zeros((100, 100))  # Urep an array of size 100 x 100
gradU = np.zeros((100, 100, 2))

U = Uatt + 7*Urep;
gradU[:, :, 0], gradU[: ,:, 1] = np.gradient(U)

fig = plt.figure() # creates space for a figure to be drawn 
# Uses a 3d prjection as model is supposed to be 3D
axes = fig.gca(projection ='3d')
# Plots the three dimensional data consisting of x, y and z 
axes.plot_surface(X, Y, U)   # plotting repulsive potential function
# show command is used to visualize data plot   
plt.show()            



dx , dy = np.gradient(U)
img = mpimg.imread('cs.png')     
gray = rgb2gray(img)    
plt.imshow(gray, cmap=plt.get_cmap('gray'), vmin=0, vmax=1)

# for coloring the quiver plot 
# n = -2
# color = np.sqrt(((dx-n)/2)*2 + ((dy-n)/2)*2)

plt.quiver(X, Y, gradU[:, :, 1], -gradU[:, :, 0], 0, alpha = 0.5) 

       

############## Gradient Descent #############################################

qs = [20, 70];

plt.plot(qs[0],qs[1], 'bo')

q = qs;
alpha = 3;
grad = [gradU[q[1]][q[0]][0], gradU[q[1]][q[0]][1]];
while np.linalg.norm(grad)  > 0.2:
    if q == q_goal :
        break;
    plt.plot(qs[0],qs[1], 'go')
    grad = [gradU[math.floor(q[1])][math.floor(q[0])][0], gradU[math.floor(q[1])][math.floor(q[0])][1]]
    
    q[0] = round(q[0] - alpha * grad[1]);
    q[1] = round(q[1] - alpha * grad[0]);


plt.show()
