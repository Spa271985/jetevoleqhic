
import matplotlib.pyplot as plt 
import math
from math import sqrt, exp
from random import random
import numpy as np
from matplotlib import pyplot

#Defining probability functions and functions containing expressions for weight

epsilon=0.0001
kap=2*(1/sqrt(epsilon)-1+sqrt(1-epsilon))
kap1=2*(1/sqrt(epsilon)-1)

def phi(x):
    return kap/sqrt(x)

def v(tau1, tau2, z, x):
    f= 1
    return f**(5/2)/(sqrt(z)+(1-z)**3/2)*exp(3.57066164/sqrt(x)*(tau1-tau2))

def s(tau, taun, x):
    return exp(3.57066164/sqrt(x)*(tau-taun))


#Defining some variables
tau0=0.0
taumax=0.15
steps = 1000

#Empty arrays for storing x-values with corresponding weight
xsol= np.array([])
w=np.array([])

#Function for analytical solution
def D_a(x):
    return taumax/(sqrt(x)*(1-x)**(3/2))*exp(-np.pi*taumax**2/(1-x))


#Computing gluon evolution and repeating n=steps times
for j in range(1,steps+1):
    R=random()
    #Setting the initial x-value to 1; the v function value for the first step to 1 
    wnum=1
    x0 = 1
    #Arrays for storing x, z, tau values for each computation
    z = [0 for a in range(10000)]
    taustr=[0 for a in range(10000)]
    x = [0 for a in range(10000)]

    #Calculation of tau1
    if (math.isclose(R, 0.0) == False):
        tau1=tau0-np.log(R)/(phi(x0))
        taustr[0]=tau0
        taustr[1]=tau1
        x[0]=x0
        
        #Step 1
        
        if tau1>taumax:
            xsol = np.append(xsol, [x0]) #store the sx value in an array
            w=np.append(w, s(taumax, taustr[0], x0)) #store the weight for the first step
            
            print(j, '\t', taustr[1])
            continue
        #Steps 2,3, ..., n+1  
        else:
            i=0
            while (taumax>=taustr[i+1]):
                i=i+1
                
                R1=random()
                R2=random()
                
                #Generating the z value
                if R1<=kap1/kap:
                    znum=float(1-epsilon/(R2+(1-R2)*sqrt(epsilon))**2)
                if R1>kap1/kap:
                    znum=float(R2**2*(1-epsilon))
              
                z[i]= znum
                x[i]=x[i-1]*z[i] #calculation of the x
                if (x[i]<1e-4):
                  break
                wnum=wnum*v(taustr[i], taustr[i-1], z[i],x[i-1]) #calculation of the v-function product
            
                taustr[i+1]=taustr[i]-np.log(random())/phi(x[i]) #generating the next tau value
                if (x[i]>=1e-4) & (taustr[i+1]>taumax) & (x[i]<0.99):
                  xsol = np.append(xsol, x[i]) #store x values in an array 
            
                  w=np.append(w, wnum*s(taumax, taustr[i], x[i])) #store the value for the weight
    
                  print(j, '\t', taustr[i+1], '\t' , x[i], '\t', wnum*s(taumax, taustr[i], x[i]))
                  break 

bins = 50  #number of bins for the histogram
                          
w = w/(steps/bins) #normalisation of the histogram

#Plotting the analytical solution
x_analytical = np.linspace(0.0001, 0.99, 100)
D_analytical = np.array([])
for i in x_analytical:
  D_analytical =np.append(D_analytical, D_a(i))
print(D_analytical)

               

#Plotting the histogram of the evolution equation
fig,ax=plt.subplots(figsize=(12,9))

#Storing the D(x, tau) values in a txt file
n, bins2, patches = plt.hist(xsol, bins, weights=w)
with open('gluon_static_values_0.15tau.txt', 'w') as file1:
    content = ['x', '\t', 'D(x,tau)', '\n']
    file1.writelines(content)
    for i in range(0, len(n)):
      print(bins2[i], '\t', n[i], file = file1)
      
    file1.close()


plt.hist(xsol, bins=bins,weights=w, color='blue', label='numerical solution')
plt.plot(x_analytical, D_analytical, color= 'red', label ='analytical solution')
pyplot.yscale('log')
pyplot.ylim(1e-4, 10e3)
pyplot.ylabel('D(x,t)')
pyplot.xlabel('x')
pyplot.title('Gluon evolution in QGP for static medium')
pyplot.legend()
ax.set_xlim([0,1.1])
textstr = '\n'.join((
    r'$\tau_{0}=%.4f$' % (tau0, ),
    r'$\tau_{max}=%.3f$' % (taumax, )))

pyplot.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=14,
        verticalalignment='top')

plt.grid()
pyplot.yscale('log')
pyplot.savefig('gluon_static_0.15tau.jpg', dpi=None, facecolor='w', edgecolor='w', orientation='portrait', papertype=None, format='jpeg', transparent=False, bbox_inches=None, pad_inches=0.1, frameon=None, metadata=None)
