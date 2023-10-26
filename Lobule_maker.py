#
import numpy as np
import matplotlib.pyplot as plt
#from matplotlib import bezier
import bezier 
a=np.array([[1,0],[0.5,1],[-0.5,1],[-1,0],[-0.5,-1],[0.5,-1], [1,0]])
print("vertices=",a)
def f(x):
    #r=np.sqrt(x)
    arr = np.array([0, 1])
    S =np.zeros(len(x)-1)
    r =np.zeros(len(x)-1)
    for i in np.arange(len(x)-1):
        for j in arr:
            S[i] += (x[i,j]-x[i+1,j])**2
            r[i]=np.sqrt(S[i])

    return r

def g(nodes):
    curve  = bezier.Curve(nodes, degree=2)
    s_vals = np.linspace(nodes[0,0],nodes[0,2],3)
    p = curve.evaluate_multi(s_vals)
    return p



if __name__ == "__main__":

    b=np.array([[1,0],[1,1.1],[0,1.3],[-1,1.1],[-1,0.2],[-1,-1.3],[0,-1.1]])
    
    print("f(b)=",f(b))
    print("--------------------------------------")
    #nodes = np.asfortranarray([[0.0, 0.625, 1.0],[0.0, 0.5, 0.5],])
    curve=np.zeros(2)
##    nodes1=np.copy(b)
##    nodes1=  nodes1+ 0.1
    rg=np.random.default_rng(1)
    rd_nodes =np.sort(rg.random((2,2,3)))
    print("rd_nodes=",rd_nodes)
    print("--------------------------------------")
    M=np.zeros((2,2,3))
    for i in np.arange(2):
         print("i=",i)
         print("rd_nodes[:,i,:]=",rd_nodes[:,i,:])
         print("--------------------------------------")
         d=g(rd_nodes[:,i,:])
         #curve[i]=bezier.Curve(rd_nodes[:,i,:], degree=2)
         print("d=",d)
         print("--------------------------------------")
         M[:,i,:]=d
         print("M[:,i,:]=",M[:,i,:])
         print("--------------------------------------")
         
print("M=",M)
print("--------------------------------------")        
##    s_vals=np.linspace(0.0,1.0,5)
##    p=curve.evaluate_multi(s_vals)
##    #x =a[:,0]
##    #y =a[:,1]
##    #plt.plot(x,y)
##    plt.plot(p[0,:],p[1,:])
##    plt.show()

    
  

      
        

 
