# -*- coding: utf-8 -*-
"""
Created on Fri May 17 05:15:11 2019

@author: Buu Duc
"""
import numpy as np
from sympy import solve
from sympy.abc import x, y, z
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
def matrix(a,b,c,d,e,f,g,h,i):
    return np.array([[a,b,c],[d,e,f],[g,h,i]])
class ungsuat:
    def __init__(self,*args):
        if ( len(args)==1):
            self.a= args[0]
        else:
            bd=args[0]
            E=float(args[1])
            v=float(args[2])
            G=E/(2*(1+v))
            K=(1/E)*np.array([[1,-v,-v],[-v,1,-v],[-v,-v,1]])
            ax,ay,az=np.linalg.solve(K,[bd[0,0],bd[1,1],bd[2,2]])
            axy=bd[0,1]*2*G
            axz=bd[0,2]*2*G
            ayz=bd[1,2]*2*G
            self.a=np.array([[ax,axy,axz],[axy,ay,ayz],[axz,ayz,az]])
    def get_us(self):
        return self.a
    def invariant(self):
        t1= np.linalg.det(np.array([[self.a[0,0],self.a[0,1]],[self.a[1,0],self.a[1,1]]]))
        t2= np.linalg.det(np.array([[self.a[1,1],self.a[1,2]],[self.a[2,1],self.a[2,2]]]))
        t3= np.linalg.det(np.array([[self.a[0,0],self.a[0,2]],[self.a[2,0],self.a[2,2]]]))
        I1=self.a[1,1]+self.a[2,2]+self.a[0,0]
        I2=t1+t2+t3
        I3=np.linalg.det(self.a)
        return I1,I2,I3
    def uschinh(self):
        return np.sort(np.roots([1,-self.invariant()[0],+self.invariant()[1],-self.invariant()[2]]))[: : -1]
    def phuongchinh(self,i):
        i=i-1
        f=self.a-np.diag([self.uschinh()[i],self.uschinh()[i],self.uschinh()[i]])
        print(f)
        k=solve([f[0,0]*x+f[0,1]*y+f[0,2]*z,f[1,0]*x+f[1,1]*y+f[1,2]*z,x**2+y**2+z**2-1],dict=True)
        return k
    def usthuytinh(self):
        return self.invariant()[0]/3
    def uscau(self):
        k=self.usthuytinh()
        return np.diag([k,k,k])
    def uslech(self):
        return self.a-self.uscau()
    def mohr(self):
        a,b,c=self.uschinh()
        A=abs(a)
        C=abs(c)
        fig, ax = plt.subplots()
        ax.set_xlim((-(A+C),(A+C)))
        ax.set_ylim((-(A+C),(A+C)))
        ax.add_artist(Circle(((a+c)/2,0),(a-c)/2))
        ax.add_artist(Circle(((b+c)/2,0),(b-c)/2,color="Red"))
        ax.add_artist(Circle(((a+b)/2,0),(a-b)/2,color="Black"))
        plt.show()
    
class biendang:
    def __init__(self,*args):
        self.args=args
        self.handle()
    def handle(self):
        if ( len(self.args)==1):
            self.a= self.args[0]
        else:
            us=self.args[0]
            E=float(self.args[1])
            v=float(self.args[2])
            G=E/(2*(1+v))
            ax=1/E*(us[0,0]-v*(us[1,1]+us[2,2]))
            ay=1/E*(us[1,1]-v*(us[0,0]+us[2,2]))
            az=1/E*(us[2,2]-v*(us[0,0]+us[1,1]))
            axy=us[0,1]/(2*G)
            axz=us[0,2]/(2*G)
            ayz=us[1,2]/(2*G)
            self.a=np.array([[ax,axy,axz],[axy,ay,ayz],[axz,ayz,az]])
    def get_bd(self):
        return self.a
    def invariant(self):
        t1= np.linalg.det(np.array([[self.a[0,0],self.a[0,1]],[self.a[1,0],self.a[1,1]]]))
        t2= np.linalg.det(np.array([[self.a[1,1],self.a[1,2]],[self.a[2,1],self.a[2,2]]]))
        t3= np.linalg.det(np.array([[self.a[0,0],self.a[0,2]],[self.a[2,0],self.a[2,2]]]))
        I1=self.a[1,1]+self.a[2,2]+self.a[0,0]
        I2=t1+t2+t3
        I3=np.linalg.det(self.a)
        return I1,I2,I3
    def bdchinh(self):
        return np.sort(np.roots([1,-self.invariant()[0],+self.invariant()[1],-self.invariant()[2]]))[: : -1]
    def phuongchinh(self,i):
        i=i-1
        f=self.a-np.diag([self.bdchinh()[i],self.bdchinh()[i],self.bdchinh()[i]])
        k=solve([f[0,0]*x+f[0,1]*y+f[0,2]*z,f[1,0]*x+f[1,1]*y+f[1,2]*z,x**2+y**2+z**2-1],dict=True)
        return k
    def bdthuytinh(self):
        return self.invariant()[0]/3
    def bdcau(self):
        k=self.bdthuytinh()
        return np.diag([k,k,k])
    def bdlech(self):
        return self.a-self.bdcau()
    def mohr(self):
        a,b,c=self.bdchinh()
        A=abs(a)
        C=abs(c)
        fig, ax = plt.subplots()
        ax.set_xlim((-(A+C),(A+C)))
        ax.set_ylim((-(A+C),(A+C)))
        ax.add_artist(Circle(((a+c)/2,0),(a-c)/2))
        ax.add_artist(Circle(((b+c)/2,0),(b-c)/2,color="Red"))
        ax.add_artist(Circle(((a+b)/2,0),(a-b)/2,color="Black"))
        plt.Text()
        plt.show()
#print(g.phuongchinh(1))
#j=g.phuongchinh(2)
#print(j[0][x])
#print(g.uslech())
#g.mohr()
#print(a.invariant())
#a.mohr()
#g=biendang(a,4,5)
#g.mohr()
a=np.array([[4,-2,1],[-2,3,5],[1,5,-1]])
g=ungsuat(a)
k=g.phuongchinh(1)
print(k[0])
f=a=np.array([[1,4,2],[4,-1,-2],[2,-2,0]])
k=biendang(f,2*(10**6),0.2)
print(k.get_bd())