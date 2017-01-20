# -*= coding: utf-8 -*-
import sys
from PyQt4 import QtCore , QtGui
import numpy as np
from math import *
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4 import NavigationToolbar2QT as NavigationToolbar
from ui import plot

def func_T(x,T):  
    beta=7.021e-4 #eV/K
    gamma=1108 #K
    y=x-((beta*T**2)/(T+gamma))
    return y
    
def absorp_coeff(L,T): #L[m],T[K]
    c=3.0e8      #m/s
    k=8.617e-5   #eV/K
    h=4.136e-15  #eVs
        
    Eg01=1.15578 #eV
    Eg02=2.5 #eV
    Egd=3.2 #eV
    Ep1=1.827e-2 #eV
    Ep2=5.773e-2 #eV
    Ep=np.array([Ep1,Ep2])

    A11=1.777e5 #m^-1/eV^2
    A12=3.98e6 #m^-1/eV^2
    A21=1.292e5 #m^-1/eV^2
    A22=2.895e6 #m^-1/eV^2
    Ad=1.052e8 #m^-1/eV^2

    A=np.array([[A11,A12],[A21,A22]])

    Eg1=func_T(Eg01,T)
    Eg2=func_T(Eg02,T)
    Egd=func_T(Egd,T)
    Eg=np.array([Eg1,Eg2])

    f=c/L
    b=h*f
    a=b-Egd
    if a>0:
        alpha=Ad*sqrt(a)
        for i in [0,1]:
            for j in [0,1]:
                alpha=alpha+A[i,j]*((h*f-Eg[j]+Ep[i])**2/(exp(Ep[i]/k/T)-1)+(h*f-Eg[j]-Ep[i])**2/(1-exp(-Ep[i]/k/T)))
    else:
        alpha=0
        for i in [0,1]:
            for j in [0,1]:
                if Eg[j]+Ep[i]<b:
                    alpha=alpha+A[i,j]*((h*f-Eg[j]+Ep[i])**2/(exp(Ep[i]/k/T)-1)+(h*f-Eg[j]-Ep[i])**2/(1-exp(-Ep[i]/k/T)))
                elif Eg[j]-Ep[i]<b and Eg[j]+Ep[i]>b:
                    alpha=alpha+A[i,j]*((h*f-Eg[j]+Ep[i])**2/(exp(Ep[i]/k/T)-1))

                else: alpha=alpha

    return alpha


class Main(QtGui.QMainWindow,plot.Ui_MainWindow):
    def __init__(self):
        super(Main,self).__init__()
        self.setupUi(self)

    def addmpl(self,fig):
        self.canvas=FigureCanvas(fig)
        self.mplvl.addWidget(self.canvas)
        self.canvas.draw()
        self.toolbar=NavigationToolbar(self.canvas,self.widget,coordinates=True)
        self.mplvl.addWidget(self.toolbar)

        
class make_figure(Figure,plot.Ui_MaintWindow):
    def __init__(self,L,T):
        super(make_figure,self).__init__()
        self.alpha=absorp_coeff(L,T)
        
    def make_fig(self,alpha):
        x=np.arange(0,200e-6,1e-6)
        y=map(lambda z:exp(-self.alpha*z),x)
        self.add_subplot(111)
        self.plot(x,y)
        form=Main
        form.addmpl(self)
        
def main():
    app=QtGui.QApplication(sys.argv)

    Make=make_figure()
    
    form=Main()
    form.show()
    
    Calic=Caliculation(form.lineEdit.text(),form.lineEdit_2.text())

    form.pushButton.clicked.connect(Make.make_fig(Calic.absorp_coeff()))
    app.exec_()

if __name__=="__main__":
    main()
