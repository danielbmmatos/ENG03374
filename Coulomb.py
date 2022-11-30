#-----------------------------------------------------------------------------#
#                      Amortecimento por atrito                               #
#                                                                             #
#-----------------------------------------------------------------------------#
# Author: Eng. Me. Daniel B M Matos
# Date  : 30/11/2022
# Last Modification: 18/10/2022
#-----------------------------------------------------------------------------#
import numpy as np
import matplotlib.pyplot as plt

def Coulomb(xo,k,m,mu):
    
    from mpl_toolkits.axisartist.axislines import SubplotZero


    fig = plt.figure(1,figsize = (9,5))
    ax  = SubplotZero(fig, 111)
    fig.add_subplot(ax)

    for direction in ["xzero", "yzero"]:
    # adds arrows at the ends of each axis
        ax.axis[direction].set_axisline_style("-|>")

    # adds X and Y-axis from the origin
        ax.axis[direction].set_visible(True)

    for direction in ["left", "right", "bottom", "top"]:
    # hides borders
        ax.axis[direction].set_visible(False)
        
    g = 9.81
    xo1 = xo
    #Definindo a frequência natural e o periodo
    
    w = np.sqrt(k/m)
    T = 2*np.pi/w
    
    #Definido o número de ciclos
    
    r = (xo - mu*m*g/k)/(2*mu*m*g/k)
    
    # Loop para o calculo dos deslocamentos

    n = (int(r)+1)   # numero de meio ciclos
    nn = 100          # numero de pontos em cada meio ciclo
    t = np.linspace(0,n*T/2,n*nn)
    d = np.zeros(len(t))
    d[0] = xo
    mul  = np.zeros(len(t))
    for i in range (0,int(n/2)+5,2):
        d[i*nn+1:(i+1)*nn+1] = (xo - mu*m*g/k)*np.cos(w*t[i*nn+1:(i+1)*nn+1]) + mu*m*g/k
        x1 = d[(i+1)*nn]
        mul[i*nn+1:(i+1)*nn+1] = mu*m*g/k
        d[(i+1)*nn:(i+2)*nn+1] = -(x1 + mu*m*g/k)*np.cos(w*t[(i+1)*nn:(i+2)*nn+1]) - mu*m*g/k
        xo = d[(i+2)*nn]
        mul[(i+1)*nn:(i+2)*nn+1] =- mu*m*g/k
            
    ax.plot(t,d,'k',label ='x(t)')
    ax.plot(t,mul*2,'r:',label ='Componente de atrito')
    ax.plot(t,(xo1-2*w*mu*m*g*t/k/np.pi),'k--',label ='Decremento linear')
    ax.plot(t,-(xo1-2*w*mu*m*g*t/k/np.pi),'k--')
    plt.xlim(0,r*T/4)
    plt.xlabel('Tempo (s)') ; plt.ylabel('x (m)')
    plt.legend()
    plt.savefig('teste.pdf')
    