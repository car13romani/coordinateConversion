# executar como: python conversionKep2Cart.py terraElemOrbit.dat

import math
import sys
import re

inputFile = open(sys.argv[1])
outputFile = open('outPosVel.dat','w')

def kep2cart(a, e, i, w, Om, theta):
  a = a*149597870.7  # au para km
  mu = 132712440041.0  # GM do sol
  i = i * math.pi / 180.0 # degree -> rad
  w = w * math.pi / 180.0
  Om = Om * math.pi / 180.0
  theta = theta * math.pi / 180.0
  
  n = (mu / (a ** 3)) ** 0.5
  
  # calcular a anomalia excentrica
  EA = 2 * math.atan((((1.0-e)/(1.0+e)) ** 0.5) * math.tan(theta / 2.0)) 

  # calcular a anomalia media
  M = EA - e * math.sin(EA)
  
  # calcular o raio
  #r_mag = a * (1.0 - e * math.cos(EA))
  r_mag = (a*(1-(e**2)))/(1+(e*math.cos(theta)))
  # calcular o momento angular
  h_mag = ((1.0 - (e ** 2)) * (a * mu)) ** 0.5
  
  p = a * (1.0 - e ** 2)
  
  # calcular as componentes de posição
  x = r_mag * (math.cos(Om) * math.cos(w + theta) - math.sin(Om) * math.sin(w + theta) * math.cos(i))
  y = r_mag * (math.sin(Om) * math.cos(w + theta) + math.cos(Om) * math.sin(w + theta) * math.cos(i))
  z = r_mag * (math.sin(i) * math.sin(w + theta))
  
  # calcular as componentes de velocidade
  xdot = (x * h_mag * e / (r_mag * p)) * math.sin(theta) - (h_mag / r_mag) * (math.cos(Om) * math.sin(w + theta) + math.sin(Om) * math.cos(w + theta) * math.cos(i))
  ydot = (y * h_mag * e / (r_mag * p)) * math.sin(theta) - (h_mag / r_mag) * (math.sin(Om) * math.sin(w + theta) - math.cos(Om) * math.cos(w + theta) * math.cos(i))  
  zdot = (z * h_mag * e / (r_mag * p)) * math.sin(theta) + (h_mag / r_mag) * (math.sin(i) * math.cos(w + theta)) 
  
  # km/s -> km/dia
  xdot = xdot*86400
  ydot = ydot*86400
  zdot = zdot*86400

  return x, y, z, xdot, ydot, zdot

a_ = 0.0
e_ = 0.0
i_ = 0.0
w_ = 0.0
Om_ = 0.0
theta_ = 0.0
l1 = 'a'

i=1
outputFile.write("*******************************************************************************\n")
outputFile.write("  Symbol meaning [1 day=86400.0 s]:\n\n")
outputFile.write("    JDTDB    Julian Day Number, Barycentric Dynamical Time\n")
outputFile.write("      X      X-component of position vector (km)                               \n")
outputFile.write("      Y      Y-component of position vector (km)                               \n")
outputFile.write("      Z      Z-component of position vector (km)                               \n")
outputFile.write("      VX     X-component of velocity vector (km/day)                           \n")
outputFile.write("      VY     Y-component of velocity vector (km/day)                           \n")
outputFile.write("      VZ     Z-component of velocity vector (km/day)      \n")
outputFile.write("*******************************************************************************\n")

for line in inputFile:

  i+=1
  if i<20: 
    continue
  
  elif i%5 == 0:
    l1 = line.split("\n")
    l1 = l1[0]
  #  print(line)
  elif i%5 == 1:
    l2 = re.split("( EC= )|( QR= )|( IN= )|( EC=)|( QR=)|( IN=)|\n", line)
    l2 = list(filter(None, l2))
    e_ = float(l2[1])
    i_ = float(l2[5])
    
  elif i%5 == 2:
    l3 = re.split("( OM= )|( W = )|( Tp=  )|( OM=)|( W =)|( Tp= )|( Tp=)|\n", line)
    l3 = list(filter(None, l3))
    Om_ = float(l3[1])
    w_ = float(l3[3])

  elif i%5 == 3:
    l4 = re.split("( N = )|( MA= )|( TA= )|( N =)|( MA=)|( TA=)|\n", line)
    l4 = list(filter(None, l4))
    theta_ = float(l4[3])

  elif i%5 == 4:
    l5 = re.split("( A = )|( AD= )|( PR= )|( A =)|( AD=)|( PR=)|                          \n", line)
    l5 = list(filter(None, l5))
    a_ = (float(l5[1]))

    cart = kep2cart(a_,e_,i_,w_,Om_,theta_)

    outputFile.write(l1)
    outputFile.write("\n")
    outputFile.write(" X = %10.15E Y = %10.15E Z = %10.15E\n" % (cart[0], cart[1], cart[2]))
    outputFile.write(" VX= %10.15E VY= %10.15E VZ= %10.15E\n" % (cart[3], cart[4], cart[5]))

inputFile.close()
outputFile.close()

