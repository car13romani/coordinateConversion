# executar como: python conversionCart2Kep.py terraPosVel.dat

import math
import sys
import re

inputFile = open(sys.argv[1])
outputFile = open('outElemOrbit.dat','w')

def cart2kep(x,y,z,xdot,ydot,zdot,t):
  mu = 132712440041.0 
  r = [x,y,z]
  v = [xdot/86400,ydot/86400,zdot/86400]   #km/dia para km/segundo
  
  # calcular o momento angular
  h = [r[1] * v[2] - r[2] * v[1], r[2] * v[0] - r[0] * v[2], r[0] * v[1] - r[1] * v[0]]
  h_mag = (h[0] ** 2 + h[1] ** 2 + h[2] ** 2) ** (1/2.0)

  # calcular o raio
  r_mag = (r[0] ** 2 + r[1] ** 2 + r[2] ** 2) ** (1/2.0)

  # calcular a velocidade
  v_mag = (v[0] ** 2 + v[1] ** 2 + v[2] ** 2) ** (1/2.0)

  # calcular a energia especifica
  E = ((v_mag ** 2) / 2) - (mu / r_mag)

  # calcular o semi-eixo maior
  a = -mu / (2 * E)

  # calcular a excentricidade
  e = (1 - ((h_mag ** 2) / (a * mu))) ** 0.5

  # calcular a inclinação
  i = math.acos(h[2] / h_mag)

  #  Compute right ascension of the ascending node 
  Om = math.atan2(h[0], -h[1])

  p = a * (1 - (e ** 2))

  # calcula a true anomaly
  theta = math.atan2(((p / mu) ** 0.5) * (v[0] * r[0] + v[1] * r[1] + v[2] * r[2] ), p - r_mag)
  
  # calcula o argumento da periapse
  w = math.atan2(r[2] / math.sin(i), r[0] * math.cos(Om) + r[1] * math.sin(Om)) - theta
  
  # calcular a anomalia excentrica
  EA = 2*(math.atan((((1-e)/(1+e))**0.5)*math.tan(theta/2)))

  # Mean motion
  n = ((mu / a ** 3) ** 0.5)*86400
  
  # calcular a anomalia media
  MA = EA - e * math.sin(EA)

  # distancia periapse e apoapse
  QR = ((h_mag**2)/mu)*(1/(1+e))/149597870.7  # km -> AU
  AD = ((h_mag**2)/mu)*(1/(1-e))/149597870.7

  # periodo sideral
  PR = 2*math.pi*(((a**3)/mu)**0.5)/86400

  # mean motion
  N = ((2*180.0)/PR)

  #Tp is time of periaptic passage, t is time step
  Tp = t - (1/n) * (EA - e * math.sin(EA))
  i = i * 180 / math.pi
  w = w * 180 / math.pi
  Om = Om * 180 / math.pi
  theta = theta * 180 / math.pi  #true anomaly ou mean_anomaly
  a = a/149597870.7

  # angulos positivos
  if w<0:
    w = 360.0+w
  if i<0:
    i = 360.0+i
  if Om<0:
    Om = 360.0+Om
  if theta<0:
    theta = 360.0+theta
  if MA<0:
    MA = 360.0+MA


  return e, QR, i, Om, w, Tp, N, MA, theta, a, AD, PR

x_ = 0
y_ = 0
z_ = 0
vx_ = 0
vy_ = 0
vz_ = 0
t_ = 0
l1 = 'a'

i=1
outputFile.write("*******************************************************************************\n")
outputFile.write("  Symbol meaning [1 au= 149597870.700 km, 1 day= 86400.0 s]\n\n")
outputFile.write("    JDTDB    Julian Day Number, Barycentric Dynamical Time\n")
outputFile.write("      EC     Eccentricity, e                                                   \n")
outputFile.write("      QR     Periapsis distance, q (au)                                        \n")
outputFile.write("      IN     Inclination w.r.t XY-plane, i (degrees)                           \n")
outputFile.write("      OM     Longitude of Ascending Node, OMEGA, (degrees)                     \n")
outputFile.write("      W      Argument of Perifocus, w (degrees)                                \n")
outputFile.write("      Tp     Time of periapsis (Julian Day Number)                             \n")
outputFile.write("      N      Mean motion, n (degrees/day)                                      \n")
outputFile.write("      MA     Mean anomaly, M (degrees)                                         \n")
outputFile.write("      TA     True anomaly, nu (degrees)                                        \n")
outputFile.write("      A      Semi-major axis, a (au)                                           \n")
outputFile.write("      AD     Apoapsis distance (au)                                            \n")
outputFile.write("      PR     Sidereal orbit period (day)           \n")
outputFile.write("*******************************************************************************\n")
outputFile.write("$$SOE\n")


for line in inputFile:

  i+=1
  if i<13: 
    continue
    #print(line)
  elif i%3 == 1:
    l1 = line.split("\n")
    l1 = l1[0]
    t_aux = l1.split(" = ")
    t_ = float(t_aux[0])
    
    #print(l1)

  elif i%3 == 2:
    pos = re.split("( X = )|( Y = )|( Z = )|( X =)|( Y =)|( Z =)|\n", line)
    pos = list(filter(None, pos))
    x_ = float(pos[1])
    y_ = float(pos[3])
    z_ = float(pos[5])
    #print(line)
        
  elif i%3 == 0:
    vel = re.split("( VX= )|( VY= )|( VZ= )|( VX=)|( VY=)|( VZ=)|\n", line)
    vel = list(filter(None, vel))
    vx_ = (float(vel[1]))  
    vy_ = (float(vel[3]))
    vz_ = (float(vel[5]))
    kep = cart2kep(x_,y_,z_,vx_,vy_,vz_,t_)

    outputFile.write(l1)
    outputFile.write("\n")
    outputFile.write(" EC= %10.15E QR= %10.15E IN= %10.15E\n" % (kep[0], kep[1], kep[2]))
    outputFile.write(" OM= %10.15E W= %10.15E Tp= %10.15E\n" % (kep[3], kep[4], kep[5]))
    outputFile.write(" N = %10.15E MA= %10.15E TA= %10.15E\n" % (kep[6], kep[7], kep[8]))
    outputFile.write(" A = %10.15E AD= %10.15E PR= %10.15E\n" % (kep[9], kep[10], kep[11]))

inputFile.close()
outputFile.close()
