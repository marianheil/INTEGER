from __future__ import print_function
from math import sqrt

def is_int(a):
  return a%1==0

def m2(p):
  return -p[0]*p[0]-p[1]*p[1]-p[2]*p[2]+p[3]*p[3]

def m(p):
  return sqrt(m2(p))

def get_E(p,m=0):
  return sqrt(+p[0]*p[0]+p[1]*p[1]+p[2]*p[2]+m*m)

def possible_mom(mass,maxp,minpt):
  results=[]
  for i0 in xrange(-maxp,maxp):
    for i1 in xrange(-maxp,maxp):
      if(i0*i0+i1*i1>minpt*minpt):
        for i2 in xrange(-maxp,maxp):
          v=[i0,i1,i2]
          E=get_E(v,mass)
          if(is_int(E)):
            v.append(E)
            results.append(v)
  return results

p_parton=possible_mom(0,100,30)
print(len(p_parton))
p_higgs=possible_mom(125,100,0)
print(len(p_higgs))

results=[]
for a in p_parton:
  for d in p_higgs:
    t = a[0]+d[0]
    if(t>0):
      break
    if(t == 0 and a[1]+d[1] == 0):
      print("POSSIBLE FOUND ", a ,d)
      results.append([a,d])
# for d in p_higgs:
#   for a in p_parton:
#     for b in p_parton:
#       for c in p_parton:
#         t = a[0]+b[0]+c[0]+d[0]
#         if(t>0):
#           break
#         if(t == 0 and a[1]+b[1]+c[1]+d[1] == 0):
#           print("POSSIBLE FOUND ", a,b,c,d)
#           results.append([a,b,c,d])

print(len(results))
