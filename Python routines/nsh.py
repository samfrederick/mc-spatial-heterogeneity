import numpy as np

def NSH(f):
  """Compute normalized spatial heterogeneity

  Combines performance speedups from refactored function `NSH_refactored()`
  and use of np.mean() in the style of the original `SH()` routine for subarray
  means instead of np.ix_() since np.ix_() is significantly more computationally 
  demanding.
  """
  m = np.mean(f)
  G = 0
  if len(f.shape) == 1:
    f = np.reshape(f, (-1, 1))
  d1 = f.shape[0]
  d2 = f.shape[1]
  N = m*(1.5*d1*d2*(d1-1)*(d2-1)+d1*(d1-1)+d2*(d2-1))/((d1*(d1-1)+1)*(d2*(d2-1)+1))
  if N == 0:
      N = 1
  f =np.hstack((f, f))
  f =np.vstack((f, f))

  for y1 in range(d2):
    for y2 in range(d2):
        if y1 == y2:
            continue
        if y2<y1:
            y2+=d2
        fbar = f[0:d1,y1:y2]
        G+= abs(np.mean(fbar)-m) 
        for x1 in range(d1):
            for x2 in range(d1):
                if x1 == x2:
                    continue
                if x2<x1:
                    x2+=d1
                fbar=f[x1:x2,y1:y2]
                G+= abs(np.mean(fbar)-m) 

  for x1 in range(d1):
    for x2 in range(d1):
        if x1 == x2:
            continue
        if x2<x1:
            x2+=d1
        fbar=f[x1:x2,0:d2]
        G+= abs(np.mean(fbar)-m) 

  return G/(N*(d1*(d1-1)+1)*(d2*(d2-1)+1))