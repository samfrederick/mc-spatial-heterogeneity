import numpy as np

def SH(f):
    """Compute spatial heterogeneity
    
    Original SH module
    """
    m=np.mean(f)
    G=0
    d1=f.shape[0]
    d2=f.shape[1]
    f =np.hstack((f, f))
    f =np.vstack((f, f))

    for y1 in range(d1+1):
        for y2 in range(y1+1,y1+d1+1):
            for x1 in range(d2+1):
                for x2 in range(x1+1,x1+d2+1):
                    fbar=f[y1:y2,x1:x2]
                    G+= abs(np.mean(fbar)-m)

    return G/(d1*d1*d2*d2)