
def ice_and_fire_colormap(N = 1000,ratio = 10,cache={}):
    if (N,ratio) in cache:
        return cache[(N,ratio)]
    import numpy as np
    from matplotlib.colors import ListedColormap
    N1 = N
    N2 = (ratio-1)*N1
    C1 = [(0.8*(1-v**0.5),0.9*(1-v),1-v**2) for v in np.linspace(0,1,N1)]
    C2 = [(v**0.25,v,v**2) for v in np.linspace(0,1,N2)]
    C3 = C1 + C2
    cache[(N,ratio)] = ListedColormap(C3)
    return cache[(N,ratio)]

def ice_colormap(N = 100000,cache={}):
    if N in cache:
        return cache[N]
    import numpy as np
    from matplotlib.colors import ListedColormap
    C1 = [(0.8*(1-v**0.5),0.9*(1-v),1-v**2) for v in np.linspace(0,1,N)]
    cache[N] = ListedColormap(C1)
    return cache[N] 

def fire_colormap(N = 100000,hue_shift = 0,cache={}):
    key = (N,hue_shift)
    if key in cache:
        return cache[key]
    import numpy as np
    import colorsys
    from matplotlib.colors import ListedColormap
    C = [(v**0.25,v,v**2) for v in np.linspace(0,1,N)]
    if hue_shift != 0:
        C = [colorsys.rgb_to_hsv(c[0],c[1],c[2]) for c in C]
        C = [colorsys.hsv_to_rgb((c[0]+hue_shift) % 1,c[1],c[2]) for c in C]
    cache[key] = ListedColormap(C)
    return cache[key]                                                                                                                                                                            


