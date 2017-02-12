

def phase_gradient(array):
    import numpy as np
    from _pypropagate import ring_derivative_2D
    dx = array.soft_copy()
    dy = array.soft_copy()
    dx.data = np.zeros(array.shape,dtype=np.float)
    dy.data = np.zeros(array.shape,dtype=np.float)
    ring_derivative_2D(np.angle(array.data),dy.data,dx.data,2*np.pi)
    return (dx,dy)

