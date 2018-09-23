def analytical_slab_waveguide(settings):
    from pypropagate.coordinate_ndarray import CoordinateNDArray
    import expresso.pycas as pc
    from pypropagate import expression_to_array
    import scipy.optimize
    import scipy.integrate
    from scipy.special import j0,j1,k0,k1
    import numpy as np
    from numpy import sqrt,ceil,pi,sign,cos,sin,exp,tan,sum,abs
    import warnings

    s = settings.symbols
    r = settings.waveguide.r
    k = float(settings.get_numeric(s.k*r))
    n1 = settings.get_as(settings.waveguide.n_1,complex)
    n2 = settings.get_as(settings.waveguide.n_2,complex)
    d  = 2
    dx = settings.get_as(s.dx/r,float)
    
    # We will determine um for all guided modes as the roots of the characteristic equation:
    def Char(kappa):
        gamma = sqrt((n1.real**2-n2.real**2)*k**2-kappa**2)
        return kappa*tan(kappa*d/2) - gamma
    kappa_max = sqrt(n1.real**2-n2.real**2)*k
        
    # maximum number of guided modes
    N = ceil(d*kappa_max/(2*pi))

    # Now find roots of the characteristic equation by searching intervals
    kappa_values = []
    segments = N
            
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        while len(kappa_values) < N:
            segments = segments*2
            if segments>10000:
                break
            kappa_values = []
            xsegs = np.linspace(0,kappa_max,segments)
            for interval in [(xi,xj) for xi,xj in zip(xsegs[1:],xsegs[:-1])]:
                va,vb = Char(np.array(interval))
                # Check if values for interval boundaries are finite and have opposite sign
                if np.isfinite(va) and np.isfinite(vb) and sign(va) != sign(vb):
                    # There might be a root in the interval. Find it using Brent's method!
                    kappa,res = scipy.optimize.brentq(Char,interval[0],interval[1],full_output=True)
                    # Check that point converged and value is small (we might have converged to a pole instead)
                    if kappa!=0 and res.converged and abs(Char(kappa))<1:
                        kappa_values.append(kappa)
    kappa_values = np.array(kappa_values)

    # Define the guided modes
    def psi(kappa,r):
        gamma = sqrt((n1.real**2-n2.real**2)*k**2-kappa**2)
        return cos(kappa*r)*(r*2<=d) + cos(kappa*d/2) * exp(-gamma*(r-d/2))*(r*2>d)
        
    # Normalize the modes by integrating the intensity of the guided modes over R
    B_values = [(2*scipy.integrate.quad(lambda r:abs(psi(kappa,r))**2,0,np.inf)[0])**(-0.5) for kappa in kappa_values]
    
    #Project the modes onto the symmetrical initial condition
    u0 = pc.numpyfy(settings.get_numeric(s.u0.subs(s.y,0).subs(s.x,s.x*r)),restype=float)
    c_values = [2*B**2*scipy.integrate.quad(lambda r:psi(kappa,r)*u0(x=r),0,np.inf)[0] for kappa,B in zip(kappa_values,B_values)]
    
    # Determine attenuation coeffcients
    mu_values = np.array([2*B**2*k*
            (scipy.integrate.quad(lambda r:abs(psi(kappa,r))**2*n1.imag,0,d/2)[0]   +
             scipy.integrate.quad(lambda r:abs(psi(kappa,r))**2*n2.imag,d/2,np.inf)[0])
                         for kappa,B in zip(kappa_values,B_values)])

    # The full solution is the superposition of the guided modes
    def field(x,r):
        solution = None
        for c,b,kappa,mu in zip(c_values,B_values,kappa_values,mu_values):
            beta = sqrt(k**2*n1.real**2 - kappa**2)
            mode = c * psi(kappa,r) * exp((mu+1j*(beta-k))*x)
            if solution is None: solution = mode
            else: solution += mode
        return solution
    
    linspace = lambda a,b,N: np.linspace(a,b,int(N))
    x_values = linspace(*settings.get_as((s.zmin/r,s.zmax/r,s.Nz),float))
    r_values = abs(linspace(*settings.get_as((s.xmin/r,s.xmax/r,s.Nx),float)))
    data = np.conjugate(field(*np.meshgrid(x_values,r_values)) )
    sx = settings.get_numeric(s.sx)
    sz = settings.get_numeric(s.sz)
    res = CoordinateNDArray(data,[(-sx/2,sx/2),(0,sz)],(s.x,s.z),settings.get_numeric_transform())
    return res


def analytical_circular_waveguide(settings):
    from pypropagate.coordinate_ndarray import CoordinateNDArray
    import expresso.pycas as pc
    import warnings
    import scipy.optimize
    import scipy.integrate
    from scipy.special import j0,j1,k0,k1
    import numpy as np
    from numpy import sqrt,ceil,pi,sign,cos,sin,exp,sum,abs

    s = settings.symbols
    r = settings.waveguide.r
    dx = float(settings.get_numeric(s.dx/r))
    kn = float(settings.get_numeric(s.k*r))
    n1 = settings.get_as(settings.waveguide.n_1,complex)
    n2 = settings.get_as(settings.waveguide.n_2,complex)
    R  = 1
    
    # We will determine um for all guided modes as the roots of the characteristic equation:
    def char(kappa):
        gamma = sqrt(kn**2*(n1.real**2 - n2.real**2) - kappa**2)
        return gamma*k1(gamma*R)*j0(kappa*R) - kappa*j1(kappa*R)*k0(gamma*R)
    
    kappa_max = kn * sqrt(n1.real**2-n2.real**2)
    
    # Dimensionless waveguide paramter
    V = kn*R*sqrt(n1.real**2 - n2.real**2)
    
    # V determines the number of guided modes
    N = ceil(V/pi)
        
    # Now find all N roots of the characteristic equation
    kappa_values = []
    segments = N
        
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

        while len(kappa_values) < N and segments < 10000:                
            segments = segments*2
            kappa_values = []
            xsegs = np.linspace(0,kappa_max,segments)
            for interval in [(xi,xj) for xi,xj in zip(xsegs[1:],xsegs[:-1])]:
                va,vb = char(np.array(interval))
                # Check if values for interval boundaries are finite and have opposite sign
                if np.isfinite(va) and np.isfinite(vb) and sign(va) != sign(vb):
                    # There might be a root in the interval. Find it using Brent's method!
                    kappa,res = scipy.optimize.brentq(char,interval[0],interval[1],full_output=True)
                    # Check that point converged and value is small (we might have converged to a pole instead)
                    if kappa!=0 and res.converged and abs(char(kappa))<1:
                        kappa_values.append(kappa)
        # Reactivate warnings
        warnings.resetwarnings()

    kappa_values = np.array(kappa_values)

    # Define the guided modes
    def psi(kappa,r):
        gamma = sqrt(kn**2*(n1.real**2 - n2.real**2) - kappa**2)
        return np.piecewise(r,[r<R,r>=R],[lambda r:j0(kappa*r),lambda r:j0(kappa*R)*k0(gamma*r)/k0(gamma*R)])
   
    # Normalize the modes by integrating the intensity of the guided modes over R^2
    B_values = [(2*pi*scipy.integrate.quad(lambda r:abs(psi(kappa,r))**2*r,0,R)[0]
                + (2*pi*scipy.integrate.quad(lambda r:abs(psi(kappa,r))**2*r,R,np.inf)[0]))**(-0.5) for kappa in kappa_values]
    
    #Project the modes onto the symmetrical initial condition
    u0 = pc.numpyfy(settings.get_numeric(s.u0.subs(s.y,0).subs(s.x,s.x*r)),restype=float)
    c_values = [B**2*2*pi*(scipy.integrate.quad(lambda r:psi(kappa,r)*u0(x=r)*r,0,R)[0]
                +scipy.integrate.quad(lambda r:psi(kappa,r)*u0(x=r)*r,R,np.inf)[0]) for kappa,B in zip(kappa_values,B_values)]
        
    # Integrate mu
    mu_values = np.array([
            2*pi*scipy.integrate.quad(lambda r:abs(psi(kappa,r)*B)**2*kn*n1.imag*r,0,R)[0] +
            2*pi*scipy.integrate.quad(lambda r:abs(psi(kappa,r)*B)**2*kn*n2.imag*r,R,np.inf)[0]
                         for kappa,B in zip(kappa_values,B_values)])
    beta_values = sqrt(n1.real**2*kn**2-kappa_values**2)
    
    # The full solution is the superposition of the guided modes
    def field(x,r):
        solution = None
        for i in range(len(kappa_values)):
            mode = c_values[i] * psi(kappa_values[i],np.abs(r)) * exp( (1j*(beta_values[i]-kn)+mu_values[i]) * x)
            if solution is None: solution = mode
            else: solution += mode
        return solution
    linspace = lambda a,b,N: np.linspace(a,b,int(N))
    x_values = linspace(*settings.get_as((s.zmin/r,s.zmax/r,s.Nz),float))
    r_values = linspace(*settings.get_as((s.xmin/r,s.xmax/r,s.Nx),float))
    data = np.conjugate(field(*np.meshgrid(x_values,r_values)))
    sx = settings.get_numeric(s.sx)
    sz = settings.get_numeric(s.sz)
    res = CoordinateNDArray(data,[(-sx/2,sx/2),(0,sz)],(s.x,s.z),settings.get_numeric_transform())
    return res