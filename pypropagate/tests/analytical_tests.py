from unittest import TestCase
from pypropagate import *

def analytical_slab_waveguide(settings):
    from pypropagate.coordinate_ndarray import CoordinateNDArray

    import scipy.optimize
    import scipy.integrate
    import numpy as np
    import warnings
    from numpy import sqrt,ceil,pi,sign,cos,sin,exp,tan,sum,abs

    s = settings.symbols

    k = settings.get_as(s.k,float)
    n1 = settings.get_as(settings.waveguide.n_1,complex).conjugate()
    n2 = settings.get_as(settings.waveguide.n_2,complex).conjugate()
    d  = 2*settings.get_as(settings.waveguide.r,float)

    # We will determine um for all guided modes as the roots of the characteristic equation:
    def Char(kappa):
        gamma = sqrt((n1.real**2-n2.real**2)*k**2-kappa**2)
        #return tan(kappa*d) - (2*kappa*gamma/(kappa**2-gamma**2))
        return tan(kappa*d/2) - gamma/kappa

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
                    kappa,r = scipy.optimize.brentq(Char,interval[0],interval[1],full_output=True)
                    # Check that point converged and value is small (we might have converged to a pole instead)
                    if kappa!=0 and r.converged and abs(Char(kappa))<1:
                        kappa_values.append(kappa)

    kappa_values = np.array(kappa_values)

    # Define the guided modes
    def psi(kappa,r):
        gamma = sqrt((n1.real**2-n2.real**2)*k**2-kappa**2)
        return cos(kappa*r)*(r*2<=d) + cos(kappa*d/2) * exp(-gamma*(r-d/2)) * (r*2>d)

    # Normalize the modes by integrating the intensity of the guided modes over R
    B_values = [(2*scipy.integrate.quad(lambda r:abs(psi(kappa,r))**2,0,np.inf)[0])**(-0.5) for kappa in kappa_values]

    #Project the modes onto the plain wave initial condition
    c_values = [2*B**2*scipy.integrate.quad(lambda r:psi(kappa,r),0,np.inf)[0] for kappa,B in zip(kappa_values,B_values)]

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
            mode = c * psi(kappa,r) * exp((-mu+1j*(beta-k))*x)
            if solution is None: solution = mode
            else: solution += mode
        return solution

    x_values = np.linspace(*settings.get_as((s.zmin,s.zmax,s.Nz),float))
    r_values = abs(np.linspace(*settings.get_as((s.xmin,s.xmax,s.Nx),float)))

    data = np.conjugate(field(*np.meshgrid(x_values,r_values)) )

    sx = settings.get_numeric(s.sx)
    sz = settings.get_numeric(s.sz)

    res = CoordinateNDArray(data,[(-sx/2,sx/2),(0,sz)],(s.x,s.z),settings.get_numeric_transform())
    return res

def analytical_circular_waveguide(settings):

    from pypropagate.coordinate_ndarray import CoordinateNDArray

    import scipy.optimize
    import scipy.integrate
    import warnings
    from scipy.special import j0,j1,k0,k1
    import numpy as np
    from numpy import sqrt,ceil,pi,sign,cos,sin,exp,sum,abs

    s = settings.symbols

    kn = settings.get_as(s.k,float)
    n1 = settings.get_as(settings.waveguide.n_1,complex).conjugate()
    n2 = settings.get_as(settings.waveguide.n_2,complex).conjugate()
    a  = settings.get_as(settings.waveguide.r,float)

    # We will determine um for all guided modes as the roots of the characteristic equation:
    def char(um):
        wm = sqrt(V**2-um**2)
        return um*j1(um)/j0(um) - wm*k1(wm)/k0(wm)

    # Dimensionless waveguide paramter
    V = kn*a*sqrt(n1.real**2 - n2.real**2)

    # V determines the number of guided modes
    N = ceil(V/pi)

    # Now find all N roots of the characteristic equation
    u_values = []
    segments = N

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

        while len(u_values) < N:
            segments = segments*2
            u_values = []
            xsegs = np.linspace(0,V,segments)
            for interval in [(xi,xj) for xi,xj in zip(xsegs[1:],xsegs[:-1])]:
                va,vb = char(np.array(interval))
                # Check if values for interval boundaries are finite and have opposite sign
                if np.isfinite(va) and np.isfinite(vb) and sign(va) != sign(vb):
                    # There might be a root in the interval. Find it using Brent's method!
                    um,r = scipy.optimize.brentq(char,interval[0],interval[1],full_output=True)
                    # Check that point converged and value is small (we might have converged to a pole instead)
                    if um!=0 and r.converged and abs(char(um))<1:
                        u_values.append(um)
    u_values = np.array(u_values)

    # Reactivate warnings
    warnings.resetwarnings()

    # Define the guided modes
    def psi(um,r):
        wm = sqrt(V**2-um**2)
        return np.piecewise(r,[r<a,r>=a],[lambda r:j0(um*r/a)/j0(um),lambda r:k0(wm*r/a)/k0(wm)])

    # Normalize the modes by integrating the intensity of the guided modes over R^2
    B_values = [(2*pi*scipy.integrate.quad(lambda r:abs(psi(um,r))**2*r,0,np.inf)[0])**(-0.5) for um in u_values]

    # Project the modes onto the plain wave initial condition
    c_values = [B**2*2*pi*scipy.integrate.quad(lambda r:psi(um,r)*r,0,np.inf)[0] for um,B in zip(u_values,B_values)]

    # Integrate mu in two parts to aviod discontinuity
    mu_values = np.array([
            2*pi*scipy.integrate.quad(lambda r:abs(psi(um,r)*B)**2*kn*n1.imag*r,0,a)[0] +
            2*pi*scipy.integrate.quad(lambda r:abs(psi(um,r)*B)**2*kn*n2.imag*r,a,np.inf)[0]
                         for um,B in zip(u_values,B_values)])

    beta_values = sqrt(n1.real**2*kn**2-u_values**2/a**2)

    # The full solution is the superposition of the guided modes
    def field(x,r):
        solution = None
        for i in range(len(u_values)):
            mode = c_values[i] * psi(u_values[i],np.abs(r)) * exp( (1j*(beta_values[i]-kn)-mu_values[i]) * x)
            if solution is None: solution = mode
            else: solution += mode
        return solution

    x_values = np.linspace(*settings.get_as((s.zmin,s.zmax,s.Nz),float))
    r_values = np.linspace(*settings.get_as((s.xmin,s.xmax,s.Nx),float))

    data = np.conjugate(field(*np.meshgrid(x_values,r_values)))

    sx = settings.get_numeric(s.sx)
    sz = settings.get_numeric(s.sz)

    res = CoordinateNDArray(data,[(-sx/2,sx/2),(0,sz)],(s.x,s.z),settings.get_numeric_transform())
    return res


class TestWaveguide(TestCase):

    def __init__(self, *args, **kwargs):
        super(TestWaveguide, self).__init__(*args, **kwargs)

        settings = presets.create_paraxial_wave_equation_settings()

        presets.set_plane_wave_initial_conditions(settings)
        s = settings.symbols

        wg = settings.create_category('waveguide')
        wg.create_symbol('n_1')
        wg.create_symbol('n_2')
        wg.create_symbol('r')

        s.n = pc.piecewise((wg.n_1,s.x**2+s.y**2<=wg.r**2),(wg.n_2,True))

        settings.wave_equation.set_energy(12*units.keV)
        settings.simulation_box.set((0.2*units.um,0.2*units.um,0.8*units.mm),(1024,1025,1024))

        wg.n_2 = presets.create_material('Ge',settings)
        wg.n_1 = 1
        wg.r = 50 * units.nm

        self.settings = settings

    def test_fresnel_2D(self):
        settings = self.settings
        s = settings.symbols
        propagator = propagators.FresnelPropagator1D(settings)
        an_field = analytical_slab_waveguide(settings)
        fd_field = propagator.run_slice()[:]
        self.assertLess(abs(an_field - fd_field)[:,s.sz/2:].max()/abs(an_field[:,s.sz/2:]).max(),0.1)

    def test_fresnel_3D(self):
        settings = self.settings
        s = settings.symbols
        propagator = propagators.FresnelPropagator2D(settings)
        an_field = analytical_circular_waveguide(settings)
        fd_field = propagator.run_slice()[:,0]
        self.assertLess(abs(an_field - fd_field)[:,s.sz/2:].max()/abs(an_field[:,s.sz/2:]).max(),0.2)

    def test_fresnel_RS(self):
        settings = self.settings.copy()
        s = settings.symbols
        wg = settings.waveguide
        s.u0 = pc.piecewise((1, abs(s.x) < wg.r), (1 - (abs(s.x) - wg.r) / (s.sx / 2 - wg.r), True))
        settings.simulation_box.set((5 * wg.r, 5 * wg.r, 0.8 * units.mm), (2048, 2048, 2048))

        propagator = propagators.FresnelPropagatorRS(settings)
        an_field = analytical_circular_waveguide(settings)
        fd_field = propagator.run_slice()[:]
        self.assertLess(abs(an_field - fd_field)[:,s.sz/2:].max()/abs(an_field[:,s.sz/2:]).max(),0.2)

    def test_finite_differences_2D(self):
        settings = self.settings
        s = settings.symbols
        propagator = propagators.FiniteDifferencesPropagator1D(settings)
        an_field = analytical_slab_waveguide(settings)
        fd_field = propagator.run_slice()[:]
        self.assertLess(abs(an_field - fd_field)[:,s.sz/2:].max()/abs(an_field[:,s.sz/2:]).max(),0.1)

    def test_finite_differences_3D(self):
        settings = self.settings
        s = settings.symbols
        propagator = propagators.FiniteDifferencesPropagator2D(settings)
        an_field = analytical_circular_waveguide(settings)
        fd_field = propagator.run_slice()[:,0]
        self.assertLess(abs(an_field - fd_field)[:,s.sz/2:].max()/abs(an_field[:,s.sz/2:]).max(),0.1)

    def test_finite_differences_RS(self):
        settings = self.settings
        s = settings.symbols
        propagator = propagators.FiniteDifferencesPropagatorRS(settings)
        an_field = analytical_circular_waveguide(settings)
        fd_field = propagator.run_slice()[:]
        self.assertLess(abs(an_field - fd_field)[:, s.sz / 2:].max() / abs(an_field[:, s.sz / 2:]).max(), 0.1)


class TestGaussian(TestCase):

    def test_gaussian_3D(self):
        settings = presets.create_paraxial_wave_equation_settings()
        s = settings.symbols
        pde = settings.partial_differential_equation

        g = settings.create_category('gaussian',info='Parameters of the gaussian beam')
        g.create_symbol('w_0',info = 'Waist size at z = 0',type=pc.Types.Real,positive=True)
        g.create_function('w_r',(s.z,),pc.sqrt(2)*pc.sqrt((g.w_0)**2/2+2*pde.A*s.z),info = 'Waist size')
        g.create_function('u',(s.x,s.y,s.z),(g.w_0)**2/g.w_r**2*pc.exp(pde.F*s.z-(s.x**2+s.y**2)/(g.w_r**2)))

        g.w_0 = 0.4*units.um

        s.n = 1
        s.u0 = g.u
        s.u_boundary = g.u

        settings.wave_equation.set_energy(12*units.keV)
        settings.simulation_box.set((2*units.um,1*units.um,20*units.mm),(500,100,1000))

        fieldan = expression_to_array(g.u.subs(s.y,0),settings)

        propagator = propagators.FiniteDifferencesPropagator2D(settings)
        fieldfd = propagator.run_slice()[:,0]

        self.assertLess(abs(fieldan - fieldfd).max(),0.01)


