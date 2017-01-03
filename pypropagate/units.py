
import expresso.pycas as pc

metric_prefixes = [
     ('Y', 24, '\\mathrm{Y}'),
     ('Z', 21, '\\mathrm{Z}'),
     ('E', 18, '\\mathrm{E}'),
     ('P', 15, '\\mathrm{P}'),
     ('T', 12, '\\mathrm{T}'),
     ('G', 9, '\\mathrm{G}'),
     ('M', 6, '\\mathrm{M}'),
     ('k', 3, '\\mathrm{k}'),
#     ('h', 2, '\\mathrm{h}'),
#     ('da', 1, '\\mathrm{da}'),
     ('',0,''),
#     ('d', -1, '\\mathrm{d}'),
#     ('c', -2, '\\mathrm{c}'),
     ('m', -3, '\\mathrm{m}'),
     ('u', -6, '\\mu'),
     ('n', -9, '\\mathrm{n}'),
     ('p', -12, '\\mathrm{p}'),
     ('f', -15, '\\mathrm{f}'),
     ('a ', -18, '\\mathrm{a }'),
     ('z', -21, '\\mathrm{z}'),
     ('y', -24, '\\mathrm{y}')
]

base_units = set()

def add_metric_prefixes(name):
    s = globals()[name]
    for p,v,l in metric_prefixes:
        globals()[p+name] = 10**v*s

def create_unit(name):
    u = pc.Symbol("SI base unit "+name, type=pc.Types.Real, positive=True ,latex=r'\mathrm{%s}' % name,repr = name)
    globals()[name] = u
    base_units.add(u)
    add_metric_prefixes(name)

create_unit('m')
create_unit('kg')
create_unit('s')
create_unit('A')
create_unit('K')
create_unit('mol')
create_unit('cd')

C = s*A
add_metric_prefixes('C')

N = kg*m/s**2
add_metric_prefixes('N')

Hz = 1/s
add_metric_prefixes('Hz')

Pa = N/m**2
add_metric_prefixes('Pa')

V = kg*m**2*s**-3*A**-1
add_metric_prefixes('V')

J = C*V
add_metric_prefixes('J')

W = J/s
add_metric_prefixes('W')

F = kg**-1*m**-2*s**4*A**2
add_metric_prefixes('F')

ohm = kg*m**2*s**-3*A**-2
add_metric_prefixes('ohm')

eV = 1.6021766208 * pc.S(10)**-19 * J
add_metric_prefixes('eV')

h = 6.626070040*pc.S(10)**-34*J*s
hbar = h/(2*pc.pi)

c = 299792458 * m / s

degrees = pc.pi/180

def contains_unit(expr):
    for e in pc.postorder_traversal(expr):
        if e in base_units:
            return True
    return False

def get_unit(expr,only_base_units = False,evaluate=True,cache = None):
    res = None
    if expr.is_symbol:
        #if only_base_units:
            if expr in base_units:
                res = expr
        #else:
        #    if pc.Type(expr).evaluate(cache = cache) == pc.Types.Unit:
        #        res = expr
    elif expr.function == pc.multiplication:
        units = []
        for arg in expr.args:
            u = get_unit(arg,only_base_units,False,cache)
            if u is not None:
                units.append(u)
        if len(units) != 0:
            res = pc.multiplication(*units)
    elif pc.negative == expr.function:
        res = get_unit(expr.args[0],only_base_units,False,cache)
    elif pc.fraction == expr.function:
        inner_unit = get_unit(expr.args[0],only_base_units,False,cache)
        if inner_unit is not None:
            res = pc.fraction(inner_unit)
    elif pc.exponentiation == expr.function:
        inner_unit = get_unit(expr.args[0],only_base_units,False,cache)
        if inner_unit is not None:
            res = pc.exponentiation(inner_unit,expr.args[1])
    if res is not None and evaluate:
        res = res.evaluate(cache = cache)
    return res



