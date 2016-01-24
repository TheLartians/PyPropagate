
# coding: utf-8

from __future__ import print_function, division

import ctypes
import tempfile
from subprocess import Popen, PIPE
import numpy as np

from sympy import *
from sympy.core import S
from sympy.core.compatibility import string_types
from sympy.printing.codeprinter import CodePrinter,Assignment
from sympy.printing.precedence import precedence


class c_complex(ctypes.Structure):
    """Class representing a 128 bit c++ complex type."""
    
    _fields_ = [('real',ctypes.c_double),('imag',ctypes.c_double)]
    def __init__(self,z,i=None):
        super(c_complex,self).__init__()
        if i:
            self.real = z
            self.imag = i
        else:
            z = complex(z)
            self.real = z.real
            self.imag = z.imag
    def __str__(self):
        return str(self.real) + '+' + str(self.imag) + 'j'
    def __repr__(self):
        return '(' + str(self.real) + ',' + str(self.imag) + ')'
    def __complex__(self):
        return complex(self.real,self.imag)

ctypes.c_complex = c_complex

def get_ctype(expr):
    if expr.is_Equality or expr.is_Relational or expr.is_Boolean: return ctypes.c_bool
    if expr.is_integer: return ctypes.c_int
    if expr.is_real: return ctypes.c_double
    return ctypes.c_complex

def get_ctype_name(type):
    
    type_names = {
        ctypes.c_bool:"bool",
        ctypes.c_int:"int",
        ctypes.c_double:"double",
        ctypes.c_complex:"complex",
        ctypes.c_bool:"bool",
        ctypes.c_bool:"bool",
        bool: "bool",
        int: "int64_t",
        float: "double",
        complex: "std::complex<double>",
        np.float64:'double',
        np.int64:'int',
        np.complex64:'std::complex<double>',
    }
    
    return type_names[type]
    
    if type == ctypes.c_bool: return "bool"
    if type == ctypes.c_int: return "int"
    if type == ctypes.c_double: return "double"
    if type == ctypes.c_complex: return "complex"
    
    if cls.array_type == bool: return "bool"
    if cls.array_type == int: return "int64_t"
    if cls.array_type == float: return "double"
    if cls.array_type == complex: return "std::complex<double>"

    raise ValueError("unknown type")

def get_evaluation_ctype_string(expr):
    if expr.is_Equality or expr.is_Relational or expr.is_Boolean: return "bool"
    if expr.is_integer: return "int"
    if expr.is_real: return "double"
    return "std::complex<double>"

class ArrayFunction(Function):
   
    @classmethod
    def eval(cls,x,y=None,z=None):
        if not cls.constant_array:
            return
        
        args = [x]
        if y != None:
            args.append(y)
        if z != None:
            args.append(z)            
        try:
            res = []
            for arg in args:
                res.append(int(arg))
            for i,arg in enumerate(args):
                if arg<0 or arg>=cls.array_data.shape[i]:
                    return 0
            return cls.array_data[tuple(res)]
        except TypeError:
            return
        
    @classmethod
    def get_ccode(cls):
        pointer, read_only_flag = cls.array_data.__array_interface__['data']

        if cls.array_type == bool: type = "bool"
        if cls.array_type == int: type = "int64_t"
        if cls.array_type == float: type = "double"
        if cls.array_type == complex: type = "std::complex<double>"

        if len(cls.array_data.shape) == 1:
            return "array_wrapper_1D<%s> data_%s(%s,%s);" % (type,cls.__name__,pointer,cls.array_data.shape[0])
        if len(cls.array_data.shape) == 2:
            return "array_wrapper_2D<%s> data_%s(%s,%s,%s);" % (type,cls.__name__,pointer,cls.array_data.shape[0],cls.array_data.shape[1])
        if len(cls.array_data.shape) == 3:
            return "array_wrapper_3D<%s> data_%s(%s,%s,%s,%s);" % (type,cls.__name__,pointer,cls.array_data.shape[0],cls.array_data.shape[1],cls.array_data.shape[2])
    
    def _eval_is_Boolean(self):
        return self.array_type == bool
    def _eval_is_integer(self):
        return self.array_type == int or self._eval_is_Boolean()
    def _eval_is_real(self):
        return self.array_type == float or self._eval_is_integer()
    def _eval_is_complex(self):
        return self.array_type == complex or self._eval_is_real()

class Array:
    """
    Sympy function representing an array. This allows numpy arrays to be used in sympy expressions.
    
    Members / Initializer Parameters
    --------------------------------
    name:   string
            Name of the array in the sympy expression. This name must be unique in the expression.
    data:   array-like
            Array to be accessed. At the moment only 1, 2 or 3 dimensional arrays are supported. Contents of this member may be modified after compilation (however, any operation that requires reallocation [such as resizing] will not influence the compiled expression).
            
    Notes
    -----
    The arguments are required be integer and are the access indices. Non-integer expressions may be rounded to integer e.g. using the sympy.floor operator. If the indices exceed the array boundaries the access will return 0.
    
    Usage Example
    -------------
    >>> array = Array("array",numpy.linspace(0,1,101))
    >>> x = sympy.Symbol("x",real=True)
    >>> f = array(sympy.floor(x*100))+x
    >>> f.subs(x,0.5)
    1.0
    """
    
    valid_types = { np.bool:bool,np.int:int,np.int32:int,np.int64:int,np.float64:float,np.complex128:complex }

    def __init__(self,name,data,is_constant = True):
        if name == "data_dummy": raise ValueError("Array name cannot be data_dummy")
        
        data_dummy = np.array(data)
        data_type = data_dummy.dtype
        
        if data_type == np.float:
            data_type = np.float64
        if data_type == np.bool:
            data_type = np.bool
        if data_type == np.int:
            data_type = np.int64
        if data_type == np.int32:
            data_type = np.int64
        if data_type == np.complex:
            data_type = np.complex128
        
        if not data_type in Array.valid_types:
            raise ValueError("Array dtype '%s' not accepted. dtype must be one of the following: %s" % (data_type,Array.valid_types.keys()))
        
        if len(data_dummy.shape) not in [1,2,3]:
            raise ValueError("Data has to be one or two or three dimensional.")
        
        name = name.lstrip()
        
        self.data = np.ascontiguousarray(data_dummy,dtype=data_type)
        self.name = name
        
        exec("class %s(ArrayFunction): nargs = %s" % (name,len(data_dummy.shape)))
        self.ArrayClass = eval(name)
        self.ArrayClass.array_type = self.valid_types[data_type]
        self.ArrayClass.array_data = self.data
        self.ArrayClass.constant_array = is_constant
        
    def __call__(self,x,y=None,z=None):
        x,y,z = (S(x),S(y),S(z))
        call_args = 1
        if y != None: call_args += 1
        if z != None: call_args += 1

        if call_args not in self.ArrayClass.nargs:
            raise ValueError("Array arguments do not match array dimension %s" % (self.ArrayClass.nargs))

        if not x.is_integer:
            raise ValueError("Non integer argument for array access")
        if y and not y.is_integer:
            raise ValueError("Non integer argument for array access")
        if z and not z.is_integer:
            raise ValueError("Non integer argument for array access")

        if call_args == 1: F = self.ArrayClass(x)
        if call_args == 2: F = self.ArrayClass(x,y)
        if call_args == 3: F = self.ArrayClass(x,y,z)

        return F        
    
class CCodeFunction(Function):
   
    @classmethod
    def eval(cls, *args):
        pass
    
    @classmethod
    def get_ccode(cls):
        return cls.__code__
    
class CCode:
    """
    Sympy function for arbitrary c++ code.
    
    Members / Initializer Parameters
    --------------------------------
    code:       string
                C++ code.
    call_name:  string
                Valid c++ identifier of function to be called with the given arguments. 
            
    Notes
    -----
    The code will be inserted into the global namespace. 
    No check is performed if the call arguments match the function signature.
    
    Usage Example
    -------------
    >>> x = sympy.Symbol('x',real=True)
    >>> cfunction = CCode("int cfunction(int x){ return 20*x; }","cfunction")
    >>> f = compile_sympy_expression(cfunction(x)+x,(x,),ctypes.c_double)
    >>> f(2)
    42.0
    """

    def __init__(self,code,call_name):
        self.call_name = call_name
        self.code = code
                
        # The following line includes a naive test to see if the array name is a valid identifier
        exec("class %s(CCodeFunction): pass" % (self.call_name))
        CCodeClass = eval(self.call_name)
        CCodeClass.__code__ = self.code
        self.CodeClass = CCodeClass
        
    def __call__(self,*args):
        return self.CodeClass(*args)

# dictionary mapping sympy function to (argument_conditions, C_function).
# Used in CCodePrinter._print_Function(self)
known_functions = {
    "Abs": "abs",
    "gamma": "tgamma",
    "sin": "sin",
    "cos": "cos",
    "tan": "tan",
    "asin": "asin",
    "acos": "acos",
    "atan": "atan",
    "atan2": "atan2",
    "exp": "exp",
    "log": "log",
    "erf": "erf",
    "sinh": "sinh",
    "cosh": "cosh",
    "tanh": "tanh",
    "asinh": "asinh",
    "acosh": "acosh",
    "atanh": "atanh",
    "floor": "floor",
    "ceiling": "ceil",
    "re": "real",
    "im": "imag",
    "Min":"std::min",
    "Max":"std::max"
}

# These are the core reserved words in the C language. Taken from:
# http://crasseux.com/books/ctutorial/Reserved-words-in-C.html

reserved_words = ['auto',
                  'if',
                  'break',
                  'int',
                  'case',
                  'long',
                  'char',
                  'register',
                  'continue',
                  'return',
                  'default',
                  'short',
                  'do',
                  'sizeof',
                  'double',
                  'static',
                  'else',
                  'struct',
                  'entry',
                  'switch',
                  'extern',
                  'typedef',
                  'float',
                  'union',
                  'for',
                  'unsigned',
                  'goto',
                  'while',
                  'enum',
                  'void',
                  'const',
                  'signed',
                  'volatile']

class CCodePrinter(CodePrinter):
    """A printer to convert python expressions to strings of c code"""
    printmethod = "_ccode"
    language = "C"

    _default_settings = {
        'order': None,
        'full_prec': 'auto',
        'precision': 15,
        'user_functions': {},
        'human': True,
        'contract': True,
        'dereference': set(),
        'error_on_reserved': False,
        'reserved_word_suffix': '_',
    }

    def __init__(self, settings={}):
        CodePrinter.__init__(self, settings)
        self.known_functions = dict(known_functions)
        userfuncs = settings.get('user_functions', {})
        self.known_functions.update(userfuncs)
        self._dereference = set(settings.get('dereference', []))
        self.reserved_words = set(reserved_words)

    def _rate_index_position(self, p):
        return p*5

    def _get_statement(self, codestring):
        return "%s;" % codestring

    def _get_comment(self, text):
        return "// {0}".format(text)

    def _declare_number_const(self, name, value):
        return "double const {0} = {1};".format(name, value)

    def _format_code(self, lines):
        return self.indent_code(lines)

    def _traverse_matrix_indices(self, mat):
        rows, cols = mat.shape
        return ((i, j) for i in range(rows) for j in range(cols))

    def _get_loop_opening_ending(self, indices):
        open_lines = []
        close_lines = []
        loopstart = "for (int %(var)s=%(start)s; %(var)s<%(end)s; %(var)s++){"
        for i in indices:
            # C arrays start at 0 and end at dimension-1
            open_lines.append(loopstart % {
                'var': self._print(i.label),
                'start': self._print(i.lower),
                'end': self._print(i.upper + 1)})
            close_lines.append("}")
        return open_lines, close_lines

    def _print_Pow(self, expr):
        if "Pow" in self.known_functions:
            return self._print_Function(expr)
        PREC = precedence(expr)
        if expr.exp == -1:
            return '1.0/%s' % (self.parenthesize(expr.base, PREC))
        elif expr.exp == 0.5:
            return 'sqrt(%s)' % self._print(expr.base)
        else:
            return 'pow(%s, %s)' % (self._print(expr.base), self._print(expr.exp))

    def _print_Min(self, expr):
        return 'std::min(%s,%s)' % (self._print(expr.args[0]),self._print(expr.args[1]))

    def _print_Max(self, expr):
        return 'std::max(%s,%s)' % (self._print(expr.args[0]),self._print(expr.args[1]))

    def _print_Equality(self, expr):
        return '(%s==%s)' % (self._print(expr.args[0]),self._print(expr.args[1]))

    def _print_Unequality(self, expr):
        return '(%s!=%s)' % (self._print(expr.args[0]),self._print(expr.args[1]))

    def _print_GreaterThan(self, expr):
        return '(%s>=%s)' % (self._print(expr.args[0]),self._print(expr.args[1]))

    def _print_LessThan(self, expr):
        return '(%s<=%s)' % (self._print(expr.args[0]),self._print(expr.args[1]))

    def _print_StrictLessThan(self, expr):
        return '(%s<%s)' % (self._print(expr.args[0]),self._print(expr.args[1]))

    def _print_StrictGreaterThan(self, expr):
        return '(%s>%s)' % (self._print(expr.args[0]),self._print(expr.args[1]))

    def _print_Sum(self, expr):
        return 'sum<%s>([&](int %s){ return %s; },%s,%s)' % (get_evaluation_ctype_string(expr.function),expr.limits[0][0].name,self._print(expr.function),self._print(expr.limits[0][1]),self._print(expr.limits[0][2]))

    def _print_Rational(self, expr):
        p, q = int(expr.p), int(expr.q)
        return '%d.0/%d.0' % (p, q)

    def _print_Exp1(self, expr):
        return "M_E"

    def _print_Pi(self, expr):
        return 'M_PI'

    def _print_Infinity(self, expr):
        return 'HUGE_VAL'

    def _print_NegativeInfinity(self, expr):
        return '-HUGE_VAL'

#    def _print_Function(self, expr):
#        return "%s(%s)" % (expr.name,',',join([self._print(c) for c in expr.args]))

    def _print_Piecewise(self, expr):
        if expr.args[-1].cond != True:
            # We need the last conditional to be a True, otherwise the resulting
            # function may not return a result.
            raise ValueError("All Piecewise expressions must contain an "
                             "(expr, True) statement to be used as a default "
                             "condition. Without one, the generated "
                             "expression may not evaluate to anything under "
                             "some condition.")
        lines = []
        if expr.has(Assignment):
            for i, (e, c) in enumerate(expr.args):
                if i == 0:
                    lines.append("if (%s) {" % self._print(c))
                elif i == len(expr.args) - 1 and c == True:
                    lines.append("else {")
                else:
                    lines.append("else if (%s) {" % self._print(c))
                code0 = self._print(e)
                lines.append(code0)
                lines.append("}")
            return "\n".join(lines)
        else:
            # The piecewise was used in an expression, need to do inline
            # operators. This has the downside that inline operators will
            # not work for statements that span multiple lines (Matrix or
            # Indexed expressions).
            ecpairs = ["((%s) ? (\n%s\n)\n" % (self._print(c), self._print(e))
                    for e, c in expr.args[:-1]]
            last_line = ": (\n%s\n)" % self._print(expr.args[-1].expr)
            return ": ".join(ecpairs) + last_line + " ".join([")"*len(ecpairs)])
        
    def _print_ArrayFunction(self,expr):
        return 'data_%s(%s)' % ( expr.__class__.__name__ ,  ','.join([self._print(arg) for arg in expr.args]))
    
    def _print_CCodeFunction(self,expr):
        return '%s(%s)' % ( expr.__class__.__name__ ,  ','.join([self._print(arg) for arg in expr.args]))

    def _print_Symbol(self, expr):

        name = super(CCodePrinter, self)._print_Symbol(expr)

        if expr in self._dereference:
            return '(*{0})'.format(name)
        else:
            return name

    def _print_sign(self, func):
        return '((({0}) > 0) - (({0}) < 0))'.format(self._print(func.args[0]))

    def _print_not_supported(self,expr):
        raise ValueError("Cannot convert to c++ code: %s" % str(expr))

    def indent_code(self, code):
        """Accepts a string of code or a list of code lines"""

        if isinstance(code, string_types):
            code_lines = self.indent_code(code.splitlines(True))
            return ''.join(code_lines)

        tab = "   "
        inc_token = ('{', '(', '{\n', '(\n')
        dec_token = ('}', ')')

        code = [ line.lstrip(' \t') for line in code ]

        increase = [ int(any(map(line.endswith, inc_token))) for line in code ]
        decrease = [ int(any(map(line.startswith, dec_token)))
                     for line in code ]

        pretty = []
        level = 0
        for n, line in enumerate(code):
            if line == '' or line == '\n':
                pretty.append(line)
                continue
            level -= decrease[n]
            pretty.append("%s%s" % (tab*level, line))
            level += increase[n]
        return pretty


def ccode(expr,simplify = False, **settings):
    if simplify: expr = expr.doit().simplify()
    return CCodePrinter(settings).doprint(expr)


def create_cfunction(name,args,expr,return_type = None):

    if return_type == None: return_type = get_ctype(expr)
    return_type_name = get_ctype_name(return_type)

    var_definitions = [ get_ctype_name(get_ctype(var)) + ' ' + ('__' if get_ctype(var) == ctypes.c_complex else '') + var.name
                       for var in args ]
    function_header = return_type_name + ' ' + name + '(' + ','.join(var_definitions) + ')'

    function_conversions = '\n'.join([ 'std::complex<double> ' + v.name + ' = create_std_complex(__' + v.name + ');'
                                      for v in args if get_ctype(v) == ctypes.c_complex ])

    return_value = ccode(expr)
    return_code = 'return ' + (return_value if not return_type == ctypes.c_complex else 'create_complex(' + return_value + ')') + ';'
    return 'extern "C" ' + '\n'.join([function_header,'{',function_conversions,return_code,'}'])

def create_c_code(body):
    includes = '\n'.join(['#include <' + name + '>' for name in ('cmath','complex','functional','stdexcept')])

    complex_operators = '''
struct complex{ double real,imag; };
inline std::complex<double> create_std_complex(complex z){ return std::complex<double>(z.real,z.imag); }
inline complex create_complex(std::complex<double> z){ complex c; c.real = z.real(); c.imag = z.imag(); return c; }
inline std::complex<double> operator*(std::complex<double> lhs, int rhs) { return lhs * double(rhs); }
inline std::complex<double> operator/(std::complex<double> lhs, int rhs) { return lhs / double(rhs); }
inline std::complex<double> operator+(std::complex<double> lhs, int rhs) { return lhs + double(rhs); }
inline std::complex<double> operator-(std::complex<double> lhs, int rhs) { return lhs - double(rhs); }
inline std::complex<double> operator*(int lhs, std::complex<double> rhs) { return double(lhs) * rhs; }
inline std::complex<double> operator/(int lhs, std::complex<double> rhs) { return double(lhs) / rhs; }
inline std::complex<double> operator+(int lhs, std::complex<double> rhs) { return double(lhs) + rhs; }
inline std::complex<double> operator-(int lhs, std::complex<double> rhs) { return double(lhs) - rhs; }
template<typename T> inline T abs(const T &v){ return v>=0?v:-v; }

template<typename T,typename F> inline T sum(F f,int begin,int end){
    T res = 0;
    for(int i=begin;i<=end;++i) res += f(i);
    return res;
}
#define I std::complex<double>(0,1)
    '''

    array_wrapper_class = '''
template <class T> struct array_wrapper_1D{
    T * data;
    size_t size;

    array_wrapper_1D(size_t data,size_t size):data((T *)data),size(size){}

    T operator()(size_t i){
        if(i>=size) return 0;
        return data[i];
    }
};

template <class T> struct array_wrapper_2D{
    T * data;
    size_t xsize,ysize;

    array_wrapper_2D(size_t data,size_t xsize,size_t ysize):data((T *)data),xsize(xsize),ysize(ysize){}

    T operator()(size_t x,size_t y){
        if(x>=xsize) return 0;
        if(y>=ysize) return 0;

        int i = x * ysize + y;
        return data[i];
    }
};

template <class T> struct array_wrapper_3D{
    T * data;
    size_t xsize,ysize,zsize;

    array_wrapper_3D(size_t data,size_t _xsize,size_t _ysize,size_t _zsize):data((T *)data),xsize(_xsize),ysize(_ysize),zsize(_zsize){}

    T operator()(size_t x,size_t y,size_t z){
        if(x>=xsize) return 0;
        if(y>=ysize) return 0;
        if(z>=zsize) return 0;

        int i = x * (ysize*zsize) + y * zsize + z;
        return data[i];
    }
};
    '''

    return '\n\n'.join([includes,complex_operators,array_wrapper_class] + body)

def compile_code(code,output_directory,name):
    object_file = output_directory+'/'+name+'.o'
    p = Popen(['g++','-o',object_file,'-c','-xc++','-std=c++11','-O3','-Ofast','-march=native','-fPIC', '-'],stdin=PIPE, stdout=PIPE, stderr=PIPE)
    p.stdin.write(code)
    p.stdin.close()
    return_code = p.wait()
    if(return_code!=0):
        raise RuntimeError("Cannot compile expression: " + p.stderr.read() + "\nFull code:\n\n" + code)

    shared_library = output_directory+'/'+name+'.so'
    p = Popen(['g++','-shared','-o',shared_library,object_file],stdin=PIPE, stdout=PIPE, stderr=PIPE)
    p.stdin.write(code)
    p.stdin.close()
    return_code = p.wait()
    if(return_code!=0):
        raise RuntimeError("Cannot convert to shared library: " + p.stderr.read())

    return shared_library

count = 0

def get_global_variables(expr,found = None):
    if found == None: found = set()
    if isinstance(expr,(CCodeFunction,ArrayFunction)): 
        found.add(expr)
    for arg in expr.args:
         get_global_variables(arg,found)
    return found

#cache = {}

def compile_sympy_expressions(name_expr_vars_ret,globals = []):
    """
    Compiles a set of sympy expressions into a c++ library callable from python.
    
    Parameters
    ----------
    name_expr_vars_ret: list containing 4-tuples: (name,expression,arguments,return_type)
                        Each tuple represents a function that will be accessable from python.
                        name:           string
                                        name of the function
                        expression:     sympy expression
                                        the return value of the function
                        arguments:      tuple of sympy symbols
                                        a tuple containing the arguments of the function
                        return_type:    ctypes type [optional, defaults to c_complex if not derivable from the expression]
                                        the return type of the expression
    globals:   list [optional]
               Global variables which can be modified after compilation
    
    Returns
    -------
    A library loaded through ctypes containing callable wrappers of the functions as members.
    
    Notes
    -----
    Each name must be a unique c++ identifier.
    The functions have an additional attribute `address` which is the memory address of the function
    
    Example
    -------
    >>> x,y = sympy.symbols("x,y",real=True)
    >>> f = x-1
    >>> g = x**2+y
    >>> lib = compile_sympy_expressions([("f",f,(x,)),("g",g,(x,y))])
    >>> lib.f(1)
    0.0
    >>> lib.g(2,1)
    5.0
    
    See also
    --------
    compile_sympy_expression:   compile a single expressions
    """
    
    #key = str(name_expr_vars_ret)
    
    #try:
    #    return cache[key]
    #except:
    #    pass
    
    code_fragments = []
    
    if not isinstance(globals,(list,tuple)):
        globals = [globals]
    
    for g in globals:
        code_fragments.append("%s %s;" % (get_evaluation_ctype_string(g),str(g)))

    global_variables = set()

    for evr in name_expr_vars_ret:
        if len(evr) == 4: name,expr,args,return_type = evr
        else:
            name,expr,args = evr
            return_type = get_ctype(expr)

        global_variables.update(get_global_variables(expr))
        code_fragments.append(create_cfunction(name,args,expr,return_type))

    code_fragments.insert(0,'\n'.join([g.get_ccode() for g in global_variables]))

    code = create_c_code(code_fragments)

    compile_dir = tempfile.mkdtemp()
    
    try:
        shared_library = compile_code(code,compile_dir,'compiled_sympy_lib')
    except:
        import shutil
        shutil.rmtree(compile_dir)
        raise

    lib = ctypes.cdll.LoadLibrary(shared_library)
    import shutil
    shutil.rmtree(compile_dir)

    class CompiledLib:
        pass

    CL = CompiledLib()

    CL.globals = { g.name:get_ctype(g).in_dll(lib,g.name) for g in globals }

    for evr in name_expr_vars_ret:
        if len(evr) == 4: name,expr,args,return_type = evr
        else:
            name,expr,args = evr
            return_type = get_ctype(expr)

        f = getattr(lib,name)
        f.argtypes = [get_ctype(arg) for arg in args]
        f.restype  = return_type
        f.address = ctypes.cast(f, ctypes.c_void_p).value
        setattr(CL,name,f)

    CL.code = code
    CL.arrays = global_variables
    
    #cache[key] = CL
    
    return CL

def compile_sympy_expression(expr,arguments,return_type = None, globals = []):
    """
    Compiles a sympy expressions into a c++ callable function.
    
    Parameters
    ----------
    expression:     sympy expression
                    the return value of the function
    arguments:      tuple of sympy symbols
                    a tuple containing the arguments of the function
    return_type:    ctypes type [optional, defaults to c_complex if not derivable from the expression]
                    the return type of the expression
    globals:        list [optional]
                    A list containing objects generating code that should be placed in the global namespace before the functions. Each element should implement the `get_ccode` member returning a string of c++ code.
    
    Returns
    -------
    A wrapped c++ function callable from python
    
    Example
    -------
    >>> x,y = sympy.symbols("x,y",real=True)
    >>> f = compile_sympy_expression(x**2-y,(x,y))
    >>> f(2,1)
    3.0
    
    See also
    --------
    compile_sympy_expressions:  simultaneously compile multiple expressions
    """
    
    if return_type: functions = [("f",expr,arguments,return_type,)]
    else: functions = [("f",expr,arguments)]

    lib = compile_sympy_expressions(functions,globals )
    lib.f.globals = lib.globals
    lib.f.arrays = lib.arrays

    return lib.f

