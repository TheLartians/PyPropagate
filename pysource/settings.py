from .categorized_dictionary import Category,CategorizedDictionary


class Settings(CategorizedDictionary):
    """
    Stores the parameters for a simulation. This class inherits from CategorizedDictionary. Initially, it contains four categories named finite_differences, simulation_box, numerics and base_units. See the info(name) member function for more information about the categories and keys.
    There is an additional symbols member containing aliases to all keys. Adding a key with the same name will overwrite the name in the symbols member.
    The `get` function provides a method of transforming an expression using all known substitutions defined by the key/value pairs of this dictionary.
    """

    def __init__(self,create_categories = True):
        import collections

        self.updaters = collections.OrderedDict()
        self._cache = dict()
        self._eval_cache = dict()
        self._initialized = True
        self._initializers = collections.OrderedDict()

        super(Settings,self).__init__()
        
        self.__doc__ = "Settings of a FiniteDifference simulation"
        
        if create_categories :
            self.create_category("numerics",info="pure numerical values")
            self.create_category("unitless",info="conversions to unitless coordinates")
            self.create_category("symbols",info="symbol aliases to frequently used symbols")

    def create_category(self,cat_name,*args,**kwargs):
        import pycas as pc
        cat = super(Settings, self).create_category(cat_name,*args,**kwargs)

        def create_symbol(name,value=None,info=None,**kwargs):
            return cat.create_key(name,pc.Symbol("%s_%s" % (name,cat_name),**kwargs),value,info)
        def create_function(name,args,value=None,info=None,**kwargs):
            return cat.create_key(name,pc.Function("%s_%s" % (name,cat_name),**kwargs)(*args),value,info)

        cat._set_attribute('create_symbol',create_symbol)
        cat._set_attribute('create_function',create_function)

        return cat

    @property
    def initializers(self):

        parent = self

        class Initializers(object):
            def __getitem__(self, item):
                return parent._initializers[item]
            def __setitem__(self, key, value):
                parent._initialized = False
                parent._initializers.__setitem__(key,value)

        return Initializers()

    def initialize(self):
        if self._initialized:
            return
        self._initialized = True

        for initializer in self._initializers.values():
            initializer(self)


    def copy(self):
        "Return a soft copy of the Settings object."
        copy = Settings(create_categories = False)
        super(Settings,self).copy(copy = copy)
        copy._initializers = (self._initializers.copy())
        copy.updaters = (self.updaters.copy())
        return copy

    def _is_numeric(self,value):
        "Tests if a value should be classified as numeric"
        if isinstance(value,(bool,int,long,float,complex)):
            return True

        import pycas as pc
        if isinstance(value,pc.Expression):
            try:
                value.N()
                return True
            except:
                pass

            import units
            if units.contains_unit(value):
                return True

        return False

    def _set_value(self,key,value):
        import pycas

        if isinstance(value,pycas.Expression):
            value = value.evaluate(cache = self.get_cache())

        super(Settings,self)._set_value(key,value)
        is_numeric = self._is_numeric(value)
        
        if key not in self.numerics.keys() and is_numeric:
            self.numerics.add_key(None,key)
        elif key in self.numerics.keys() and not is_numeric:
            self.numerics.remove_key(key)

        self._eval_cache = {}
        self._initialized = False

    def _get_evaluator(self,numeric = False,unitless = False):

        key = (numeric,unitless)
        try:
            return self._eval_cache[key]
        except:
            pass


        from pycas import Expression,RewriteEvaluator,ReplaceEvaluator,MultiEvaluator,Wildcard,S,Tuple

        replacement_evaluator = ReplaceEvaluator(recursive=True)
        rule_evaluator = RewriteEvaluator(recursive=True)
        evaluator = MultiEvaluator(recursive=True)

        evaluator.add_evaluator(replacement_evaluator)
        evaluator.add_evaluator(rule_evaluator)

        numeric_keys = self.numerics.keys()
        unitless_keys = self.unitless.keys()

        for s,r in self.data.iteritems():

            if (not numeric and (s in numeric_keys)) or (not unitless and (s in unitless_keys)):
                continue

            if isinstance(s,Expression) and isinstance(r,(Expression,int,float,complex)):
                sr = S(r)
                if s==sr:
                    continue
                if s.is_function:
                    wc_args = {arg:Wildcard(arg.name) for arg in s.args}
                    rule_evaluator.add_rule(s.subs(wc_args),sr.subs(wc_args))
                else:
                    replacement_evaluator.add_replacement(s,r)

        self._eval_cache[key] = evaluator

        return evaluator

    def clear_cache(self,*args):
        self._cache = {}

    def get_cache(self, *args):
        return self._cache

    def get(self,expr,numeric = False,unitless = False,evaluate = True):
        self.initialize()
        evaluator = self._get_evaluator(numeric, unitless)

        res = evaluator(expr,cache = self._eval_cache)
        if evaluate == True:
            res = res.evaluate(cache = self.get_cache())
        return res

    def get_numeric(self,expr,**kwargs):
        return self.get(expr,numeric=True,**kwargs)

    def get_unitless(self,expr,**kwargs):
        return self.get(expr,numeric=True,unitless=True,**kwargs)

    def get_optimized(self,expr,**kwargs):
        from pycas.evaluators.optimizers import optimize_for_compilation
        return optimize_for_compilation(self.get_unitless(expr,**kwargs),cache = self.get_cache())


    def get_definition(self,expr):
        return self.data[expr]

    def get_as(self,expr,type):
        import pycas as pc
        res = self.get_unitless(expr)
        if res.function == pc.Tuple:
            return [type(e) for e in res]
        return type(res)

    def get_numeric_transform(self):
        return lambda x: self.get_numeric(x)









