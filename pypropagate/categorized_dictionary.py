

class Category(object):

    """Category class for CategorizedDictionary."""
    
    def __init__(self, parent):
        self._parent = parent
        self._keys = {}
        self.subcategories = {}
        self.attributes = {}
        self.methods = {}
        self._locked = False
        self.__doc__ = ""
        self._key_doc = {}
        self._anon_count = 0
        self.__initialised = True

    def copy(self,parent=None,copy=None):
        "Return a copy of the Category."
        if parent == None:
            parent = self._parent
        if copy == None:
            copy = Category(parent)

        copy._keys = self._keys.copy()
        copy._locked = self._locked
        copy._key_doc = self._key_doc.copy()
        copy.__doc__ = self.__doc__
        
        for name,key in copy._keys.iteritems():
            copy.__dict__[name] = key
        
        for name,cat in self.subcategories.iteritems():
            catcop = cat.copy(parent = copy)
            copy.subcategories[name] = catcop
            copy.__dict__[name] = catcop

        for name,method in self.methods.iteritems():
            copy.add_method(name,method)

        for name,attribute in self.attributes.iteritems():
            copy.add_attribute(name,attribute)

        return copy
        
    def lock(self,name = None,reason = None):
        """If no name parameter provided this prevents new keys and subcategories to be added. Otherwise associated the key is locked."""
        if name == None:
            self._locked = True
        else:
            self._lock(self.get_key(name),reason)
    
    def unlock(self,name = None):
        "reverses the Category.locked"
        if name == None:
            self._locked = False
        else:
            self._unlock(self.get_key(name))
    
    def _get_top(self):
        if self._parent != None:
            return self._parent._get_top()
        return self
    
    def __str__(self):
        res = []
        if len(self.__doc__)>0:
            res += [self.__doc__]
        if len(self._keys) > 0:

            def shorten_value(v):
                s = str(v)
                if len(s) < 13: return s
                return '%s...' % s[0:10]

            res += ["Keys:\n-----"]
            res += ["%s:\t%s" % 
                    (("%s = %s" % (name,shorten_value(self.get_value(name)))) if self.get_value(name)!=None else name,
                     self._key_doc[name] if name in self._key_doc 
                                         else self._keys[name]) for name in self._keys]
        if len(self.subcategories) > 0: res += ["Subcategories:\n--------------"]
        res += ["%s:\t%s" % (name,self.subcategories[name].__doc__) for name in self.subcategories.keys()]
        return '\n'.join(res)
    
    def set_info(self,name,info):
        """Set info for a key name"""
        if info != None: 
            self._key_doc[name] = info
        elif name in self._key_doc:
            del self._key_doc[name]
    
    def info(self,name = None):
        """Get info for the category or for a key name"""
        if name==None: return str(self)
        return self._key_doc[name] if name in self._key_doc else None
    
    def add_category(self,name,category):
        """Adds a subcategory."""
        if self._locked:
            raise AttributeError("modifying locked category")
        if name in self.subcategories:
            import warnings
            warnings.warn("overwriting subcategory %s" % name, UserWarning)
            del self.subcategories[name]
            if name in self.__dict__:
                del self.__dict__[name]
        self.subcategories[name] = category
        self._set_attribute(name,category)
    
    def create_category(self,name,info=None):
        """Creates a named subcategory."""
        new_category = Category(self)
        if info != None:
            new_category.__doc__ = info
        self.add_category(name,new_category)
        return new_category
    
    def get_category(self,name):
        """Gets the named subcategory."""
        return self.subcategories[name]
    
    def has_category(self,name):
        """Check if the category has a subcategroy with the name."""
        return name in self.subcategories
    
    def has_name(self,name):
        return name in self._keys
    
    def get_key(self,name):
        try:
            return self._keys[name]
        except KeyError:
            raise NameError('name %s not in category' % name)
    
    def get_value(self,name):
        return self._get_value(self.get_key(name))
    
    def _get_value(self,key):
        return self._parent._get_value(key)
    
    def _key_exisits(self,key):
        return self._parent._key_exisits(key)
    
    def create_key(self,name,key,value = None,info=None):
        """Creates a new key in the categorized dictionary."""
        if self._locked:
            raise AttributeError("modifying locked category")
        if name in self.__dict__: 
            import warnings
            warnings.warn("overwriting attribute %s" % name, UserWarning)
        if self._key_exisits(key):
            import warnings
            warnings.warn("overwriting key %s" % key, UserWarning)
            self._unlock(key)
        
        self._create_key(name,key,value,info,self)
        self.add_key(name,key,info)

        return key
    
    def add_key(self,name,key,info=None,warn=True):
        """Adds an existing key to the category"""
        if self._locked:
            raise AttributeError("modifying locked category")
        if not self._key_exisits(key):
            raise KeyError("adding undefined key %s" % key)
        if warn and name in self.__dict__: 
            import warnings
            warnings.warn("overwriting attribute %s" % name, UserWarning)
        if name == None:
            self._anon_count += 1
            name = "_anonymous_key_%s" % self._anon_count
        self._keys[name] = key
        if info!=None:
            self._key_doc[name] = info
        self._set_attribute(name,key)

        return key

    def add_attribute(self, attr, value):
        if attr in self.__dict__:
            import warnings
            warnings.warn("overwriting attribute %s" % name, UserWarning)
        self.attributes[attr] = value
        self.__dict__[attr] = value

    def add_method(self, attr, value):
        if attr in self.__dict__:
            import warnings
            warnings.warn("overwriting attribute %s" % name, UserWarning)
        self.methods[attr] = value

        from types import MethodType
        self.__dict__[attr] = MethodType(value,self)

    def remove_key(self,key):
        """Removes a key from the category. The global key still remains valid."""
        for name, comp in self._keys.iteritems():
            if key == comp:
                self.remove_name(name) 
                return
        raise ValueError("attempting to remove undefined key %s" % key)

    def remove_name(self,name):
        """Removes a name from the category. The associated key still remains valid."""
        
        if name not in self._keys:
            import warnings
            warnings.warn("attempting to remove undefined name %s" % name, UserWarning)
            return

        del self._keys[name]
        del self.__dict__[name]
        if name in self._key_doc:
            del self._key_doc[name]
    
    def _create_key(self,name,key,value,info,sender):
        self._parent._create_key(name,key,value,info,sender)
    
    def _set_value(self,key,value):
        "Sets the value of a key in the dictionary"
        self._parent._set_value(key,value)
        
    def _get_value(self,key):
        "Gets the value of a key in the dictionary"
        return self._parent._get_value(key)
            
    def _lock(self,key,reason):
        self._parent._lock(key,reason)
   
    def _unlock(self,key):
        self._parent._unlock(key)
        
    def _set_attribute(self,attr,value):
        super(Category,self).__setattr__(attr, value)
        
    def __setattr__(self, item, value):
        if not self.__dict__.has_key('_Category__initialised'):
            return self._set_attribute(item, value)
        is_key = self._keys.has_key(item)
        if not is_key:
            if hasattr(self,item):
                return super(Category,self).__setattr__(item, value)
            raise AttributeError('Category does not contain name %s' % item)
        else:
            self._set_value(self._keys[item],value)
                            
    def keys(self):
        return self._keys.values()
    
    def names(self):
        return self._keys.keys()
    
    def dictionary(self,keys=None):
        """Return a dictionary with the defined keys of this category"""
        if keys == None:
            keys = self.keys()
        return self._parent.dictionary(keys)
    
    def export(self,dictionary,**kwargs):
        """Copies name/key pairs into a target dictionary. Can be used to load variables into the local namespace by passing the globals() dictionary. Note that exisiting variables will be overwritten."""
        if isinstance(dictionary,Category):
            for name in self.names():
                info = self.info(name)
                dictionary.add_key(name,self.get_key(name),info=info,**kwargs)
        else:
            dictionary.update(self._keys,**kwargs)


class CategorizedDictionary(Category):
    
    """
    Auto-complete compatible dictionary wrapper where keys can be grouped into categories. Keys can be documented using the named argument info during creation. Useful for storing settings.
    
    Members / Initializer Parameters:
    --------------------------------
    data: Dictionary containing key/value pairs
    
    Example:
    --------
    >>> settings = CategorizedDictionary()
    >>> settings.create_category("variables")
    >>> settings.create_category("paramters")
    >>> settings.paramters.create_key("i","index")
    >>> settings.variables.create_key("x","Var(x)",1)
    >>> settings.variables.create_key("y","Var(y)")
    >>> settings.variables.y = 42
    >>> settings.dictionary()
    {'Var(x)': 1, 'Var(y)': 42}
    >>> settings.undefined_keys()
    {'index'}
    """
    
    def __init__(self,data=None):
        if data == None:
            data = {}
        self.data = data
        self.locked_keys = {}
        super(CategorizedDictionary,self).__init__(self)
    
    def copy(self,copy = None):
        "Return a copy of the CategorizedDictionary."
        if copy == None:
            copy = CategorizedDictionary()
        copy.data = self.data.copy()
        copy.locked_keys = self.locked_keys.copy()
        super(CategorizedDictionary,self).copy(copy=copy)
        return copy
    
    def _lock(self,key,reason):
        self.locked_keys[key] = reason
    
    def _unlock(self,key):
        try:
            del self.locked_keys[key]
        except KeyError:
            pass
    
    def _set_value(self,key,value):
        if key in self.locked_keys:
            raise ValueError('key %s is locked. reason: %s' % (key,self.locked_keys[key]))
        self.data[key] = value 
        
    def _get_value(self,key):
        try:
            return self.data[key]
        except KeyError:
            raise KeyError("Value for key %s undefined" % key)
    
    def _key_exisits(self,key):
        return key in self.data
    
    def _create_key(self,name,key,value,info,sender):
        if self._key_exisits(key) and value == None:
            return
        self._set_value(key,value)
        
    def all_keys(self):
        return self.data.keys()
    
    def __setitem__(self, key, value):
        self._set_value(key,value)
    
    def __getitem__(self,key):
        return self.data[key]
        
    def undefined_keys(self):
        """Get all keys with value equal to None"""
        res = set()
        for key in self.data:
            if self.data[key] == None:
                res.add(key)
        return res

    def is_defined(self,key):
        return self.data[key] != None

    def defined_keys(self):
        """Get all keys with value not equal to None"""
        return set(self.get_all_keys()) - self.get_undefined_keys()
    
    def __repr__(self):
        return "<CategorizedDictionary: %s,Subcategories: %s>" % (self._keys.keys(),self.subcategories.keys())
    
    def dictionary(self,keys=None):
        """Return a dictionary with all defined keys and values of this category"""
        if keys == None:
            keys = self.data.keys()
        return {key:self.data[key] for key in keys if self.data[key]!=None}
