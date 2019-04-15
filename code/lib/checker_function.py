#First, you can use special formal parameter syntax *. 
#If the function definition has a formal parameter 
#preceded by a single *, then Python populates that 
#parameter with any positional parameters that are not 
#matched by preceding formal parameters (as a tuple). 
#If the function definition has a formal parameter 
#preceded by **, then Python populates that parameter 
#with any keyword parameters that aren't matched by 
#preceding formal parameters (as a dict). The function's
#implementation can check the contents of these 
#parameters for any "optional parameters" of the sort 
#you want.

#For instance, here's a function opt_fun which takes two 
#positional parameters x1 and x2, and looks for another 
#keyword parameter named "optional".

#>>> def opt_fun(x1, x2, *positional_parameters, **keyword_parameters):
#...     if ('optional' in keyword_parameters):
#...         print 'optional parameter found, it is ', keyword_parameters['optional']
#...     else:
#...         print 'no optional parameter, sorry'
#... 
#>>> opt_fun(1, 2)
#no optional parameter, sorry
#>>> opt_fun(1,2, optional="yes")
#optional parameter found, it is  yes
#>>> opt_fun(1,2, another="yes")
#no optional parameter, sorry

#Second, you can supply a default parameter value of some 
#value like None which a caller would never use. If the 
#parameter has this default value, you know the caller did 
#not specify the parameter. If the parameter has a non-default 
#value, you know it came from the caller.



class CheckerFunction(object):
    def __init__(self, function, **defaults):
        self.function = function
        self.defaults = defaults

    def __call__(self, **kwargs):
        for key in self.defaults:
            if(key in kwargs):
                if(kwargs[key] == self.defaults[key]):
                    print 'passed default'
                else:
                    print 'passed different'
            else:
                print 'not passed'
                kwargs[key] = self.defaults[key]

        return self.function(**kwargs)

    
    
#def f(a):
#    print a

#check_f = CheckerFunction(f, a='z') 
#check_f(a='z') # Passed default
#check_f(a='b') # Passed different
#check_f() # Not passed