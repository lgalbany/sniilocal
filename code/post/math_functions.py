def gaussian(x, vec = []):  # IDL version
    '''
    IDL-like input: 
    x = data array
    vec[0] = amplitude (maximum value at x = mean)
    vec[1] = mean
    (Optional) 
    vec[2] = sigma
    vec[3] = constant offset
    
    The output is a number. To evaluate a symbolic value, do
    import sympy as sp
        
    ga = gaussian(Symbol('x'),[1,2])
    print ga
    '''
        
    import numpy as np
        
    A = vec[0]
    mu = vec[1]
    
    
    if len(vec) == 3:
        
        sigma = vec[2]
        
    else:
        
        sigma = 1./(np.sqrt(2.*np.pi)*vec[0])
        
        
    if len(vec) == 4:
        
        return A*np.exp(-np.power((x - mu)/sigma, 2.)/2.) + vec[3]
    
    else:
        
        return A*np.exp(-np.power((x - mu)/sigma, 2.)/2.)
    
    
    
    
def polynomial(x = [], c = []): # Python has np.poly1d, but I want the opposite, as in IDL, C_0 + C_1*x + C_2*x**2 +...
    '''
    Return c0 + c1*x +c2*x^2 + ...
    Arguments:
        X. The variable. This value can be a scalar, vector or array.
        C. The vector of polynomial coefficients. The degree of the polynomial is len(c) -1.
        
    The output is a number value. To evaluate a symbolic value, do
    import sympy as sp
    
    pl = poly(sp.Symbol('x'), [1,2,3])
    print pl
    '''
        
    import numpy as np
        
    # Beware of numpy.float32/64
    if (type(c) is float or type(c) is np.float32 or type(c) is np.float64 or type(c) is int) : c = [c]
    
    pl = 0.
    
    for cidx, cel in enumerate(c):
        
        pl += cel * x**(cidx) 
        
    return pl