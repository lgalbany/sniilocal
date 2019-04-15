def multipeak(p, x, y, err, fjac=None, fit=False, npeak=0):
    
    '''Parameter values are passed in "p"
       If fjac==None then partial derivatives should not be computed.
       It will always be None if MPFIT is called with default flag.
       Non-negative status value means MPFIT should continue, negative means
       stop the calculation.'''
    
    # See this for exception handling: https://doughellmann.com/blog/2009/06/19/python-exception-handling-techniques/
        
    import numpy as np
    import matplotlib.pyplot as plt
    from math_functions import gaussian, polynomial
    

    model = np.zeros(len(y)) # Double precision in the original IDL code by default in Python (dtype = float64)
    n = len(x)
    
    
    for i in range(npeak):
        
        ins = np.where(abs(x - p[3*i+1]) <= 5*p[3*i+2])[0]   # Points within 5 sigmas (x-p are the residuals)
        cn = len(ins) # Equivalent no np.count_nonzero(ins)
        
        if cn > 0:
            
            ins1 = np.where((ins >= 0) & (ins < n))[0]
            cn1 = len(ins1)
            
            if cn1 > 0:
                
                # CAUTION!!!!!!!!!!!: ORIGINAL IDL READS
                #model[ins[ins1]] = model[ins[ins1]] + gaussian(x[ins[ins1]], p[3*i:3*i+2])
                # BUT PYTHON HAS A DIFFERENT COLON INDEXING, AND THE UPPER LINE
                # WOULD NOT TAKE SIGMA INTO ACCOUNT
                model[ins[ins1]] = model[ins[ins1]] + gaussian(x[ins[ins1]], p[3*i:3*i+3])
                
                
                # plt.plot(x[ins[ins1]], model[ins[ins1]])
        # plt.plot(x, y, color='black')
        
    # Add a line with the last two parameters from parinfo
    # Intercept: np.min(y), slope: 0 (assuming a flat continuum)
    if len(p) - 3*npeak > 0: 
        model += polynomial(x, p[3*npeak:])

        
    res = (y - model)/err

    if fit:  # User has set 'fit'
        res = model
    
        
    status = 0
    
    # Strictly necessary in the Python version of MPFIT. 
    # See the MPFIT example in the documentation
    return [status, res] 









def multipeak_fit(x, y, initpars, npeak, parinfo, fit=False, \
                  err=None, maxiter=200, ftol=1.E-10, gtol=1.E-10, xtol=1.E-10, plot=False, silent=False): 
    
    import numpy as np
    import matplotlib.pyplot as plt
    from mpfit import mpfit
    from math_functions import gaussian, polynomial
    
    
    if err is None:
        err = [1.0]*len(x)
        
    quiet = 0
    
    if silent: # quiet = 1 will print no output 
        quiet = 1
    
    
    p0 = initpars #p.params
    
    fa = {'x':x, 'y':y, 'err':err, 'npeak':npeak}
    # xtol and ftol are set to 1.E-10 by default in MPFIT
    p = mpfit(multipeak, p0, functkw = fa, parinfo = parinfo, maxiter=maxiter, ftol=ftol, gtol=gtol, xtol=xtol, quiet=quiet) 
    
    #print 'p.niter = ',p.niter # Number of iterations
    #print 'p.status = ',p.status # Status. Look at ftol and xtol
    #print 'p.params', p.params # Fitted values for the initial parameters
    #print 'p.error = ',p.perror # Covariance matrix with formal 1-sigma errors
    #print 'p.fnorm = ',p.fnorm # Chi2
    #print 'DOF = ',len(x) - len(p.params) # Degrees of freedom
    #print 'p.covar = ', p.covar # Covariance matrix (11x11 in this case)
    #print 'p = ', p
    
    ft = multipeak(p.params, x = x, y = y, err = err, npeak = npeak, fit=fit)[1] 
    # Where does the [1] come from? See last line of multipeak
    
    
    if plot: # Plot data and final model
        
        plt.figure()
        plt.plot(x, y, color='black', ls='dashed', lw=1.5)
        plt.plot(x, ft, color='purple', ls='solid', lw=1.5)

        # Considers a "fourth" peak for Halpha + NII_a + NII_b + Halpha_SNR
        color_array = np.array(['blue', 'red', 'darkorange', 'turquoise']) 
        for np in range(npeak):
            plt.plot(x, gaussian(x, p.params[3*np:3*np+3]) + polynomial(x, p.params[3*npeak:]),\
                color=color_array[np])


        
    return [p, ft]