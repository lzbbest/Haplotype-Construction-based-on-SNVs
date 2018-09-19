import numpy as np
import midaco

def ACOopt(ratio,p):
    
    # mathematical formulation of the general multi-objective MINLP
    
    def problem_function(x):
        # x[0] = M             x[1] = m
        f = [0.0]*1 # Initialize array for objectives F(X)
        g = [0.0]*3 # Initialize array for constraints G(X)
    
        # Objective functions F(X)
    
        y = lambda p,M,m: (p*M+1-p)/(p*m+1-p)
        #f[0] = abs(y(p,x[0],x[1])-ratio[4]) + abs(y(p,x[0],x[1])-ratio[5])
        for i in range(len(ratio)): f[0] = f[0] + abs(y(p,x[0],x[1])-ratio[i])
        
        #  Equality constraints G(X) = 0 MUST COME FIRST in g[0:me-1]
        #g[0] = x[0] - 1.0
        # Inequality constraints G(X) >= 0 MUST COME SECOND in g[me:m-1] 
        g[0] = x[0] -x[1]       
        g[1] = x[0]
        g[2] = x[1]
        
        return f,g
    
    
    #  MAIN PROGRAM  
    
    key = b'MIDACO_LIMITED_VERSION___[CREATIVE_COMMONS_BY-NC-ND_LICENSE]'
    
    problem = {} # Initialize dictionary containing problem specifications
    option  = {} # Initialize dictionary containing MIDACO options
    
    problem['@'] = problem_function # Handle for problem function name
    
    # Step 1: Problem definition
    
    # STEP 1.A: Problem dimensions
    ##############################
    problem['o']  = 1  # Number of objectives 
    problem['n']  = 2  # Number of variables (in total) 
    problem['ni'] = 2  # Number of integer variables (0 <= ni <= n) 
    problem['m']  = 3  # Number of constraints (in total) 
    problem['me'] = 0  # Number of equality constraints (0 <= me <= m) 
    
    # STEP 1.B: Lower and upper bounds 'xl' & 'xu'  
    ##############################################  
    problem['xl'] = [ 0,0 ]
    problem['xu'] = [ 100,2 ]
    
    # STEP 1.C: Starting point 'x'  
    ##############################  Initialize
    problem['x'] = problem['xl'] # Here for example: starting point = lower bounds
        
    ########################################################################
    ### Step 2: Choose stopping criteria and printing options    ###########
    ########################################################################
       
    # STEP 2.A: Stopping criteria 
    #############################
    option['maxeval'] = 10000     # Maximum number of function evaluation (e.g. 1000000) 
    option['maxtime'] = 60*60*24  # Maximum time limit in Seconds (e.g. 1 Day = 60*60*24) 
    
    # STEP 2.B: Printing options  
    ############################ 
    option['printeval'] = 1000  # Print-Frequency for current best solution (e.g. 1000) 
    option['save2file'] = 1     # Save SCREEN and SOLUTION to TXT-files [0=NO/1=YES]
    
    ########################################################################
    ### Step 3: Choose MIDACO parameters (FOR ADVANCED USERS)    ###########
    ########################################################################
    
    option['param1']  = 0.0  # ACCURACY  
    option['param2']  = 0.0  # SEED  
    option['param3']  = 0.0  # FSTOP  
    option['param4']  = 0.0  # ALGOSTOP  
    option['param5']  = 0.0  # EVALSTOP  
    option['param6']  = 0.0  # FOCUS  
    option['param7']  = 0.0  # ANTS  
    option['param8']  = 0.0  # KERNEL  
    option['param9']  = 0.0  # ORACLE  
    option['param10'] = 0.0  # PARETOMAX
    option['param11'] = 0.0  # EPSILON  
    option['param12'] = 0.0  # BALANCE
    option['param13'] = 0.0  # CHARACTER
    
    ########################################################################
    ### Step 4: Choose Parallelization Factor   ############################
    ########################################################################
    
    option['parallel'] = 0 # Serial: 0 or 1, Parallel: 2,3,4,5,6,7,8...
    
    ########################################################################
    ############################ Run MIDACO ################################
    ########################################################################
    solution = midaco.run( problem, option, key )
    var_list = solution['x']
    return var_list

'''
if __name__ == '__main__': 
    ratio = np.array([1.1,1.4,1.6,1.8,2.1,1.6,1.7])
    p = 0.65
    var_list = ACOopt(ratio,p)
'''

         
