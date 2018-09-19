#import numpy as np
import pulp
import logging
#import pandas as pd
logger = logging.getLogger("AppName")
class MILPsolver(object):
    def __init__(self,dp4,c_i):
        self.dp4 = dp4
        self.c_i = c_i

    def optimal_one_seg(self,ratio,ci):  # optimize in all segs
        # pho=tumor_percent  in k-th segments
        # ratio: each locus M/m ratio type:np.array
        n = len(ratio)
        model = pulp.LpProblem("Major-minor ploidy problem", pulp.LpMinimize)
        # variables
        M = pulp.LpVariable('M', lowBound=0,upBound = 10,cat='Integer')
        m = pulp.LpVariable('m', lowBound=0,upBound = 3,cat='Integer')
        lalpha = [pulp.LpVariable('e' + str(i), lowBound=0) for i in range(n)] 
        cc = pulp.LpVariable('cc',lowBound=0)
        # Objective
        result = 0.0
        result += cc*n
        # dp4
        for i in range(n):
            result += lalpha[i]
        model += result, 'obj'
        # Constraints
        model += M + m -ci >= -cc
        model += M + m -ci <= cc
        for i in range(n): #Ci
            model += M - m*ratio[i]  >= -lalpha[i]
            model += M - m*ratio[i]  <= lalpha[i]
        model.solve()
        '''
        try:
            model.solve()
        except Exception:
            logger.debug('Problem infeasible')
        '''
        pulp.LpStatus[model.status]
        M = pulp.value(M)
        m = pulp.value(m)
        return M,m
    
    # optimize all segs
    def optimal(self):
        MM,mm = [],[]
        nsegs = len(self.c_i)
        for i in range(nsegs):
            ci = self.c_i[i]
            ratio = self.dp4[i]
            M,m = self.optimal_one_seg(ratio,ci)
            MM.append(M)
            mm.append(m)
        return MM,mm
        

'''
x1 = pulp.LpVariable("x1", 0, 500)
x2 = pulp.LpVariable("x2", 0, 6)
x3 = pulp.LpVariable("x3", 0, 10)
x4 = pulp.LpVariable("x4", 0, 8)
prob = pulp.LpProblem("myProblem", pulp.LpMinimize)
prob += 0.5*x1 + 0.2*x2 + 0.3*x3 + 0.8*x4
prob += 400*x1 + 200*x2 + 150*x3 + 500*x4 >= 500
prob += 3*x1 + 2*x2 + 0 + 0 >= 6
prob += 2*x1 + 2*x2 + 4*x3 + 4*x4 >= 10
prob += 2*x1 + 4*x2 + 1*x3 + 5*x4 >= 8
status = prob.solve()
pulp.LpStatus[status]
print(pulp.value(x1))
pulp.value(x2)
pulp.value(x3)
pulp.value(x4)

def w(dp4file):
    dp4,i = [],0
    with open(dp4file) as r:
        for line in r:
            info = line.split()
            info = list(map(float,info))
            temp = np.array(info)
            dp4.append(temp)
            i+=1
            #print(i)
    return dp4


if __name__ == '__main__': 
    #dp4 = [np.array([2.1,1.6,1.7]),np.array([3.9,3.6,3.7]),np.array([2.1,2.6])]
    #pho = np.array([0.99,0.99,0.99])
    #c_i = np.array([10,10,10])
    dp4 = w('dp4.csv')
    #dp4 = dp4[1:3]
    c_i = np.loadtxt('c_i.txt')
    #c_i = c_i[1:3]
    solver = MILPsolver(dp4,c_i)
    M,m = solver.optimal()
    a = pd.DataFrame([])
    a['M'] = M
    a['m'] = m
    a.to_csv('out2.csv',index=False)
    #M,m = solver.optimal_one_seg(dp4[0],c_i[0])
'''



