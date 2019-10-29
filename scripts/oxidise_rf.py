import numpy as np

class Reaction(object):
    def __init__(self, rate, first, second):
        self.rate = rate
        self.first = first
        self.second = second

reactions = []
with open('params/oxidise.data','r') as f:
    count = 0
    for line in f:
        count += 1
    N = count /3 
    f.seek(0)
    for i in range(N):
        rate = float(f.readline())
        first = [ int(j) for j in f.readline().split() ]
        second = [ int(j) for j in f.readline().split() ]

        a = Reaction(rate,first,second)
        a.lograte = np.log10(rate)
        a.first_alc_above = sum([i == 1 for i in a.first])
        a.first_alc_below = sum([i ==-1 for i in a.first])
        a.first_epo_above = sum([i == 2 for i in a.first])
        a.first_epo_below = sum([i ==-2 for i in a.first])
        a.second_alc_above = sum([i == 1 for i in a.second])
        a.second_alc_below = sum([i ==-1 for i in a.second])
        a.second_epo_above = sum([i == 2 for i in a.second])
        a.second_epo_below = sum([i ==-2 for i in a.second])
        reactions += [a]

attributes = ['first_alc_above','first_alc_below','first_epo_above',
              'first_epo_below','second_alc_above','second_alc_below',
              'second_epo_above','second_epo_below']

X = np.zeros((len(reactions),len(attributes)))
for ri,r in enumerate(reactions):
    for ai,a in enumerate(attributes):
        X[ri,ai] = getattr(r,a)
Y = np.zeros(len(reactions))
for ri,r in enumerate(reactions):
    Y[ri] = getattr(r,'lograte')


def init_random_forest():
    from sklearn.ensemble import RandomForestRegressor
    rf_reg = RandomForestRegressor(n_estimators=500, max_depth = 4)
    rf_reg.fit(X,Y)
    return rf_reg



def fit_empirical(reactions, a,b,c,d,e,f,g,h):
    m = [a,b,c,d,e,f,g,h]
    rates = []
    try:
        len(reactions) 
    except TypeError: reactions = [reactions]
    for reaction in reactions:
        steric = 0
        polar = 0
        hbond = 0
        edge = 0

        for state in reaction.first:
            if state == 1: steric += 1
            elif state == 2: steric += m[6]
            if abs(state) == 1: polar += 1
            if abs(state) == 2: polar += m[7]
            
            #if state == 1: hbond =
        
        for state in reaction.second:
            if state == 1: hbond += 1
            if state == 3: edge = 1
        
        steric = (  m[0]    * steric 
                  + m[1] * steric * steric)
        polar  = (  m[2] * polar
                  + m[3] * polar * polar)
        hbond  = (  m[4] * hbond
                  + m[5] * hbond * hbond) 

        if edge:
            rate = 1
        else:
            #rate = 10 ** (steric + polar + hbond)
            rate = steric + polar + hbond
        rates += [rate]
    return rates


if __name__ == "__main__":
    
    rf = init_random_forest()
    
    p0 = [ -3.867, 0.185,
              23.169, -5.138,
              11.648, -4.413,
              1, 0.633]

    from scipy.optimize import curve_fit

    print p0
    popt, pcov = curve_fit(fit_empirical, reactions, Y)
    print popt
    print pcov

    for i in range(len(reactions)):
        print Y[i], fit_empirical(reactions[i],*popt)[0], rf.predict([X[i]])[0]

    fi = rf.feature_importances_
    for i,j in zip(attributes,fi):
       print i,j
