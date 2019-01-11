import numpy as np

class Reaction(object):
    def __init__(self, rate, first, second):
        self.rate = rate
        self.first = first
        self.second = second

with open('params/oxidise.data','r') as f:
    count = 0
    for line in f:
        count += 1
    N = count /3 
    f.seek(0)
    reactions = []
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



