from sklearn.metrics import r2_score

def linear(num, start, slope):
    #formatter = '{0:.3f}'
    if num < start:
        p_weight = 1.0
        d_weight = 0.0
    elif (num - start) * slope >= 1:
        p_weight = 0.0
        d_weight = 1.0
    else:
        d_weight = (num - start) * slope
        p_weight = 1 - d_weight
    #p_weight = float(formatter.format(p_weight))
    #d_weight = float(formatter.format(d_weight))
    return p_weight, d_weight

def geometric(num, start, decay_rate):
    #formatter = '{0:.3f}'
    if num < start:
        p_weight = 1.0
        d_weight = 0.0
    else:
        p_weight = (1.0 - decay_rate) ** (num - start)
        d_weight = 1 - p_weight
    #p_weight = float(formatter.format(p_weight))
    #d_weight = float(formatter.format(d_weight))
    return p_weight, d_weight

for i in range(1):
    pl, dl = linear(i+1,100,1/400)
    pg, dg = geometric(i+1,0,1/400)
    #print(f'{i+1}\t{pl:.4f}\t{pg}')
    #print('L', i+1, linear(i+1,10,0.05))
    #print('G', i+1, geometric(i+1,10,0.05))

ime = [12, 10, 4, 0]
exp = [6, 5, 2, 0]
print(r2_score(ime, exp))

'''
def weights(start,end,total):
    weights = {}
    inc = 1/(end - start)  #increment for calculating the linear weights
    for num in range(start):
        weights[num + 1] = [1,0]
    for num in range(end-start):
        num = num + 1 + start
        weights[num] = [1 - (inc * (num - start)),0 + (inc * (num - start))]
    for num in range(total - end):
        weights[num + 1 + end] = [0,1]
    return weights

    elif (num - start) * decay_rate >= 1:
        p_weight = 0.0
        d_weight = 1.0
'''
'''
length - 100
offset - 20
decay - 20 - 60
prox/dist - 60 - 100
'''
'''
def linear(num, start, end):
    if num < start:
        p_weight = 1
        d_weight = 0
    elif num > end:
        p_weight = 0
        d_weight = 1
    elif num >= start and num <= end:
        a = 1/(end-start)
        p_weight = 1 - (a*(num-start))
        d_weight = 0 + (a*(num-start))
    return p_weight, d_weight
'''
