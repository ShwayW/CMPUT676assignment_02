# The code for question 01 of assignment 02
# Shway Wang
# Nov. 27, 2022

## Imports
import numpy as np
from random import uniform

# initialize the input sequence
def init_input_seq(p_min, p_max, epsilon, n):
	seq = []
	for i in range(n):
		# exclude 0 when randomly initializing w
		w = uniform(0, epsilon)
		while (w == 0):
			w = uniform(0, epsilon)
			
		# initialize v randomly
		ratio = uniform(p_min, p_max)
		v = w * ratio
		
		# append v and w to seq
		assert(v/w >= p_min and v/w <= p_max)
		seq.extend([v, w])
	return seq


# the threshold function from pp. 25 of lec. 7
def threshold(p_min, p_max, y):
	beta = 1/(1 + np.log(p_max/p_min))
	if (y >= 0 and y < beta):
		return p_min
	elif (y >= beta and y <= 1):
		return p_min * np.exp(y / beta - 1)


# the threshold-based algorithm from pp. 17 of lec. 7
def threshold_based_alg(p_min, p_max, input_seq):
	y = 0
	decision_seq = []
	for tI in range(0, len(input_seq), 2):
		v = input_seq[tI]
		w = input_seq[tI + 1]
		p = threshold(p_min, p_max, y)
		if (v/w < p):
			x = 0
		elif (v/w >= p):
			x = 1
		y += w * x
		y = min(1, y)
		decision_seq.append(x)
	return decision_seq


if (__name__ == "__main__"):
	# the lower bound on threshold
	p_min = 0.1
	
	# the upper bound on threshold
	p_max = 2
	
	# the epsilon value
	epsilon = 0.005
	
	# number of iterations
	n = 1000
	
	# initialize the input sequence
	input_seq = init_input_seq(p_min, p_max, epsilon, n)
	
	# compute the answer by the threshold-based algorithm
	ans = threshold_based_alg(p_min, p_max, input_seq)
	
	print(ans)
	
	
	
	
