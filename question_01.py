# The code for question 01 of assignment 02
# Shway Wang
# Nov. 27, 2022

## Imports
import numpy as np

# the threshold function from pp. 25 of lec. 7
def threshold(p_min, p_max, y):
	beta = 1/(1 + np.log(p_max/p_min))
	if (y >= 0 and y < beta):
		return p_min
	elif (y >= beta and y <= 1):
		return p_min * np.exp(y / beta - 1)

# the threshold-based algorithm from pp. 17 of lec. 7
def threshold_based_alg(p_min, p_max, T):
	for t in range(T):
		if (v/w < p):
			x_t = 0
		else:
			x_t = 1
		y_t = 

if (__name__ == "__main__"):
	# the lower bound on threshold
	p_min = 0
	
	# the upper bound on threshold
	p_max = 2
	
	# number of iterations
	T = 100
	
	# compute the answer by the threshold-based algorithm
	ans = threshold_based_alg(p_min, p_max, T)
	
	
	
	
