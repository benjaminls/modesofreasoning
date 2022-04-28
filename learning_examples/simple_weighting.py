#%%
import numpy as np
import matplotlib.pyplot as plt

# This example is transcribed from the video:
# https://www.youtube.com/watch?v=3dn7-vaIP1g&ab_channel=GarethTribello
# 
# I do not own, nor no claim to be the original author, of this code. 
# All credit goes to Gareth Tribello and the video listed above. 
# 

nsamples, wtot, xvals, yvals = 200, 0, np.linspace(1, nopts, nopts), np.zeros(opts)
# yvals -> array to hold histogram, initially all zeros
# 
for i in range(nsamples):
    # generate random variables (rv) and their corresponding weights
    sample, w = weighted_discrete_rv()
    # accumulate total weight in bin that sample falls within
    # by doing this, we are using yvals to accumulate the total
    #   weight that falls in that bin rather than the total number
    #   of points that fall in that bin, which is what we would have
    #   accumulated if the histogram was unweighted
    yvals[int(sample)] = yvals[int(sample)] + w
    wtot = wtot + w

# normalize yvals using total weight across all bins
for i in range(len(yvals)):
    yvals[i] = yvals[i] / wtot
plt.bar(xvals, yvals, 0.1)
plt.show()