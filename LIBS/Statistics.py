# alternatively use Numpy or SciPy 
# this is useful, because it does not rely on neither of them (good for salviati)

## {{{ http://code.activestate.com/recipes/511478/ (r1)
import math
import functools

def percentile(N, percent, key=lambda x:x):
    """
    Find the percentile of a list of values.

    @parameter N - is a list of values. Note N MUST BE already sorted.
    @parameter percent - a float value from 0.0 to 1.0.
    @parameter key - optional key function to compute value from each element of N.

    @return - the percentile of the values
    """
    if not N:
        return None
    k = (len(N)-1) * percent
    f = math.floor(k)
    c = math.ceil(k)
    if f == c:
        return key(N[int(k)])
    d0 = key(N[int(f)]) * (c-k)
    d1 = key(N[int(c)]) * (k-f)
    return d0+d1

# median is 50th percentile.
median = functools.partial(percentile, percent=0.5)
median([1,2,3,4,5]) # output: 3
## end of http://code.activestate.com/recipes/511478/ }}}
