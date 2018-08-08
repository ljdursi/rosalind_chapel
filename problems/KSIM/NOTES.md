For this one had to greatly trim the timing case, as otherwise
python code even with numba takes forever.

Note that string subsetting operations are quite expensive, 
so turned the strings into arrays of ascii values.
