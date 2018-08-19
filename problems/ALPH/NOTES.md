Again, trying to return any sort of non-primitive data type from a
recursive function requires creating (and then accessing) extraneous
records - here, returning an associative array in Tree.sequences()

Also, having to fight against promotion flattening is just absurd in
a technical computing language
