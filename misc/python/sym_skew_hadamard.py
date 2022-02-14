#!/usr/bin/env python3

# Script to get the symmetric skew hadamard matrices by the method described
# on the link below without doing errors.

# https://www.rangevoting.org/SkewHad.html

N = 19; # N has to be a prime. Will print a matrix 1 dimension bigger

squares_mod = set([]);

for v in range (N):
  squares_mod.add ((v * v) % N)

print (squares_mod)

mtx = [];
for _ in range (N):
  mtx.append ([False] * N)

# Build matrix
# row = i, col = j
for i in range (1, N + 1):
  for j in range (1, N + 1):
    mtx[i - 1][j - 1] = (j - i) in squares_mod

# Reflect (brute force, as it has to be a skew matrix)
for i in range (0, N):
  for j in range (i + 1, N):
        mtx[j][i] = not mtx[i][j]

# Print
for col in range (N + 1):
  print ('+', end = ' ')
print ('')

for i in range (0, N):
  print ('-', end = ' ')
  for j in range (0, N):
    print ('+' if mtx[i][j] == True else '-', end = ' ')
  print ('')
