
from __future__ import print_function
import os
import sys
import math

TOLERANCE = 1e-08


# different compilers give different results so
# read in the reference file from the command line
if len(sys.argv) != 2 or not os.path.exists(sys.argv[1]):
    print("Usage: python BM.py <REFERENCE_FILE>")
    sys.exit(1)

ref_file = sys.argv[1]
new_file = "BMtest.txt"

print("Comparing files:")
print("  Ref: %s" % ref_file)
print("  New: %s" % new_file)

fo = open(ref_file, 'r')
ft = open(new_file, 'r')
lo = fo.readlines()
lt = ft.readlines()
fo.close()
ft.close()

column_tag = ['timestep',
              'last section',
              'photon number a',
              'photon number b',
              'probability of state g',
              'probability of state e',
              'probability of state f',
              'probability of atom decay',
              'probability of transmission a',
              'probability of reflection a',
              'probability of transmission b',
              'probability of reflection b']

count_diff = 0
for j in range(len(lo)):
    tempo = [float(k) for k in lo[j].split()]
    tempt = [float(k) for k in lt[j].split()]
    for l in range(2,12):
        if (tempo[l] != tempt[l]):
            # check if percentage difference is greater than tolerance
            if tempo[l] == 0 or math.fabs(100.0*(tempo[l] - tempt[l]) / tempo[l]) > TOLERANCE:
                print('*********************')
                print('Deviation detected!')
                print('*********************')
                print('Time step = ' + str(tempt[0]))
                print('Last section visited = ' + str(tempt[1]))
                print('Affected value = ' + column_tag[l])
                print('Percentage difference = '
                      + str(100.0*(tempo[l] - tempt[l])/tempo[l]) + '%')
                print(' ')
                count_diff += 1

if count_diff > 0:
    sys.exit(1)

print("Files match")
