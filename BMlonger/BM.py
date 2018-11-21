#from pylab import *

fo = open('./BMoriginal.txt','r')
ft = open('./BMtest.txt','r')
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

for j in range(len(lo)):
    tempo = [float(k) for k in lo[j].split()]
    tempt = [float(k) for k in lt[j].split()]
    for l in range(2,12):
        if (tempo[l] != tempt[l]):
            print('*********************')
            print('Deviation detected!')
            print('*********************')
            print('Time step = ' + str(tempt[0]))
            print('Last section visited = ' + str(tempt[1]))
            print('Affected value = ' + column_tag[l])
            print('Percentage difference = '
                  + str(100.0*(tempo[l] - tempt[l])/tempo[l]) + '%')
            print(' ')
