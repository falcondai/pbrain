#!/usr/bin/python

import sys
import pickle

from scipy import zeros

def get_xform(filename):
    ascii_file = file(filename, 'r')
    print ascii_file
    while 1:
        ascii_line = ascii_file.readline().strip()
        print "ascii_line:" , ascii_line
        if (ascii_line[0:7] == 'Inverse'):
            # next few lines are matrix
            matrix_line1  = ascii_file.readline().strip().split()
            matrix_line2  = ascii_file.readline().strip().split()
            matrix_line3  = ascii_file.readline().strip().split()
            matrix_line4  = ascii_file.readline().strip().split()

            print matrix_line1, matrix_line2, matrix_line3
            xform = zeros((4,4), 'd')
            for i in range(0,4):
                print "i is ", i
                xform[0][i] = float(matrix_line1[i])
                xform[1][i] = float(matrix_line2[i])
                xform[2][i] = float(matrix_line3[i])
                xform[3][i] = float(matrix_line4[i])

            break
        elif not ascii_line:
            print "didn't find \'Inverse\'!"
            sys.exit(1)
    return xform
    

def run(filename1,filename2):

    xform = get_xform(filename1)

    print "xform is\n", xform

    scipy_mat = xform
    scipy_mat = scipy_mat.transpose()
    fh = file(filename2, 'w')
    pickle.dump(scipy_mat, fh)
    fh.close()
