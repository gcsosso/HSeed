import os.path
import math
import sys
import numpy

def Lateral_Prob( infilename,  \
                  outfilename, \
                  xmin,        \
                  xmax,        \
                  ymin,        \
                  ymax,        \
                  dx,          \
                  dy):
    """ 
    [CODE FROM STEVE]

    infilename  : name of input file. Two column whitespace
                  separated. Accepts "#" as comment.
    outfilename : name of output file, given as x, y, z data.
    xmin, ymin  : specifies lower limit of distribution.
    xmax, ymax  : specifies upper limit of distribution.
    dx, dy      : bin widths.

    """

    if os.path.isfile(outfilename):
        # don't allow overwriting of files
        print "Not allowed to overwrite files!"
        print "attempted file to overwrite: %s" % outfilename
        sys.exit(1)
    else:
        infile  = open(infilename ,'r')
        outfile = open(outfilename,'w')
        
        nbin_x = math.trunc((xmax-xmin)/dx) + 1
        nbin_y = math.trunc((ymax-ymin)/dy) + 1

        prob = numpy.empty( (nbin_x, nbin_y) )
        prob.fill(10**-10)
        
        ndat = 0.0

        outfile.write("# xmin = " + str(xmin) + ", ymin = " + str(ymin) + 
                      ", xmax = " + str(xmax) + ", ymax = " + str(ymax) +
                      ", dx = "   + str(dx)   + ", dy = "   + str(dy) + " Angstrom\n")

        # Loop over all lines in the input file
        for line in infile:
            info = line.split()

            print "BEEP"
            # check for comment line
            if info[0][0] == "#":
                continue
            
            # increase the total number of data points by 1
            ndat += 1.0

            x = float(info[0])
            y = float(info[1])

            # calculate the appropriate bin
            xbin = math.trunc( (x-xmin)/dx )
            ybin = math.trunc( (y-ymin)/dy )

            print xbin, ybin

            if xbin >= nbin_x-1:
              print "[ERROR]. xbin out of bounds."
              print xbin, nbin_x
              sys.exit()
            if ybin >= nbin_y-1:
              print "[ERROR]. ybin out of bounds."
              print ybin, nbin_y
              sys.exit()

            prob[xbin][ybin] += 1.0
        
        # normalize the distribution
        prob = prob/ndat

        logP = numpy.empty( (nbin_x, nbin_y) )
        logP = (-1)*numpy.log(prob)

        zero_logP = numpy.amin(logP)

        # print the probability distribution to the outfile
        for i in range(nbin_x):

            bc_x = xmin + i*dx

            for j in range(nbin_y):

                bc_y = ymin + j*dy

                logP[i][j] -= zero_logP 

                outfile.write(str(bc_x) + " " + str(bc_y) + " " + str(prob[i][j]) + " " + str(logP[i][j]) + "\n")

            outfile.write("\n")
