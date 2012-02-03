#!/usr/bin/python

import sys

import fractions

import snappy

if not len(sys.argv) == 5:
    print "Usage: SnapReprDehnFill trig_file targetDir maxDehnFillCoefficientP maxDehnFillCoefficientQ"
    sys.exit(1)

trig_file, targetDir, pDehn, qDehn = sys.argv[1:5]

pDehn = int(pDehn)
qDehn = int(qDehn)

M = snappy.Manifold(trig_file)

baseFile = targetDir + '/' + trig_file.split('/')[-1][:-5]

for c in range(M.num_cusps()):
    for i in range(-pDehn, pDehn+1):
        for j in range(0, qDehn+1):
            if j > 0 or i > 0:
                if fractions.gcd(i,j) == 1:
                    M = snappy.Manifold(trig_file)
                    M.dehn_fill((i,j),c)
                    k = M.solution_type()
                    if not ('flat' in k or 'degenerate' in k):
                        newFilename = baseFile + "".join(
                            ["(%d,%d)" % (i,j) if x == c else "(0,0)" for x in range(M.num_cusps())])
                        print newFilename + '.trig'
                        M.save(newFilename + ".trig")
