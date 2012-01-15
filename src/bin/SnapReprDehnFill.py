#!/usr/bin/python

import sys

import fractions

import snappy

if not len(sys.argv) == 4:
    print "Usage: SnapReprDehnFill trig_file targetDir maxDehnFillCoefficient"
    sys.exit(1)

trig_file, targetDir, d = sys.argv[1:4]

d = int(d)

M = snappy.Manifold(trig_file)

baseFile = targetDir + '/' + trig_file.split('/')[-1][:-5]

for c in range(M.num_cusps()):
    for i in range(-d, d+1):
        for j in range(0, d+1):
            if j > 0 or i > 0:
                if fractions.gcd(i,j) == 1:
                    M = snappy.Manifold(trig_file)
                    M.dehn_fill((i,j),c)
                    k = M.solution_type()
                    if not ('flat' in k or 'degenerate' in k):
                        newFilename = baseFile + "(%d,%d)" % (i,j)
                        print newFilename + '.trig'
                        M.save(newFilename + ".trig")
