#!/usr/bin/python

import os
import sys

files = os.listdir('.')

for f in files:
    
    # a trig file that has not _ordered or _unorderable

    if '.trig' in f:
        if not 'order' in f:

            print f
            
            prefix = f[:-5]
            unorderableF = prefix + "_unorderable.trig"

            for j in files:
                if j[:len(prefix)] == prefix:
                    if not '.trig' in j:
                        if not 'order' in j:
                            if unorderableF in files:
                                if len(sys.argv) > 1 and sys.argv[1] == 'move':
                                    os.rename(j, prefix+'_unorderable'+j[len(prefix):])
                                print j, "->", prefix+'_unorderable'+j[len(prefix):]
                            else:
                            
                                print "delete", j
                                if len(sys.argv) > 1 and sys.argv[1] == 'delete':
                                    os.remove(j)
