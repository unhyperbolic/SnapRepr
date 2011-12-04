#!/usr/bin/python

import sys
import os

if not len(sys.argv) == 2:
    print "Usage: detect_files_that_should_have_unorderable.py TRIANGULATION_FILE"
    sys.exit(1)

trig_file = sys.argv[1]

print trig_file
if "order" in trig_file:
    print "This file already went through order_triangulation.py"
    sys.exit(0)

if not trig_file[-5:] == '.trig':
    print "Not a .trig extenstion"
    sys.exit(1)

trig_file_base = trig_file[:-5]

if os.path.isfile(trig_file_base+"_ordered.trig"):
    print "Triangulation can be ordered"
else:
    print "Triangulation is unorderable"
    print "New triangulation file: ", trig_file_base+"_unorderable.trig"

    open(trig_file_base+"_unorderable.trig","w").write(open(trig_file,"r").read())

    trig_file_path = "/".join(trig_file.split("/")[:-1])

    trig = trig_file_base.split("/")[-1]

    for i in os.listdir("./"+trig_file_path):
        
        if trig == i[:len(trig)]:
            if not i[-5:] == ".trig":
                if not "order" in i:
                    target = trig + "_unorderable" + i[len(trig):]
                    src = i
                    if trig_file_path:
                        src = trig_file_path + "/" + src
                        target = trig_file_path + "/" + target

                    print "                         ", src, target
                    os.rename(src,target)

