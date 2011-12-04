#!/usr/bin/python

import subprocess
import sys
import time
import random

class process:
    def __init__(self, psOutLine):
        pid, rss, self.state, coresize = psOutLine.split()
        self.pid = int(pid)
        self.rss = int(rss)
        self.coresize = int(coresize)
    def __repr__(self):
        return "process(pid = %d, rss = %d,, coresize = %d state = '%s')" % (self.pid, self.rss, self.coresize, self.state)

def get_swap_rate():
    l = subprocess.Popen("vmstat 1 2", shell = True, stdout = subprocess.PIPE)
    k = l.stdout.read()

    k = int(k.split()[44]) + int(k.split()[45])
    return k

def get_processes():
    l = subprocess.Popen("ps -C magma.exe -o pid,rss,s,sz",
                         shell = True, stdout = subprocess.PIPE)

    k = l.stdout.read().split('\n')[1:-1]
    m = [process(l) for l in k]
    return m
    
def get_running_processes(p):
    return [x for x in p if 'R' in x.state]

def get_disc_processes(p):
    return [x for x in p if 'D' in x.state]

def get_stopped_proccesses(p):
    return [x for x in p if 'T' in x.state]

def get_running_or_disc_processes(p):
    return [x for x in p if ('R' in x.state or 'D' in x.state)]

def get_process_most_ram(p):
    m = [x for x in p]
    m.sort(key = lambda x: x.rss)
    if len(m):
        return m[-1]
    return None

def get_total_memory(p):
    return sum([x.rss for x in p])

def STOP_process(p):
    assert isinstance(p, process)
    print "                Stopping", p.pid
    subprocess.Popen("kill -STOP %d" % p.pid, shell = True)

def CONT_process(p):
    assert isinstance(p, process)
    print "                Running", p.pid
    subprocess.Popen("kill -CONT %d" % p.pid, shell = True)

d = 10 * [ 0 ]

while True:
    p = get_processes()
    s = get_swap_rate()
    m = get_total_memory(p)
    print "============================"
    print "Swap rate: %d      Memory: %d" % (s, m)
    print "Processes:"
    for x in p:
        print "    ", x

    dp = get_disc_processes(p)

    d = d[1:] + [len(dp)]

    print d, sum(d)

    rdp = get_running_or_disc_processes(p)

    if sum(d) > 4:
        if len(rdp) > 1 and m > 2500000:
            STOP_process(rdp[-1])
    elif (not random.randint(0,15)) or m < 1000000:
        stp = get_stopped_proccesses(p)
        print stp
        if stp:
            CONT_process(stp[-1])
    time.sleep(10)


while True:
    p = get_processes()
    s = get_swap_rate()
    m = get_total_memory(p)
    print "============================"
    print "Swap rate: %d      Memory: %d" % (s, m)
    print "Processes:"
    for x in p:
        print "    ", x
    if s > 200:
        print "==> Swapping"
        p1 = get_process_most_ram(p)
        print "          MOST RAM PID", p1
        for x in p:
            if x == p1:
                CONT_process(x)
            else:
                STOP_process(x)
    if m < 2500000 or len(p) < 2:
        print "==> Safe"
        for x in p:
            CONT_process(x)
    if not get_running_processes(p):
        print "==> No process running"
        CONT_process(p[0])

    time.sleep(2)

