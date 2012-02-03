import sys
import signal
import subprocess

_pariWrapper = "SnapRepr_PariProcess.py"
_pariProcess = None

def _alarmHandler(signum, frame):
    _destroyPariProcess()
    from algebra import pari
    raise pari.TimeoutAlarm()

# destroy the pari process
def _destroyPariProcess():
    global _pariProcess
    _pariProcess.kill()
    _pariProcess = None

# start a new pari process
def _createNewPariProcess():
    global _pariProcess
    _pariProcess = subprocess.Popen(
        _pariWrapper,
        stdin=subprocess.PIPE, stdout=subprocess.PIPE)

# get a working pari process
def _getPariProcess():
    global _pariProcess

    # if no pari process exists, start one
    if not _pariProcess:
        _createNewPariProcess()
    else:
        # if a pari process exists but terminated, start one
        _pariProcess.poll()
        if not _pariProcess.returncode is None:
            _createNewPariProcess()
            
    return _pariProcess

# run a pari command in a remote pari process

def remotePariEval(s, timeout = 60):

    # get a working pari process
    pariProcess = _getPariProcess()

    if not timeout is "process":
        signal.signal(signal.SIGALRM, _alarmHandler)
        signal.alarm(timeout)

    pariProcess.stdin.write(repr(s)+'\n')
    pariProcess.stdin.flush()
    pariProcess.stdout.flush()
    resultStr = pariProcess.stdout.readline().strip()
    if not resultStr:
        from algebra import pari
        raise pari.PariError
    result = eval(resultStr)
    if not timeout is "process":
        signal.alarm(0) # reset alarm
    return result
