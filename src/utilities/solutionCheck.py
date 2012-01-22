import globalsettings

def _solutionCheck(dict1, dict2):
    
    maxErr = globalsettings.getSetting("maximalError")

    assert set(dict1.keys()) == set(dict2.keys())

    for k in dict1.keys():
       if not abs(dict1[k] - dict2[k]) < maxErr:
           return False
    return True

def solutionCheck(sola, solb):
    
    assert len(sola) == len(solb)
    
    for a in sola:
        foundMatch = False
        for b in solb:
            if _solutionCheck(a,b):
                foundMatch = True
        assert foundMatch

