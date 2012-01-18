import re
from linearCombinations import binarySearch
from utilities.basicAlgorithms import safeDictLookup
from linearCombinations.multiple import Multiple

def matchClosed(row):
    match = re.match('[msvt]\d+(\(-?\d+,-?\d+\))+$', row['Name'])

    if not match:
        return None
    
    degree = safeDictLookup(row,"InvariantTraceFieldDegree", 1000)
    if not degree:
        degree = 1000

    return (row["Tetrahedra"], degree, row["Name"])

def matchCusped(row):
    match = re.match('[msvt]\d+$', row['Name'])

    if not match:
        return None

    degree = safeDictLookup(row,"InvariantTraceFieldDegree", 1000)
    if not degree:
        degree = 1000

    return (row["Tetrahedra"], degree, row["Name"])

def matchLink(row):
    crossings = None
    components = None
    index = None

    preference = 1000

    matchRolfsen = re.match("(\d+)(\^(\d+))?_(\d+)", row['Name'])
    if matchRolfsen:
        crossings, s, components, index = matchRolfsen.groups()
        index = int(index)
        preference = 0

    matchKnot = re.match("(\d+)([an])(\d+)", row['Name'])
    if matchKnot:
        crossings, alt, index = matchKnot.groups()
        index = int(index)
        components = 1
        if alt == "a":
            preference = 1
        else:
            preference = 2

    matchMorwen = re.match("(\d+)(\^(\d+))?_(DT\[\w+\])", row['Name'])
    if matchMorwen:
        crossings, dum, components, index = matchMorwen.groups()
        preference = 3
    
    if not crossings:
        return None

    crossings = int(crossings)
    if components:
        components = int(components)
    components = 1

    degree = safeDictLookup(row, "InvariantTraceFieldDegree", 1000)
    if not degree:
        degree = 1000
    
    return (components, crossings, degree, row["Tetrahedra"], index, preference)

def findOneCanonicalRepresentative(listOfMultiples, matchFunction):
    
    # make into triple containing value of matchFunction
    matches = [
        Multiple(factor = multiple.factor, 
                 row = multiple.row,
                 key = matchFunction(multiple.row))
        for multiple in listOfMultiples]

    # filter those with a match
    matched = [m for m in matches if m.key]

    # sort by match
    matched.sort(key = lambda x:x.key)

    return matched[0:1]

# returns pairs

def defaultFindCanonicalRepresentatives(listOfPairs):

    result = []

    for matchFunction in [matchClosed, matchCusped, matchLink]:
        
        result += findOneCanonicalRepresentative(listOfPairs, matchFunction)

    return result

def defaultPrintRepresentatives(multiples):
    
    strs = [str(multiple) for multiple in multiples]
    if len(strs) == 1:
        return strs[0]
    else:
        return "(%s)" % ("=".join(strs))

# Assumes that listOfDicts is sorted by Volume

def filterRepresentativeMfds(listOfDicts, 
                             findCanonicalRepresentativesFunction = 
                                      defaultFindCanonicalRepresentatives,
                             formatFunction = defaultPrintRepresentatives):

    assert isinstance(listOfDicts, list)
    assert isinstance(listOfDicts[0], dict)

    # The largest volume
    maxVolume = listOfDicts[-1]['Volume']

    # Volumes and multiple of Volumes already covered
    markAsMultiple = len(listOfDicts) * [ False ]

    # the result
    representatives = [ ]

    # iterate through all manifolds not covered yet
    for i in range(len(listOfDicts)):
        if not markAsMultiple[i]:

            # current volume
            volume = listOfDicts[i]['Volume']

            # manifolds with that volume or a multiple of that volume
            mfdsWithThatVolumeOrMultiple = []

            # iterate through multiples of the volume
            for m in range(1, int(maxVolume/volume) + 2):

                firstIndex, lastIndex = binarySearch.matchingRowIndices(
                    listOfDicts, 'Volume', m * volume)

                # add manifolds of that volume to the list
                for j in range(firstIndex, lastIndex):
                    mfdsWithThatVolumeOrMultiple.append(
                        Multiple(factor = m, row = listOfDicts[j]) )
                    markAsMultiple[j] = True

            # find the canonical representatives

            canonicalMfdsWithThatVolumeOrMultiple = (
                findCanonicalRepresentativesFunction(
                    mfdsWithThatVolumeOrMultiple))

            # add the result to the list

            if canonicalMfdsWithThatVolumeOrMultiple:
                representatives.append(
                    { 'Volume' : volume ,
                      'Representatives' : 
                         formatFunction(
                            canonicalMfdsWithThatVolumeOrMultiple)})
            else:
                raise Exception, "Unmatched manifold for " + str(mfdsWithThatVolumeOrMultiple)

    return representatives

