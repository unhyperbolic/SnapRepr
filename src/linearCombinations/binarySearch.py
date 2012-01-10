import globalsettings

def _binarySearchIndex(listOfDicts, key, value, firstIndex, lastIndex):
    
    # invariant: the index i of the last entry x in listOfDics 
    #            with x[key] < value
    #            is in [firstIndex, lastIndex)

    if firstIndex + 1 >= lastIndex:
        return firstIndex
    
    middleIndex = (firstIndex + lastIndex) / 2

    if listOfDicts[middleIndex][key] < value:
        return _binarySearchIndex(listOfDicts, key, value,
                                  middleIndex, lastIndex)
    else:
        return _binarySearchIndex(listOfDicts, key, value,
                                  firstIndex, middleIndex)

# returns the index of the last row x with x['key'] <= value
# (adding maximalError to it)

def binarySearchIndex(listOfDicts, key, value):
    return _binarySearchIndex(
        listOfDicts, key,
        value = value + globalsettings.getSetting("maximalError"),
        firstIndex = 0, lastIndex = len(listOfDicts))

def binarySearch(listOfDicts, key, value):
    
    index = binarySearchIndex(listOfDicts, key, value)

    maxError = globalsettings.getSetting("maximalError")

    if (index < len(listOfDicts) and
        abs(listOfDicts[index][key] - value) < maxError):

        return listOfDicts[index]

    return None

def matchingRowIndices(listOfDicts, key, value):

    index = binarySearchIndex(listOfDicts, key, value) + 1
    firstIndex = index
    lastIndex = index

    maxError = globalsettings.getSetting("maximalError")
    
    while (firstIndex > 0 and 
           abs(listOfDicts[firstIndex - 1][key] - value) < maxError):
        firstIndex -= 1

    return firstIndex, lastIndex

def matchingRows(listOfDicts, key, value):

    firstIndex, lastIndex = matchingRowIndices(listOfDicts, key, value)

    return [ listOfDicts[index] 
             for index in range(firstIndex, lastIndex) ]
