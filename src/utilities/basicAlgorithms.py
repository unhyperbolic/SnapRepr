import operator

def safeDictLookup(d, k, default = None):
    if d.has_key(k):
        return d[k]
    else:
        return default

def indexedIterable(iterable):
    """
    >>> for index, item in indexedIterable(["Zero", "One", "Two", "Three"]):
    ...     print index, item
    0 Zero
    1 One
    2 Two
    3 Three
    """

    currentIndex = 0
    for item in iterable:
        yield currentIndex, item
        currentIndex += 1

# take a dictionary and turn it into a tuple of pairs sorted by keys

def dictToOrderedTupleOfPairs(d):
    """
    dictToOrderedTupleOfPairs(
    ...      { 'key3':'value3', 'key1':'value1', 'key2':'value2' })
    (('key1', 'value1'), ('key2', 'value2'), ('key3', 'value3'))
    """

    l = d.items()
    l.sort(key = lambda x:x[0])
    return tuple(l)

# given a list of dictionaries, combine values of the different
# dictionaries having the same key using combineFunction.

def combineDicts(listOfDicts, combineFunction):
    """
    >>> d = combineDicts(
    ...      [ {'key1': 1, 'key2': 2},
    ...        {'key1': 1} ],
    ...      combineFunction = operator.add)
    >>> d['key1']
    2
    >>> d['key2']
    2
    """

    result = {}
    for aDict in listOfDicts:
        for k, v in aDict.items():
            if result.has_key(k):
                result[k] = combineFunction(result[k], v)
            else:
                result[k] = v
    return result
