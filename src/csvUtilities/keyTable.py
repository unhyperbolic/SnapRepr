
# given [{'name': 'a', 'value' : 1}, {'name':'b', 'value' : 2}]
# returns {'a':{'name':'a', 'value':1}, 'b':{'name':'b', 'value':2}}

def keyTableUnique(listOfDicts, key):

    keyedTable = {}

    for d in listOfDicts:
        if d.has_key(key):
            keyedTable[d[key]] = d

    return keyedTable

# given [{'name': 'a', 'value' : 1}, {'name':'a', 'value' : 2}]
# returns {'a':set([{'name':'a', 'value':1}, {'name':'a', 'value':2}])}

def keyTable(listOfDicts, key):

    keyedTable = {}

    for d in listOfDicts:
        if d.has_key(key):
            if keyedTable.has_key(key):
                keyedTable[key].add(d)
            else:
                keyedTable[key] = set([d])

    return keyedTable
