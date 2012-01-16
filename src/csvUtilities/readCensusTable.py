import csv
import mpmath

conversionDict = {
    
    'Volume' : mpmath.mpf,
    'Tetrahedra' : int,
    'InvariantTraceFieldDegree' : int,
    'SL(N,C)' : int

    }

header = ["Manifold",
          "File",
          "Ordered",
          "Tetrahedra",
          "Cusps",
          "SL(N,C)",
          "Obstruction Class",
          "Index Component",
          "Number Components",
          "Dimension",
          "Additional Eqns",
          "Number Solutions",
          "Index Solution",
          "Warning",
          "Volume",
          "CS",
          "LinearCombinations",
          "CPUTIME"]

class CensusTable:
    def __init__(self, listOfDicts, header):
        assert isinstance(listOfDicts, list)
        assert isinstance(header, list)
        if len(listOfDicts) > 0:
            assert isinstance(listOfDicts[0], dict)

        self.listOfDicts = listOfDicts
        self.header = header

def readCensusTable(csvFile,
                    sortKey = None,
                    readHeaderFromFile = True):

    def convertDict(d):

        global conversionDict

        for k, v in d.items():
            if conversionDict.has_key(k):
                try:
                    d[k] = conversionDict[k](v)
                except:
                    d[k] = None

        return d

    if isinstance(csvFile, str):
        csvFile = open(csvFile, 'rb')

    if readHeaderFromFile:
        censusVolumesReader = csv.DictReader(csvFile)
    else:
        censusVolumesReader = csv.DictReader(csvFile,
                                             fieldnames = header)

    censusVolumesDict = [convertDict(d) for d in censusVolumesReader]

    if sortKey:
        censusVolumesDict.sort(key = lambda d: d[sortKey])
        
    return CensusTable(listOfDicts = censusVolumesDict,
                       header = censusVolumesReader.fieldnames)
