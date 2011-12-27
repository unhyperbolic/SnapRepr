import csv
import mpmath

conversionDict = {
    
    'Volume' : mpmath.mpf,
    'Tetrahedra' : int

    }

def readCensusTable(csvFile, sortKey = "Volume"):

    if isinstance(csvFile, str):
        csvFile = open(csvFile, 'rb')

    header = csv.reader(csvFile).next()

    csvFile.seek(0)

    censusVolumesReader = csv.DictReader(csvFile)
    
    def convertDict(d):

        global conversionDict

        for k, v in d.items():
            if conversionDict.has_key(k):
                d[k] = conversionDict[k](v)

        return d

    censusVolumesDict = [convertDict(d) for d in censusVolumesReader]
    
    if sortKey:
        censusVolumesDict.sort(key = lambda d: d[sortKey])
        
    return censusVolumesDict, header
