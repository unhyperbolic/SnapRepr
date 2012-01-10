import mpmath
import globalsettings
from libc.stdlib cimport malloc, free

cdef extern from "math.h":
    double fabs(double x)

_epsilon = mpmath.mpf(1e-10)

# takes a list of dictionaries encoding a table of census manifolds with one key being volume
# saves an object with that list and a C array of doubles

cdef class censusTable:
    
    # C-style array of doubles
    cdef double *volumesArray
    cdef int leng

    # Python-style variable owned by Python
    cdef readonly censusTableDicts

    def __cinit__(self, dicts):
    
        # skip rows in dicts with zero volume and sort
        self.censusTableDicts = [d for d in dicts if d['Volume'] > _epsilon]
        self.censusTableDicts.sort(key = lambda x : x['Volume'])

        # allocate memory
        self.leng = len(self.censusTableDicts)
        self.volumesArray = <double *> malloc(self.leng * sizeof(double))

        for i, d in zip(range(self.leng),
                        self.censusTableDicts):
            self.volumesArray[i] = d['Volume']

    def __dealloc__(self):
        if self.volumesArray is not NULL:
            free(self.volumesArray)
	    
    def findAsTwoTermCombination(self, target, factorRange = 20):
        results = []
        for factor1 in range(1, factorRange + 1):
            for factor2 in range(1, factorRange + 1):
                results += self._findAsSumOrDifference(target,
                                                       factor1, factor2)
        return results

    def _findAsSumOrDifference(self, target, factor1, factor2):
        return (
              self._findAsSum(target, factor1, factor2)
            + self._findAsDifference(target, factor1, factor2))

    def _findAsSum(self, target, factor1, factor2):

        cdef int firstIdx, secondIdx
        cdef double cTarget
        cdef int cFactor1, cFactor2

        results = []

        cTarget = target
        cFactor1 = factor1
        cFactor2 = factor2

        firstIdx = 0
        secondIdx = self.leng - 1

        while (firstIdx < self.leng and secondIdx >= 0
               and 
               ((cFactor1 * self.volumesArray[firstIdx] + 1E-9)
                 < cFactor2 * self.volumesArray[secondIdx])):
            
            # test if we have match with C doubles
            if fabs(  cFactor1 * self.volumesArray[firstIdx]
                    + cFactor2 * self.volumesArray[secondIdx]
                    - cTarget) < 1E-9:

                # test if we have a match with mpmath.mpfloat
                if (abs(  factor1 * self.censusTableDicts[firstIdx]['Volume']
                        + factor2 * self.censusTableDicts[secondIdx]['Volume']
                        - target)
                    < globalsettings.getSetting("maximalError")):
                    results.append(
                        [ (factor1, self.censusTableDicts[firstIdx]),
                          (factor2, self.censusTableDicts[secondIdx])])
            
            # increase firstIdx or decrease secondIdx using C doubles
            if (  cFactor1 * self.volumesArray[firstIdx]
                  + cFactor2 * self.volumesArray[secondIdx] > cTarget):
                secondIdx -= 1
            else:
                firstIdx +=1

        return results

    def _findAsDifference(self, target, factor1, factor2):

        cdef int firstIdx, secondIdx
        cdef double cTarget
        cdef int cFactor1, cFactor2

        results = []

        cTarget = target
        cFactor1 = factor1
        cFactor2 = factor2

        firstIdx = 0
        secondIdx = 0

        while (firstIdx < self.leng and
               secondIdx < self.leng):

            # test if we have a match with C doubles
            if ((not firstIdx == secondIdx) and 
                fabs(  cFactor1 * self.volumesArray[firstIdx]
                     - cFactor2 * self.volumesArray[secondIdx]
                     - cTarget) < 1E-9):
                    
                # test if we have a match with mpmath.mpfloat
                if (abs(  factor1 * self.censusTableDicts[firstIdx]['Volume']
                        - factor2 * self.censusTableDicts[secondIdx]['Volume']
                        - target)
                    < globalsettings.getSetting("maximalError")):
                    
                        results.append(
                            [ ( factor1, self.censusTableDicts[firstIdx]),
                              (-factor2, self.censusTableDicts[secondIdx])])

            if (  cFactor1 * self.volumesArray[firstIdx]
                - cFactor2 * self.volumesArray[secondIdx] > cTarget):
                secondIdx += 1
            else:
                firstIdx += 1

        return results
