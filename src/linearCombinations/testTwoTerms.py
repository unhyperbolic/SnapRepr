from csvUtilities.readCensusTable import readCensusTable
import globalsettings
import mpmath
from linearCombinations import binarySearch, filterRepresentativeMfds, twoTerms
from linearCombinations.formatLinearCombination import formatLinearCombination

mpmath.mp.dps = 70
globalsettings.setSetting("maximalError", mpmath.mpf("0.1") ** 50)
globalsettings.setSetting("maximalErrorDigits", 50)

testCensusTablePath = globalsettings.getSetting("testCensusTablePath")
censusTableFile = testCensusTablePath + "/exampleCensusTable.csv"
censusTable = readCensusTable(censusTableFile, sortKey = "Volume")

representatives = filterRepresentativeMfds.filterRepresentativeMfds(censusTable.listOfDicts)

def checkBinarySearchResult(vol, resultNames):

    if isinstance(vol, str):
        vol = mpmath.mpf(vol)

    rows = binarySearch.matchingRows(censusTable.listOfDicts, "Volume", vol)
    names = [ row['Name'] for row in rows ]

    assert set(resultNames) == set(names), Exception(
        "Expected: %s\nGot %s" % (resultNames, names))

def testBinarySearch():

    checkBinarySearchResult(
        "2.0298832128193072500424051085490405718833786150605995840349782135",
        ['4_1', '4a1', '4_DT[dadbcda]',
         'm010(4,1)', 'm010(-4,3)', 'm011(-3,1)', 'm015(6,1)', 'm016(-2,3)',
         'm019(3,2)', 'm036(-3,2)',
         'm003', 'm004'])

    checkBinarySearchResult(
        "4.0568602242368201441819241840176594471606643770938095602265600154",
        ['m199', 'm200', 'm201'])

    checkBinarySearchResult(
        "4.0597664256386145000848102170980811437667572301211991680699564271",
        ['11^2_DT[kbchdEfcHiJKaBG]', '6^2_2', '6^2_DT[fbccdefacb]',
         'm202', 'm203', 'm206', 'm207', 'm208'])
    
    checkBinarySearchResult(
        "4.9108330387398016772330757358262743608287703174737661838818269790",
        ["m395"])

    checkBinarySearchResult(
        "5.0747080320482681251060127713726014297084465376514989600874455338",
        ["m410", "m412"])

    checkBinarySearchResult("0.1",[])

    checkBinarySearchResult("4.12312313132", [])
    
    checkBinarySearchResult("10000.0", [])

def testRepresentatives():

    assert representatives[0:5] == [
        { 'Volume': mpmath.mpf('0.9427073627769277209212996030922116475903271057668831590145067757529341787'),
          'Representatives': '(m003(-3,1)=m015/3=5a1/3)'},
        { 'Volume': mpmath.mpf('0.9813688288922320880914521897944270682381643219063124386426041997774204779'),
          'Representatives': '(m003(-2,3)=m019/3)'}, 
        { 'Volume': mpmath.mpf('1.014941606409653625021202554274520285941689307530299792017489106776597468'),
          'Representatives': '(m007(3,1)=m003/2=4_1/2)'}, 
        { 'Volume': mpmath.mpf('1.263709238658043655884716346808091104479410748485054933690581587601325651'),
          'Representatives': '(m003(-4,3)=m155/3)'}, 
        { 'Volume': mpmath.mpf('1.284485300468354442460337084868711258912418773183928092149437015169645086'),
          'Representatives': '(m004(6,1)=m006/2=7_4/4)'}
        ]

def testTwoTerms():
    
    table = twoTerms.censusTable(representatives)

    vol = representatives[3]['Volume'] * 3 + representatives[5]['Volume'] * 7

    assert '3 * (m003(-4,3)=m155/3) +7 * m004(1,2)' in [
        formatLinearCombination(l) for l in table.findAsTwoTermCombination(vol)]

    vol = representatives[5]['Volume'] * 4 - representatives[7]['Volume'] * 5

    assert '4 * m004(1,2) -5 * m003(-4,1)' in [
        formatLinearCombination(l) for l in table.findAsTwoTermCombination(vol)]
    

def testAll():
    testBinarySearch()
    testRepresentatives()
    testTwoTerms()
