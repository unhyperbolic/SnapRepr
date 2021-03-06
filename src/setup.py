from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

ext_modules=[
    Extension("algebra._pari",
              ["algebra/_pari.pyx"],
              libraries=["pari"]),
    Extension("linearCombinations.twoTerms",
              ["linearCombinations/twoTerms.pyx"])
    ]

#    Extension("linear_combinations.find_as_multiple_sum_diff",
#              ["linear_combinations/find_as_multiple_sum_diff.pyx"]),
#    Extension("linearCombinations._findAsSumOrDiff",
#              ["linearCombinations/_findAsSumOrDiff.pyx"])
#]

setup(
    name = "SnapRepr",
    cmdclass = {"build_ext": build_ext},
    ext_modules = ext_modules
)
