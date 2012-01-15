def formatLinearCombination(l):

    # contains a list of pairs (factor, row of the database)

    result = ""

    assert isinstance(l, list)
    for factor, row in l:
        if factor == +1:
            result += "+"
        elif factor == -1:
            result += "-"
        else:
            result += "%+d * " % factor
        result += row['Representatives'] + " "

    if result[0] == "+":
        return result[1:].strip()
    else:
        return result.strip()
