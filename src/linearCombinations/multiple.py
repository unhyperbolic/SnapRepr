# represents Volume of the manifold row / factor

class Multiple:
    def __init__(self, factor, row, key = None):
        self.factor = factor
        self.row = row
        self.key = key

    def asPair(self):
        return (self.factor, self.row)

    def getRepresentedVolume(self):
        return self.row['Volume'] / self.factor

    def __str__(self):
        if self.factor > 1:
            return self.row['Name'] + "/" + str(self.factor)
        else:
            return self.row['Name']

    def __repr__(self):
        return str(self)
