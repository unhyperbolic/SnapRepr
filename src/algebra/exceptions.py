class NumericalError(Exception):
    def __init__(self, val, msg):
        self.val = val
        self.msg = msg
    def __str__(self):
        return "NumericalError(val = %s) : %s" % (self.val, self.msg)

