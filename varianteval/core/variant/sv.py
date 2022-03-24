class SVCallset:
    def __init__(self):
        self.calls = []
        # interval tree index etc

    def add(self, sv):
        self.calls.append(sv)

    def __len__(self):
        return len(self.calls)

    def __str__(self):
        pass

    def report(self):
        pass

class SV:
    def __init__(self, sv_type):
        self.sv_type = sv_type

    def __str__(self):
        pass

class Adjacency:
    def __init__(self):
        pass

    def __str__(self):
        pass

class Breakpoint:
    def __init__(self):
        pass

    def __str__(self):
        pass



