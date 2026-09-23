
import group_methods as grp

#rewriting rules for word problem

shortenings = {
    # -------------------------------------------------
    # Cancellation
    # -------------------------------------------------
    (1, -1): [],
    (2, -2): [],
    (3, -3): [],
    (4, -4): [],
    (5, -5): [],
    (6, -6): [],
    (7, -7): [],
    (-1, 1): [],
    (-2, 2): [],
    (-3, 3): [],
    (-4, 4): [],
    (-5, 5): [],
    (-6, 6): [],
    (-7, 7): [],
    # -------------------------------------------------
    # Positive: [i,j] -> [-k]
    # -------------------------------------------------
    (1, 4): [-2],
    (4, 2): [-1],
    (2, 1): [-4],
    #
    (3, 3): [-1],
    (3, 1): [-3],
    (1, 3): [-3],
    #
    (5, 5): [-4],
    (5, 4): [-5],
    (4, 5): [-5],
    #
    (6, 6): [-2],
    (6, 2): [-6],
    (2, 6): [-6],
    #
    (7, 1): [-6],
    (1, 6): [-7],
    (6, 7): [-1],
    #
    (7, 2): [-5],
    (2, 5): [-7],
    (5, 7): [-2],
    #
    (7, 4): [-3],
    (4, 3): [-7],
    (3, 7): [-4],
     # -------------------------------------------------
    # Negative: [-i,-j] -> [k]
    # -------------------------------------------------
    (-4, -1): [2],
    (-2, -4): [1],
    (-1, -2): [4],
    #
    (-3, -3): [1],
    (-1, -3): [3],
    (-3, -1): [3],
    #
    (-5, -5): [4],
    (-4, -5): [5],
    (-5, -4): [5],
    #
    (-6, -6): [2],
    (-2, -6): [6],
    (-6, -2): [6],
    #
    (-1, -7): [6],
    (-6, -1): [7],
    (-7, -6): [1],
    #
    (-2, -7): [5],
    (-5, -2): [7],
    (-7, -5): [2],
    #
    (-4, -7): [3],
    (-3, -4): [7],
    (-7, -3): [4],
}

rewritings = {
    (-1, 2): [6, -5],
    (-1, 3): [3, -1],
    (-1, 4): [6, -3],
    (-1, 5): [4, -7],
    (-1, 6): [4, -6],
    (-1, 7): [3, -4],
    #
    (-2, 1): [5, -6],
    (-2, 3): [1, -7],
    (-2, 4): [5, -3],
    (-2, 5): [1, -5],
    (-2, 6): [6, -2],
    (-2, 7): [6, -1],
    #
    (-3, 1): [1, -3],
    (-3, 2): [7, -1],
    (-3, 4): [3, -2],
    (-3, 5): [7, -5],
    (-3, 6): [3, -7],
    (-3, 7): [1, -4],
    #
    (-4, 1): [3, -6],
    (-4, 2): [3, -5],
    (-4, 3): [2, -3],
    (-4, 5): [5, -4],
    (-4, 6): [2, -7],
    (-4, 7): [5, -2],
    #   
    (-5, 1): [7, -4],
    (-5, 2): [5, -1],
    (-5, 3): [5, -7],
    (-5, 4): [4, -5],
    (-5, 6): [7, -6],
    (-5, 7): [4, -2],
    #
    (-6, 1): [6, -4],
    (-6, 2): [2, -6],
    (-6, 3): [7, -3],
    (-6, 4): [7, -2],
    (-6, 5): [6, -7],
    (-6, 7): [2, -1],
    #
    (-7, 1): [4, -3],
    (-7, 2): [1, -6],
    (-7, 3): [4, -1],
    (-7, 4): [2, -5],
    (-7, 5): [2, -4],
    (-7, 6): [1, -2],
}

class A2Word:
    
    def __init__(self, letters, normal=False):
        self.letters = list(letters)
        self.normal = normal
        if not self.normal:
            self.normal_form()

    def __repr__(self):
        return f"A2Word({self.letters})"

    @classmethod
    def identity(cls):
        return cls([], normal=True)

    def normal_form(self):
        if self.normal:
            return self
        else:
            changed = True
            while changed:
                changed = False
                i = 0
                while i < len(self.letters) - 1:
                    pair = (self.letters[i], self.letters[i + 1])
                    if pair in shortenings:
                        self.letters[i:i + 2] = shortenings[pair]
                        changed = True
                        i = max(i - 1, 0)
                    elif pair in rewritings:
                        self.letters[i:i + 2] = rewritings[pair]
                        changed = True
                        i = max(i - 1, 0)
                    else:
                        i += 1
            self.normal = True
            return self

    def __mul__(self, other):
        return A2Word(self.letters + other.letters)

    def __eq__(self, other):
        return self.letters == other.letters
    
    def __lt__(self, other):
        return (len(self.letters),tuple(self.letters)) < (len(other.letters),tuple(other.letters))

    def __hash__(self):
        return hash(tuple(self.letters))
    
    def inverse(self):
        return A2Word([-x for x in reversed(self.letters)])

################

geometric_gens = [A2Word([i]) for i in [-7,-6,-5,-4,-3,-2,-1,1,2,3,4,5,6,7]]

def geometric_ball(radius):
    return grp.cayley_ball(geometric_gens, radius)