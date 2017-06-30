class Pair:
    PTYPE_MUST = 1
    PTYPE_NO = 2
    PTYPE_CAN = 3
    
    PTYPE_STR = {1:"MUST", 2:"NO", 3:"CAN"}

    def __init__(self, id1, id2, ptype ):
        self.id1 = id1
        self.id2 = id2
        self.ptype = ptype