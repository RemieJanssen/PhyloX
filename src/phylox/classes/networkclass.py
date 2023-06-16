from enum import Enum

class DiNetworkClass(str, Enum):
    TC = "tree-child"
    TB = "tree-based"
    OR = "orchard"
    SF = "stack-free"
    BI = "binary"