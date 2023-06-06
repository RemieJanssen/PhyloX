from enum import Enum

class NetworkClass(str, Enum):
    TC = "tree-child"
    TB = "tree-based"
    OR = "orchard"
    SF = "stack-free"
    BI = "binary"