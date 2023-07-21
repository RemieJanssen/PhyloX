from enum import Enum


class DiNetworkClass(Enum):
    TC = "tree-child"
    TB = "tree-based"
    OR = "orchard"
    SF = "stack-free"
    BI = "binary"
