from enum import Enum


class MoveType(str, Enum):
    NONE = "NONE"
    TAIL = "TAIL"
    HEAD = "HEAD"
    RSPR = "RSPR"
    VPLU = "VPLU"
    VMIN = "VMIN"
    VERT = "VERT"
    ALL = "ALL"
