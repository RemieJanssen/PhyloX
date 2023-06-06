from enum import Enum


class MoveType(str, Enum):
    NONE = "NONE"
    TAIL = "TAIL"
    HEAD = "HEAD"
    RSPR = "RSPR"
    VPLU = "VPLU"  # not currently in use
    VMIN = "VMIN"  # not currently in use
    VERT = "VERT"