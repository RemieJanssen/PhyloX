from phylox.exceptions import InvalidMoveDefinition
from phylox.rearrangement.movetype import MoveType


class Move(object):
    def __init__(self, *args, **kwargs):
        try:
            self.move_type = kwargs["move_type"]
        except KeyError:
            raise InvalidMoveDefinition("Missing move_type.")

        # None type move
        if self.move_type == MoveType.NONE:
            return

        # TAIL/HEAD move (i.e. RSPR/horizontal)
        if self.move_type == MoveType.RSPR:
            raise InvalidMoveDefinition(
                "rSPR moves must be defined as moves of type tail or head."
            )
        if self.move_type in [MoveType.TAIL, MoveType.HEAD]:
            try:
                self.origin = kwargs["origin"]
                self.moving_edge = kwargs["moving_edge"]
                self.target = kwargs["target"]
            except KeyError:
                raise InvalidMoveDefinition(
                    "Missing one of origin, moving_edge, or target."
                )

            if self.move_type == MoveType.TAIL:
                self.moving_node = self.moving_edge[0]
            else:
                self.moving_node = self.moving_edge[1]

            if self.moving_edge == self.target:
                raise InvalidMoveDefinition(
                    "Moving edge must not be the target edge."
                )

            return

        # TODO Write vertical move parsing
        if self.move_type == MoveType.VPLU:
            try:
                self.start_edge = kwargs["start_edge"]
                self.end_edge = kwargs["end_edge"]
                self.start_node = kwargs["start_node"]
                self.end_node = kwargs["end_node"]
            except KeyError:
                raise InvalidMoveDefinition(
                    "Missing one of start_edge, end_edge, start_node, or end_node."
                )
            return

    def is_type(self, move_type):
        if (
            self.move_type == MoveType.NONE
            or (
                move_type == MoveType.RSPR
                and self.move_type in [MoveType.TAIL, MoveType.HEAD]
            )
            or (
                move_type == MoveType.VERT
                and self.move_type in [MoveType.VPLU, MoveType.VMIN]
            )
        ):
            return True
        return move_type == self.move_type