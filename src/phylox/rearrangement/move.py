from phylox.exceptions import InvalidMoveDefinitionException, InvalidMoveException
from phylox.rearrangement.movetype import MoveType
from copy import deepcopy
from phylox.rearrangement.movability import check_valid
from phylox.base import find_unused_node

def apply_move(network, move):
    """
    Apply a move to the network, not in place.
    returns True if successful, and False otherwise.
    """
    check_valid(network, move)
    new_network = deepcopy(network)

    if move.move_type in [MoveType.TAIL, MoveType.HEAD]:
        if move.moving_node in move.target:
            # move does not impact the network
            return True
        new_network.remove_edges_from(
            [
                (move.origin[0], move.moving_node),
                (move.moving_node, move.origin[1]),
                move.target,
            ]
        )
        new_network.add_edges_from(
            [
                (move.target[0], move.moving_node),
                (move.moving_node, move.target[1]),
                move.origin,
            ]
        )
        return new_network
    elif move.move_type in [MoveType.VPLU]:
        new_network.remove_edges_from([move.start_edge, move.end_edge])        
        new_network.add_edges_from([(move.start_edge[0], move.start_node), (move.start_node, move.start_edge[1]),(move.end_edge[0], move.end_node), (move.end_node, move.end_edge[1]), (move.start_node,move.end_node)])
        return new_network
    elif move.move_type in [MoveType.VMIN]:
        parent_0 = network.parent(move.removed_edge[0], exclude=[move.removed_edge[1]])
        child_0 = network.child(move.removed_edge[0], exclude=[move.removed_edge[1]])
        parent_1 = network.parent(move.removed_edge[1], exclude=[move.removed_edge[0]])
        child_1 = network.child(move.removed_edge[1], exclude=[move.removed_edge[0]])
        new_network.remove_edges_from([(parent_0, move.removed_edge[0]), (move.removed_edge[0], child_0),(parent_1, move.removed_edge[1]), (move.removed_edge[1], child_1), move.removed_edge])
        new_network.add_edges_from([(parent_0, child_0), (parent_1, child_1)])
        return new_network
    elif move.move_type in [MoveType.NONE]:
        return network
    raise InvalidMove("only tail or head moves are currently valid.")



def apply_move_sequence(network, seq_moves):
    for move in seq_moves:
        network = apply_move(network, move)
    return network


class Move(object):
    def __init__(self, *args, **kwargs):
        try:
            self.move_type = kwargs["move_type"]
        except KeyError:
            raise InvalidMoveDefinitionException("Missing move_type.")

        # None type move
        if self.move_type == MoveType.NONE:
            return

        # TAIL/HEAD move (i.e. RSPR/horizontal)
        elif self.move_type == MoveType.RSPR:
            raise InvalidMoveDefinitionException(
                "rSPR moves must be defined as moves of type tail or head."
            )
        elif self.move_type in [MoveType.TAIL, MoveType.HEAD]:
            try:
                self.origin = kwargs["origin"]
                self.moving_edge = kwargs["moving_edge"]
                self.target = kwargs["target"]
            except KeyError:
                raise InvalidMoveDefinitionException(
                    "Missing one of origin, moving_edge, or target."
                )

            if self.move_type == MoveType.TAIL:
                self.moving_node = self.moving_edge[0]
            else:
                self.moving_node = self.moving_edge[1]

            if self.moving_edge == self.target:
                raise InvalidMoveDefinitionException(
                    "Moving edge must not be the target edge."
                )

            return
        # VERT move (i.e. SPR/vertical)
        elif self.move_type == MoveType.VERT:
            raise InvalidMoveDefinitionException(
                "vertical moves must be defined as moves of type VPLU or VMIN."
            )
        # VPLU/VMIN move
        elif self.move_type == MoveType.VPLU:
            try:
                self.start_edge = kwargs["start_edge"]
                self.end_edge = kwargs["end_edge"]
                self.start_node = kwargs.get("start_node", None)
                self.end_node = kwargs.get("end_node", None)
                network = kwargs.get("network", None)
                if (self.start_node is None or self.end_node is None) and network is None:
                    raise InvalidMoveDefinitionException(
                        "Either a start_node and end_node, or a network must be given."
                    )
                if self.start_node is None:
                    self.start_node = find_unused_node(network)
                if self.end_node is None:
                    self.end_node = find_unused_node(network, exclude=[self.start_node])

                if self.start_edge == self.end_edge:
                    raise InvalidMoveDefinitionException(
                        "Start edge must not be the end edge."
                    )
            except KeyError:
                raise InvalidMoveDefinitionException(
                    "Missing one of start_edge, end_edge."
                )
        elif self.move_type == MoveType.VMIN:
            try:
                self.edge = kwargs["removed_edge"]
            except KeyError:
                raise InvalidMoveDefinitionException(
                    "Missing removed_edge in definition."
                )
        else:
            raise InvalidMoveDefinitionException("Invalid move type.")



    def is_type(self, move_type):
        if self.move_type == MoveType.ALL:
            return True
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