from phylox.rearrangement.move import Move
from phylox.rearrangement.movetype import MoveType

def apply_move(network, move):
    """
    Apply a move to the network.
    returns True if successful, and False otherwise.
    """
    check_valid(network, move)

    if move.move_type in [MoveType.TAIL, MoveType.HEAD]:
        if move.moving_node in move.target:
            # move does not impact the network
            return True
        network.remove_edges_from(
            [
                (move.origin[0], move.moving_node),
                (move.moving_node, move.origin[1]),
                move.target,
            ]
        )
        network.add_edges_from(
            [
                (move.target[0], move.moving_node),
                (move.moving_node, move.target[1]),
                move.origin,
            ]
        )
        return
    if move.move_type in [MoveType.NONE]:
        return
    # TODO implement vertical moves
    raise InvalidMove("only tail or head moves are currently valid.")

def apply_move_sequence(network, seq_moves):
    for move in seq_moves:
        network.apply_move(move)