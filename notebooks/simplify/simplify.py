from .utils import Pangraph
from .circularize import remove_transitive_edges


def simplify(graph: Pangraph, focal_paths: list[str]):
    """Remove all paths that are not in focal_paths, and simplify the graph removing all transitive edges."""

    # remove unwanted paths
    path_ids_to_remove = [
        pid for pid in graph.paths.keys() if graph.paths[pid].name not in focal_paths
    ]

    for pid in path_ids_to_remove:
        graph.remove_path(pid)

    # simplify graph
    remove_transitive_edges(graph)
