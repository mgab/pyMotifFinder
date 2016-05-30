# -*- coding: utf-8 -*-
"""Find network motifs using the Enumerate SUbgraph (ESU) algorithm."""

from __future__ import division, print_function

from operator import itemgetter
import random

import networkx as nx
import networkx.algorithms.isomorphism as iso


def find_pattern(graph, pattern, sign_sensitive=False):
    """Find all the subgraphs isomorphic to the given pattern.

    Return an iterator over all the subgraphs in `graph` that are isomorphic to
    the specified `pattern`. It works for both directed and undirected graphs,
    but the type of `graph` and `pattern` arguments must be coherent.
    If `sign_sensitive` is False, only the topology is considered for the
    isomorphism test. Otherwise, the "sign" attribute of each edge must have
    the same value.

    Based on the code from
    https://zulko.wordpress.com/2012/10/13/finding-a-subnetwork-with-a-given-topology-in-e-coli/
    """
    edge_match = None
    if sign_sensitive:
        def edge_match(e1, e2):
            return (e1["sign"] == e2["sign"])

    if graph.is_directed() and pattern.is_directed():
        matcher = iso.DiGraphMatcher(graph, pattern, edge_match=edge_match)
    elif not graph.is_directed() and not pattern.is_directed():
        matcher = iso.GraphMatcher(graph, pattern, edge_match=edge_match)
    else:
        raise TypeError("type of `graph` and `pattern` arguments is not "
                        "coherent!")

    return matcher.subgraph_isomorphisms_iter()


def enumerate_subgraphs(graph, seed_nodes=None, size=3):
    """Find all connected subgraphs of the given size.

    Return an iterator over all the connected subgraphs of the specified
    `size`. If `seed_nodes` is not None, it lists only subgraphs containing
    those nodes. Note that if the `seed_nodes` are not connected between them,
    disconnected subgraphs may be returned.

    It uses the algorithm EnumerateSUbgraphs or ESU used by FANMOD and
    described in:
    S. Wernicke, "Efficient Detection of Network Motifs," IEEE/ACM Transactions
    on Computational Biology and Bioinformatics, 2006. doi:10.1109/TCBB.2006.51
    """
    if seed_nodes is None:
        seed_nodes = []

    # Change the node names for numerical labels (required by ESU)
    node_mapping = dict(zip(graph.nodes(), range(graph.number_of_nodes())))
    graph = nx.relabel_nodes(graph, node_mapping, copy=True)
    seed_nodes = [node_mapping[n] for n in seed_nodes]

    # Use a recurrent function to get each subgraph and change the names back
    rev_node_mapping = {v: k for k, v in node_mapping.iteritems()}
    del node_mapping
    for subgraph in _enumerate_subgraphs(graph, seed_nodes=seed_nodes,
                                         size=size):
        yield [rev_node_mapping[n] for n in subgraph]


def _enumerate_subgraphs(graph, seed_nodes, size):
    """Private method to Enumerate SUbgraphs using a recursive approach."""
    if len(seed_nodes) == size:
        yield seed_nodes
    else:
        if len(seed_nodes) == 0:
            neighbors = xrange(graph.order())
        else:
            max_node = max(seed_nodes)
            neighbors = (ngbr for n in seed_nodes
                         for ngbr in graph.predecessors(n)+graph.successors(n)
                         if max_node < ngbr and ngbr not in seed_nodes)

        for ngbr in neighbors:
            for result in _enumerate_subgraphs(graph, seed_nodes+[ngbr], size):
                yield result


def count_unique_topologies(topologies):
    """Count the number of instances of each unique topology.

    Accept a list of graph objects. Return a list containing one list for each
    isomorphic topology found, with the graph object as the first element and
    the number of found instances as the second.
    """
    unique_topologies = []
    for topo in topologies:
        for i in range(len(unique_topologies)):
            if nx.is_isomorphic(unique_topologies[i][0], topo):
                unique_topologies[i][1] += 1
                break
        else:
            unique_topologies.append([topo, 1])

    return unique_topologies


def randomize_graph(graph, swap_steps=None, prng=None, maxsteps=None):
    u"""Randomize a graph preserving the in and out degree of each node.

    Return a randomized copy a graph by crossing pairs of edges. It randomly
    selects two edges, removes them, and connects the source nodes of the first
    and second edges with the target nodes of the second and first edge,
    respectively.

    Parameters
    ----------
    swap_steps : int, optional (default=3*graph.number_of_edges())
        Number of edge swaps before considering the network is successfully
        randomized.
    maxsteps : int, optional (default=swap_steps*10)
        Some edge swaps are not allowed (because one of the new edges already
        exists). This indicates the total number of edge swap trials before
        raising a `RuntimeError`.

    This is the randomization algorithm used in mfinder and described in:
    Yeger-Lotem E, Sattath S, Kashtan N, et al. "Network motifs in integrated
    cellular networks of transcription–regulation and protein–protein
    interaction." Proc. Natl. Acad. Sci. USA, 2004. doi:10.1073/pnas.0306752101
    """
    graph = graph.copy()
    prng = prng or random.Random()
    swap_steps = swap_steps or graph.number_of_edges()*3
    maxsteps = maxsteps or swap_steps * 10

    step = swaps = 0
    while swaps < swap_steps:
        step += 1
        if step > maxsteps:
            raise RuntimeError("Reached max number of steps in the "
                               "randomization process.")
        (s1, t1), (s2, t2) = prng.sample(graph.edges(), 2)
        if not graph.has_edge(s1, t2) and not graph.has_edge(s2, t1):
            swaps += 1
            graph.add_edges_from(((s1, t2), (s2, t1)))
            graph.remove_edges_from(((s1, t1), (s2, t2)))

    return graph


def find_motifs_slow(graph, size=3, min_occurrences=5,  rand_networks=1000,
                     prng=None):
    u"""Count all motifs of a given size and its statistical relevance.

    Identify all motifs of a given size in a network asses if they are enriched
    in a statistically significant way (p-value) by comparing with
    randomizations of that network.

    It will return a list with one tuple for each of the found motifs. Each of
    the tuples will contain a graph object with the motif topology, the number
    of occurrences in the original network, and the probability of finding that
    motive at least that number of times in randomized networks (p-value of
    enrichment).

    WARNING: A much faster version of this function is `find_motifs`. This one
        looks for every motif one by one in every randomized network, while the
        other uses the ESU algorithm to list all motifs also for the randomized
        networks.

    Parameters
    ----------
    size : int, optional (default=3)
        Size of the motifs to search.
    min_occurrences : int, optional (default=5)
        Indicates the minimum number of times a motif has to be found in the
        original network to be considered.
    rand_networks : int, optional (default=1000)
        Number of randomized networks used to compute the statistical
        significance of the enrichment.
    prng : int, optional
        Pseudo Random Number Generator. If not provided is taken from the
        `random` module.
    """
    all_subgraphs_iter = enumerate_subgraphs(graph, size=size)
    motifs = count_unique_topologies((graph.subgraph(sub)
                                      for sub in all_subgraphs_iter))
    motifs = [[topo, n, 0] for topo, n in motifs if n >= min_occurrences]

    # Count how many times we find the same motifs in random networks
    for i in range(rand_networks):
        if i % 10 == 0:
            print(i, end=" ")
        rgraph = randomize_graph(graph, prng=prng)
        for i in range(len(motifs)):
            rmotifs = find_pattern(rgraph, motifs[i][0])
            # Avoid counting permutations of the same graph
            if motifs[i][1] <= len(set((tuple(sorted(k)) for k in rmotifs))):
                motifs[i][2] += 1

    return sorted([(topo, n, float(m)/rand_networks) for topo, n, m in motifs],
                  key=itemgetter(2))


def find_motifs(graph, size=3, min_occurrences=5,  rand_networks=1000,
                ping_every=0, prng=None):
    u"""Count all motifs of a given size and its statistical relevance.

    Identify all motifs of a given size in a network asses if they are enriched
    in a statistically significant way (p-value) by comparing with
    randomizations of that network.

    It will return a list with one tuple for each of the found motifs. Each of
    the tuples will contain a graph object with the motif topology, the number
    of occurrences in the original network, and the probability of finding that
    motive at least that number of times in randomized networks (p-value of
    enrichment).

    Parameters
    ----------
    size : int, optional (default=3)
        Size of the motifs to search.
    min_occurrences : int, optional (default=5)
        Indicates the minimum number of times a motif has to be found in the
        original network to be considered.
    rand_networks : int, optional (default=1000)
        Number of randomized networks used to compute the statistical
        significance of the enrichment.
    prng : int, optional
        Pseudo Random Number Generator. If not provided is taken from the
        `random` module.
    ping_every : int, optional (default=0)
        Print the number of randomized networks processed every `ping_every`
        networks. If equal to 0 or False nothing is printed.
    """
    all_subgraphs_iter = enumerate_subgraphs(graph, size=size)
    motifs = count_unique_topologies((graph.subgraph(sub)
                                      for sub in all_subgraphs_iter))
    motifs = [[topo, n, 0] for topo, n in motifs if n >= min_occurrences]

    # Count how many times we find the same motifs in random networks
    for i in range(rand_networks):
        if ping_every and i % ping_every == 0:
            print(i, end=" ")

        rgraph = randomize_graph(graph, prng=prng)
        rand_subgraphs_iter = enumerate_subgraphs(rgraph, size=size)
        rand_motifs = count_unique_topologies((rgraph.subgraph(sub)
                                               for sub in rand_subgraphs_iter))
        for rmotif, rn in rand_motifs:
            if rn < min_occurrences:
                continue
            for i in range(len(motifs)):
                if rn >= motifs[i][1] and nx.is_isomorphic(motifs[i][0],
                                                           rmotif):
                    motifs[i][2] += 1
    if ping_every:
        print()  # Line break in the output
    return sorted([(topo, n, float(m)/rand_networks) for topo, n, m in motifs],
                  key=itemgetter(2))
