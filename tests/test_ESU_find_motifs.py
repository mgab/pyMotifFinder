#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Test suite for the BruteForce module."""

import unittest as ut
import random

import networkx as nx

from ESU import *

class EnumerateSubgraphsTests(ut.TestCase):
    def setUp(self):
        self.tree_graph = nx.balanced_tree(2, 2, nx.DiGraph())

    def test_tree_graph(self):
        self.assertEqual(list(enumerate_subgraphs(self.tree_graph, size=3)),
                         [[0, 1, 2], [0, 1, 3], [0, 1, 4], [0, 2, 5],
                          [0, 2, 6], [1, 3, 4], [2, 5, 6]])
        self.assertEqual(list(enumerate_subgraphs(self.tree_graph, size=4)),
                         [[0, 1, 2, 3], [0, 1, 2, 4], [0, 1, 2, 5],
                          [0, 1, 2, 6], [0, 1, 3, 4], [0, 2, 5, 6]])

    def test_null_subgraoh(self):
        self.assertEqual(list(enumerate_subgraphs(self.tree_graph, size=0)),
                         [[]])


class RandomizeGraphTests(ut.TestCase):
    def setUp(self):
        self.tree_graph = nx.balanced_tree(2, 2, nx.DiGraph())

    def test_tree_graph(self):
        res = randomize_graph(self.tree_graph, prng=random.Random(1))
        self.assertFalse(nx.is_isomorphic(self.tree_graph, res))
        self.assertEqual(self.tree_graph.in_degree(), res.in_degree())
        self.assertEqual(self.tree_graph.out_degree(), res.out_degree())


if __name__ == '__main__':
    ut.main()
