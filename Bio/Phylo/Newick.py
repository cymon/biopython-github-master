# Copyright (C) 2009 by Eric Talevich (eric.talevich@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

"""Classes corresponding to Newick trees, also used for Nexus trees.

See classes in Bio.Nexus: Trees.Tree, Trees.NodeData, and Nodes.Chain.
"""
__docformat__ = "epytext en"

import warnings

import BaseTree


class Tree(BaseTree.Tree):
    """Newick Tree object."""

    def __init__(self, root=None, rooted=False, id=None, name='', weight=1.0):
        BaseTree.Tree.__init__(self, root=root or Clade(),
                rooted=rooted, id=id, name=name)
        self.weight = weight


class Clade(BaseTree.Subtree):
    """Newick Clade (subtree) object."""

    def __init__(self, branch_length=1.0, name=None, clades=None,
            support=None, comment=None):
        BaseTree.Subtree.__init__(self, branch_length=branch_length,
                name=name, clades=clades)
        self.support = support
        self.comment = comment

