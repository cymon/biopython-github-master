# Copyright (C) 2010 by Eric Talevich (eric.talevich@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

"""Short helper functions (syntax sugar) used in Bio.Phylo.

The amount of code in this file should be kept to a minimum.
"""
__docformat__ = "epytext en"


def trim_str(text, maxlen=60):
    """Truncate a string to maxlen characters, including ellipsis."""
    assert isinstance(text, basestring), \
            "%s should be a string, not a %s" % (text, type(text))
    if len(text) > maxlen:
        return text[:maxlen-3] + '...'
    return text
