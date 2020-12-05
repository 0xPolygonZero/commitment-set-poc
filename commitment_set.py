#!/usr/bin/env python3

import abc
import collections
import heapq
import random
import sys

class LivenessMask:
    """Tracks which commitments are live. Essentially a sparse bit vector."""

    def __init__(self, mask):
        self.mask = mask

    def random(total_commitments, live_commitments):
        consumed_commitments = total_commitments - live_commitments
        mask = [1] * live_commitments + [0] * consumed_commitments
        random.shuffle(mask)
        return LivenessMask(mask)

    def to_rle(self):
        return [len(zero_run) for zero_run in repr(self).split('1')]

    def from_rle(rle):
        bit_string = '1'.join(['0' * zero_run for zero_run in rle])
        return LivenessMask.from_bit_string(bit_string)

    def compress(self):
        return huffman_compress(self.to_rle())

    def decompress(bits):
        return LivenessMask.from_rle(huffman_decompress(bits))

    def __repr__(self):
        """Represents this mask as a string of zeros and ones."""
        return ''.join(str(int(bit)) for bit in self.mask)

    def from_bit_string(bit_string):
        return LivenessMask([bool(int(bit)) for bit in bit_string])

def huffman_compress(symbols):
    code_tree = HuffmanCodeTree.for_symbols(symbols)
    symbol_prefix_map = code_tree.symbol_prefix_map()
    bits = code_tree.serialize()
    for symbol in symbols:
        bits += symbol_prefix_map[symbol]
    return bits

def huffman_decompress(bits):
    code_tree = HuffmanCodeTree.parse(bits)
    symbols = []
    while bits:
        node = code_tree
        while isinstance(node, HuffmanCodeNonLeaf):
            branch = bits.pop(0)
            node = [node.left, node.right][branch]
        symbols.append(node.symbol)
    return symbols

"""A tree (not necessarily the root) of Huffman codes."""
class HuffmanCodeTree(abc.ABC):
    def __lt__(self, other):
        return self.count < other.count

    def for_symbols(symbols):
        """Creates a Huffman code tree for the given sequence of symbols."""

        # Make a leaf for each symbol, and organize them as a min-heap
        counts = collections.Counter(symbols).items()
        trees = [HuffmanCodeLeaf(count, symbol) for symbol, count in counts]
        heapq.heapify(trees)

        # Repeatedly merge the two lowest-weight subtrees until we're left with a single tree.
        while len(trees) > 1:
            left = heapq.heappop(trees)
            right = heapq.heappop(trees)
            heapq.heappush(trees, HuffmanCodeNonLeaf(left, right))

        return trees[0]

    def symbol_prefix_map(self):
        """Generates a map from symbols to their prefixes in this tree."""
        result = {}
        self.symbol_prefix_map_helper([], result)
        return result

    @abc.abstractmethod
    def symbol_prefix_map_helper(self, prefix, result):
        pass

    @abc.abstractmethod
    def serialize(self):
        pass

    def parse(bits):
        leaf = bits.pop(0)
        if leaf:
            return HuffmanCodeLeaf.parse_leaf(bits)
        else:
            return HuffmanCodeNonLeaf.parse_non_leaf(bits)

def int_to_bit_list(n):
    bits = []
    while n:
        bits.append(bool(n & 1))
        n >>= 1
    return bits

def bit_list_to_int(bits):
    return sum(2**i * b for i, b in enumerate(bits))

"""A leaf node in a Huffman code tree."""
class HuffmanCodeLeaf(HuffmanCodeTree):
    def __init__(self, count, symbol):
        self.count = count
        self.symbol = symbol

    def symbol_prefix_map_helper(self, prefix, result):
        result[self.symbol] = prefix

    def serialize(self):
        bits = [True] # is-leaf flag

        # To encode the symbol length (in bits) n, we will append n zeros followed by a one.
        symbol_bits = int_to_bit_list(self.symbol)
        for bit in symbol_bits:
            bits.append(False)
        bits.append(True)

        # Append the symbol itself.
        bits += symbol_bits

        return bits

    def parse_leaf(bits):
        # Parse the symbol length n, encoded as n zeros followed by a one.
        symbol_len = 0
        while not bits.pop(0):
            symbol_len += 1

        # Parse the symbol.
        symbol = bit_list_to_int(bits[:symbol_len])
        del bits[:symbol_len]

        return HuffmanCodeLeaf(None, symbol)

    def __repr__(self):
        return "{}x{}".format(self.symbol, self.count)

"""A non-leaf node in a Huffman code tree."""
class HuffmanCodeNonLeaf(HuffmanCodeTree):
    def __init__(self, left, right):
        if left.count is not None and right.count is not None:
            self.count = left.count + right.count
        else:
            self.count = None

        self.left = left
        self.right = right

    def symbol_prefix_map_helper(self, prefix, result):
        self.left.symbol_prefix_map_helper(prefix + [False], result)
        self.right.symbol_prefix_map_helper(prefix + [True], result)

    def serialize(self):
        bits = [False] # is-leaf flag
        bits += self.left.serialize()
        bits += self.right.serialize()
        return bits

    def parse_non_leaf(bits):
        left = HuffmanCodeTree.parse(bits)
        right = HuffmanCodeTree.parse(bits)
        return HuffmanCodeNonLeaf(left, right)

    def __repr__(self):
        return "({}, {})".format(self.left, self.right)

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print('Usage: {} [total commitments] [live commitments]'.format(sys.argv[0]))
        exit()

    total_commitments = int(sys.argv[1])
    live_commitments = int(sys.argv[2])

    lm = LivenessMask.random(total_commitments, live_commitments)
    bits = lm.compress()
    lm_parsed = LivenessMask.decompress(bits[:])

    print('Total commitments: {}'.format(total_commitments))
    print('Live commitments: {}'.format(live_commitments))
    print('Fraction live: {}%'.format(100 * live_commitments / total_commitments))
    print('Compressed size: {} bits'.format(len(bits)))
    print('Bits per live commitment: {}'.format(len(bits) / live_commitments))
    print('Serialization test: {}'.format('pass' if lm.mask == lm_parsed.mask else 'fail'))
