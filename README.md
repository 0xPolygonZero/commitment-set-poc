# Mir Commitment Set Proof of Concept

This is a proof-of-concept implementation of Mir's [Commitment Set](https://mirprotocol.org/blog/Reducing-state-size-on-Mir) storage model, which is based on the [Modified Huffman coding](https://en.wikipedia.org/wiki/Modified_Huffman_coding). The focus is on demonstrating space efficiency, so this code is not very performant and should not be used in production.


## Usage

```
python3 commitment_set.py [total commitments] [live commitments]
```


## Example

```
$ python3 commitment_set.py 1000000 100000
Total commitments: 1000000
Live commitments: 100000
Fraction live: 10.0%
Compressed size: 473763 bits
Bits per live commitment: 4.73763
Serialization test: pass
```
