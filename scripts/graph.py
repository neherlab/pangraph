from numpy import np

from block import Block

class Graph(object):
	"""docstring for Graph"""
	def __init__(self, sequence_name, sequence):
		super(Graph, self).__init__()
		tmp = Block()
		self.blocks = [tmp.from_sequence(sequence_name, sequence)]
		self.sequences = {sequence_name: self.blocks[0]}


