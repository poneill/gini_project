"""
Assume a file is of the form:

tail -n5 */*.csv > filename

Parse the resulting file and return the energy matrices.
"""

from utils import *

def chunk2matrix(chunk):
        relevant_chunk = chunk[1:]
        rows = [map(float,line.split(","))
                for line in relevant_chunk
                if not line == "\n"]
        matrix = transpose(rows)
        assert all(len(col) == 4 for col in matrix),chunk
        return matrix
        

def parse_energy_matrices(filename):
    with open(filename) as f:
        lines = f.readlines()
    matrix_chunks = split_on(lines,lambda line:line.startswith("="))
    matrices = map(chunk2matrix,matrix_chunks)
    return matrices
    
