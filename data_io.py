import numpy as np

def read_npy_file(file_path):
    """Reads a .npy file and returns the data."""
    data = np.load(file_path)
    return data

def read_dat_file(file_path):
    """Reads a .dat file and returns the data as a list."""
    with open(file_path, 'r') as file:
        data = file.readlines()
    return [line.strip() for line in data]