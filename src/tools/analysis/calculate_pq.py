import numpy as np

def calculate_pq(positions, q_values):
    """
    Calculate the structure factor P(q) for a series of frames.

    Parameters:
    positions: list of np.array, each of shape (N, 3) where N is number of particles
    q_values: np.array of q values

    Returns:
    np.array of shape (num_frames, num_q) with P(q) for each frame and q
    """
    pq_all = []
    for frame in positions:
        N = len(frame)
        Sq = np.ones(len(q_values)) / N  # initialize with 1/N for self-overlap
        for i in range(N - 1):
            for j in range(i + 1, N):
                r = np.linalg.norm(frame[i] - frame[j])
                if r > 0:  # avoid division by zero, though unlikely
                    for k, q in enumerate(q_values):
                        Sq[k] += 2.0 / N / N * np.sin(q * r) / (q * r)
        pq_all.append(Sq)
    return np.array(pq_all)


def calculate_gr(positions, r_bins, box_length):
    """
    Calculate the radial distribution function g(r) for a series of frames.

    Parameters:
    positions: list of np.array, each of shape (N, 3) where N is number of particles
    r_bins: np.array of r bin edges
    box_length: float, length of the simulation box for periodic boundaries

    Returns:
    np.array of shape (num_frames, len(r_bins)-1) with g(r) for each frame and r bin
    """
    gr_all = []
    for frame in positions:
        N = len(frame)
        rho = N / (box_length ** 3)
        distances = []
        for i in range(N - 1):
            for j in range(i + 1, N):
                dr = frame[i] - frame[j]
                # Apply minimum image convention
                dr = dr - box_length * np.round(dr / box_length)
                r = np.linalg.norm(dr)
                distances.append(r)
        hist, _ = np.histogram(distances, bins=r_bins)
        r_mid = (r_bins[:-1] + r_bins[1:]) / 2
        dr = r_bins[1] - r_bins[0]  # assuming uniform bins
        expected = 4 * np.pi * r_mid**2 * dr * rho * (N / 2)
        gr = hist / expected
        gr_all.append(gr)
    return np.array(gr_all)
