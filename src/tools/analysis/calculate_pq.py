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
