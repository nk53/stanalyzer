from math import ceil


def chunk_frames(n_frames: int, n_workers: int):
    """
    Split frame indices into roughly equal chunks.

    Example:
        chunk_frames(1001, 4)

        ->

        [
            (0, 251),
            (251, 501),
            (501, 751),
            (751, 1001),
        ]
    """

    if n_frames < 0:
        raise ValueError("n_frames must be non-negative")

    if n_workers < 1:
        raise ValueError("n_workers must be at least 1")

    if n_frames == 0:
        return []

    chunk_size = ceil(n_frames / n_workers)

    chunks = []

    for start in range(0, n_frames, chunk_size):
        stop = min(start + chunk_size, n_frames)
        chunks.append((start, stop))

    return chunks
