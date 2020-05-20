from numpy import sqrt

def total_number_of_segments(D, R):
    if (D == 1):
        n = 1
        i = 1
        while (i <= R):
            n += 2
            i += 1
        return n
    n = total_number_of_segments(D - 1, R)
    i = 1
    while (i <= R):
        round_correction = 0.0000005
        n += 2 * total_number_of_segments(D - 1, sqrt(R ** 2 - i ** 2) - round_correction)
        i += 1
    return n