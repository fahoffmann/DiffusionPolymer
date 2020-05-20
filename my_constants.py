from my_geometry import total_number_of_segments

class Constants:

    def __init__(self):
        self.number_of_iterations = 8
        self.dimension = 3
        self._radius = 1
        self.concentration = 0.2
        self.D_r = 1
        self.MC_steps = 10000000
        self.min_walk_length = 1
        self.max_walk_length = 2
        self._volume = None
        self._number_of_walkers = None
        self.diffusion = True

    @property
    def volume(self):
        return self._volume

    @property
    def number_of_walkers(self):
        return self._number_of_walkers

    @property
    def radius(self):
        return self._radius

    @volume.setter
    def volume(self, value):
        self._number_of_walkers = int(self.concentration * value)
        self._volume = value

    @number_of_walkers.setter
    def number_of_walkers(self, value):
        self._number_of_walkers = value

    @radius.setter
    def radius(self, value):
        self._volume = total_number_of_segments(self.dimension, value)
        self._number_of_walkers = int(self.concentration * self.volume)
        self._radius = value