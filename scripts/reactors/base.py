
class Reactor(object):
    """
    Base class that performs a reaction on sim object
    """
    def __init__(self, conversion=1):
        self.conversion = conversion

    def react(self, sim):
        """
        Perform reaction on sim, return with changes
        """
        raise NotImplementedError
