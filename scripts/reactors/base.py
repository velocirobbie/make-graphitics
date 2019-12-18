
class Reactor(object):
    """
    Base class that performs a reaction on sim object
    """

    def react(self, sim):
        """
        Perform reaction on sim, return with changes
        """
        raise NotImplementedError
