class Scenarios(dict):
    """
    Scenarios is the class that stores and define multiple scenarios that will be simulated through the renewable
    simulator. It is inherited from the python built-in class, dict.

    The dictionary values only accept tuples or lists.

    A single scenario will be defined by a certain combination of values from all keys. For example, scenario #1 is
    defined using the first value in all keys; scenario #2 is defined using the second value in the first key and
    the first value in all other keys.
    """

    def total_scenarios(self):
        """
        Counts the total number of scenarios that could be defined using the current object of the class Scenarios.
        The total number is calculated from the multiplication of the length of each dictionary entry.

        :return: the number of total scenarios
        """
        total = 1;

        for key in self:
            if isinstance(self[key], (list, tuple)):
                total *= len(self[key])
            else:
                raise Exception('Scenarios only allow lists or tuples as values')
        return total

    def get_scenario(self, index):

        import math

        # Check out-of-bound issues
        total_scenarios = self.total_scenarios()
        if index >= total_scenarios:
            raise Exception("Index {} is out of bound. The toal number of scenarios is {}".format(index, total_scenarios))

        # Initialize variables
        accumulate_len = 1
        current_scenario = {}

        for key, value in self.items():

            # Calculate the current index for this configuration
            if accumulate_len == 1:
                accumulate_len = len(value)
                current_index = index % accumulate_len
            else:
                current_index = math.floor(index / accumulate_len) % len(value)
                accumulate_len *= len(value)

            current_scenario[key] = value[current_index]

        return current_scenario

    def print(self):
        from pprint import pprint
        print("There are {} scenarios defined by the current dictionary.".format(self.total_scenarios()))
        pprint(self)
