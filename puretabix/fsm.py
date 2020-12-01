"""
Finite State Machine
--------------------

Simple implementation of a finite state machine, used for parsing in a
character-by-character manner


based on https://gist.github.com/brianray/8d3e697dbbf150f725291d74ac0cee8b
"""

import re


class RegexTransition:
    __slots__ = ["dst", "condition", "callback"]

    def __init__(self, destination_state, condition, callback):
        self.dst = destination_state
        self.condition = re.compile(condition)
        self.callback = callback

    def match(self, _input):
        return self.condition.match(_input)


class SetInTransition:
    __slots__ = ["dst", "condition", "callback"]

    def __init__(self, destination_state, condition, callback):
        self.dst = destination_state
        self.condition = frozenset(condition)
        self.callback = callback

    def match(self, _input):
        return _input in self.condition


class SetNotInTransition:
    __slots__ = ["dst", "condition", "callback"]

    def __init__(self, destination_state, condition, callback):
        self.dst = destination_state
        self.condition = frozenset(condition)
        self.callback = callback

    def match(self, _input):
        return _input not in self.condition


class FSMachine:
    def __init__(self):
        self.transitions = {}

    def add_transition(
        self, start_state, end_state, transition_class, condition, callback
    ):
        """

        callback should be a function that accepts the current character of the state
        machine as the last positional arguments i.e. it is called like:
        callback(*args, current_char, **kwargs)

        because callback uses the last positional argument as the character, the first
        positional argment can be used by a class instance as self
        """
        if start_state not in self.transitions:
            self.transitions[start_state] = []
        self.transitions[start_state].append(
            transition_class(end_state, condition, callback)
        )

    def run(self, inputs, initial_state, *args, **kwargs):
        self.current_state = initial_state
        for c in inputs:
            self.process_next(c, args, kwargs)
            # if state is None, early exist from inputs
            if not self.current_state:
                break

        assert not self.current_state, f"Unexpected ending in {self.current_state}"

    def process_next(self, _input, callback_args, callback_kwargs):
        frozen_state = self.current_state
        for transition in self.transitions[frozen_state]:
            if transition.match(_input):
                # found a transition that matches
                # update the state
                self.current_state = transition.dst
                # call the callback, if it exists
                # because callback uses the last positional argument as the input, the first
                # positional argment can be used by a class instance as self
                if transition.callback:
                    transition.callback(*callback_args, _input, **callback_kwargs)
                # say that we matched this input
                return True
        raise ValueError(f"Unrecognized input {_input} in state {frozen_state}")
