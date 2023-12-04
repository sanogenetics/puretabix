"""
Finite State Machine
--------------------

Simple implementation of a finite state machine, used for parsing in a
character-by-character manner


based on https://gist.github.com/brianray/8d3e697dbbf150f725291d74ac0cee8b
"""

import logging
import re
from typing import (
    Any,
    Callable,
    Dict,
    FrozenSet,
    Iterable,
    Mapping,
    Optional,
    Pattern,
    Type,
    Union,
)

logger = logging.getLogger(__name__)

# note: Optional[Callable[[Any, str], None]] isn't the optimal type hint for the callbacks
# In theory, they can take an arbitrary number of parameters as long as the last position is a str
# and the others are matched via callback_args/callback_kwargs
# This might be able to be captured in newer Python versions with Protocol and/or ParamSpec


class RegexTransition:
    dst: Any
    condition: Pattern[str]
    callback: Optional[Callable[[Any, str], None]]
    __slots__ = ["dst", "condition", "callback"]

    def __init__(
        self,
        destination_state: Any,
        condition: str,
        callback: Callable[[Any, str], None],
    ):
        self.dst = destination_state
        self.condition = re.compile(condition)
        self.callback = callback

    def match(self, _input: str) -> bool:
        return bool(self.condition.match(_input))


class SetInTransition:
    dst: Any
    condition: FrozenSet[str]
    callback: Optional[Callable[[Any, str], None]]
    __slots__ = ["dst", "condition", "callback"]

    def __init__(
        self,
        destination_state: Any,
        condition: Iterable[str],
        callback: Optional[Callable[[Any, str], None]],
    ):
        self.dst = destination_state
        self.condition = frozenset(condition)
        self.callback = callback

    def match(self, _input: Union[str, bytes]) -> bool:
        return bool(_input in self.condition)


class SetNotInTransition:
    dst: Any
    condition: FrozenSet[str]
    callback: Optional[Callable[[Any, str], None]]

    __slots__ = ["dst", "condition", "callback"]

    def __init__(
        self,
        destination_state: Any,
        condition: Iterable[str],
        callback: Optional[Callable[[Any, str], None]],
    ):
        self.dst = destination_state
        self.condition = frozenset(condition)
        self.callback = callback

    def match(self, _input: Union[str, bytes]) -> bool:
        return bool(_input not in self.condition)


class FSMachine:
    transitions: Dict[Any, Any]

    def __init__(self) -> None:
        self.transitions = {}

    def add_transition(
        self,
        start_state: Any,
        end_state: Any,
        transition_class: Type[Any],
        condition: Any,
        callback: Optional[Callable[[Any, str], None]],
    ) -> None:
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

    def run(
        self,
        inputs: Iterable[str],
        initial_state: Any,
        *args: Any,
        **kwargs: Dict[Any, Any],
    ) -> None:
        self.current_state = initial_state
        for c in inputs:
            self.process_next(c, args, kwargs)
            # if state is None, early exit
            if not self.current_state:
                break

        # process that we reached the end of the input
        if self.current_state:
            self.process_next(None, args, kwargs)

        # check at a valid end state
        assert not self.current_state, f"Unexpected ending at {self.current_state}"

    def process_next(
        self,
        _input: Optional[str],
        callback_args: Any,
        callback_kwargs: Mapping[Any, Any],
    ) -> bool:
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
