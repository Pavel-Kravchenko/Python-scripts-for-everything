"""
Deterministic Finite Automaton (DFA) pattern matching.

Two implementations are provided:

1. Functional approach (avtomat_dictionary_maker / avtomat_search):
   Builds a transition table as a dict-of-dicts keyed by state index,
   then scans the text through the automaton.

2. Class-based approach (DFAMatcher):
   Each state is a node object with a transition dict; the automaton is
   constructed by the same prefix-function trick and searching is done
   via an instance method.

Both rely on the same helper: prefix_length, which computes how many
characters of the pattern are matched by the suffix of a given probe string.
This is equivalent to running the prefix (KMP) function.
"""
from __future__ import annotations


def prefix_length(pattern: str, probe: str) -> int:
    """Return the length of the longest prefix of pattern that is also a
    suffix of probe, using the KMP prefix function approach.

    This is used to determine automaton transition targets.
    """
    combined = pattern + "#" + probe
    j = 0
    i = 1
    sp = [0] * len(combined)
    combined = combined + "$"
    while combined[i] != "$":
        while j > 0 and combined[i] != combined[j]:
            j = sp[j - 1]
        if combined[i] == combined[j]:
            j += 1
        sp[i] = j
        i += 1
    return sp[-1]


# ---------------------------------------------------------------------------
# Functional approach
# ---------------------------------------------------------------------------

def build_automaton(pattern: str, alphabet: str) -> dict[int, dict[str, int]]:
    """Build a DFA transition table for the given pattern over alphabet.

    Returns a dict mapping state index -> {character -> next_state}.
    State 0 is the start state; state len(pattern) is the accept state.
    """
    automaton: dict[int, dict[str, int]] = {}
    for i in range(len(pattern) + 1):
        prefix = pattern[:i]
        transitions: dict[str, int] = {}
        for char in alphabet:
            probe = prefix + char
            transitions[char] = prefix_length(pattern, probe)
        automaton[i] = transitions
    return automaton


def automaton_search(text: str, automaton: dict[int, dict[str, int]],
                     pattern_len: int) -> list[int]:
    """Scan text through the prebuilt automaton and return all match positions.

    Args:
        text: The text to search in.
        automaton: Transition table from build_automaton.
        pattern_len: Length of the pattern (accept state index).

    Returns:
        List of 0-based start positions where pattern was found.
    """
    state = 0
    matches = []
    states = list(automaton.keys())
    for i, char in enumerate(text):
        state = automaton[state][char]
        if state == pattern_len:
            matches.append(i - pattern_len + 1)
    return matches


# ---------------------------------------------------------------------------
# Class-based approach
# ---------------------------------------------------------------------------

class _DFANode:
    """Single state in the DFA; stores transitions as a dict."""

    def __init__(self, transitions: dict[str, int]):
        self.dest = transitions


class DFAMatcher:
    """DFA-based exact pattern matcher for a fixed DNA alphabet (A, T, G, C).

    Build once per pattern, then call search() on any number of texts.
    """

    ALPHABET = ("A", "T", "G", "C")

    def __init__(self, pattern: str):
        self.pattern = pattern
        self.nodes: list[_DFANode] = []
        for i in range(len(pattern) + 1):
            prefix = pattern[:i]
            transitions = {
                char: prefix_length(pattern, prefix + char)
                for char in self.ALPHABET
            }
            self.nodes.append(_DFANode(transitions))

    def search(self, text: str) -> list[int]:
        """Return all 0-based start positions of pattern in text."""
        state = 0
        accept = len(self.nodes) - 1
        matches = []
        for i, char in enumerate(text):
            state = self.nodes[state].dest[char]
            if state == accept:
                matches.append(i - len(self.pattern) + 1)
        return matches


if __name__ == "__main__":
    alphabet = "ATGC"
    pattern = "ATTCTGATTT"
    text = "AATGCCGTATTCTATTCTGATTTCTGAATTCTGATTTTTAGT"

    # Functional approach
    automaton = build_automaton(pattern, alphabet)
    matches = automaton_search(text, automaton, len(pattern))
    print(f"Functional DFA matches of '{pattern}': {matches}")

    # Class-based approach
    dfa = DFAMatcher(pattern)
    print(f"Class DFA matches of '{pattern}': {dfa.search(text)}")

    # Another example
    dfa2 = DFAMatcher("ATCCG")
    print(f"Class DFA matches of 'ATCCG' in 'ATTGATCCGACGATCCTG': {dfa2.search('ATTGATCCGACGATCCTG')}")
