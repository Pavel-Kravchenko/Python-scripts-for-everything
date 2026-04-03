"""
Rabin-Karp pattern matching algorithm with a custom rolling hash.

The hash function maps each character to its digit by subtracting the ASCII
offset for lowercase letters (ord(c) - 96), then concatenates those single
digits into an integer.  Rolling is done by dropping the leading digit via
modulo and appending the new character's digit.

This is a simplified, educational implementation — not a production-grade
polynomial rolling hash.
"""
from __future__ import annotations


def make_hash(string: str, rehash: list | None = None) -> int:
    """Compute an integer hash for a string, or roll the hash one step.

    Args:
        string: When rehash is None, the full substring to hash.
                When rehash is provided, the single new character to append.
        rehash: [old_hash_value, pattern_length] used for rolling.

    Returns:
        Integer hash of the (new) window.
    """
    if rehash is None:
        # Build hash by concatenating each character's digit
        out_str = ""
        for ch in string:
            out_str += str(ord(ch) - 96)[0]
        return int(out_str)
    else:
        old_hash, pattern_len = rehash[0], rehash[1]
        # Drop the leading digit: keep only the last (pattern_len - 1) digits
        if pattern_len > 1:
            modulus = int("1" + "0" * (pattern_len - 1))
        else:
            modulus = 1
        out_str = str(old_hash % modulus) + str(ord(string) - 96)[0]
        return int(out_str)


def rabin_karp(pattern: str, text: str) -> int:
    """Return the start index of the first occurrence of pattern in text,
    or -1 if not found.

    Uses a rolling hash to avoid recomputing the hash from scratch at each
    position.
    """
    if len(pattern) > len(text) or len(pattern) == 0 or len(text) == 0:
        return -1

    pattern_hash = make_hash(pattern)
    text_hash = make_hash(text[0:len(pattern)])

    for i in range(0, len(text) - len(pattern) + 1):
        if text_hash == pattern_hash:
            return i
        # Roll the hash: drop the leftmost character, append the next one
        if i + len(pattern) < len(text):
            text_hash = make_hash(text[i + len(pattern)], rehash=[text_hash, len(pattern)])

    return -1


if __name__ == "__main__":
    # Note: the hash function uses ord(c) - 96, which maps lowercase ASCII
    # letters a-z to digits 1-26.  Uppercase letters and digits produce
    # negative or out-of-range offsets, so inputs should use lowercase.
    print(rabin_karp("gc", "aaaaatgcatgc"))           # expected: 5
    print(rabin_karp("atgc", "ttttatgcgggg"))         # expected: 4
    print(rabin_karp("aaaa", "tttgatgcgggg"))         # expected: -1
    # Original examples from source
    print(rabin_karp("gc", "aaaaatgcatgc"))
    print(rabin_karp("ccccccggagagagcgc",
                     "aaagcgcgcgcgccccccggagagagcgcagcagcaatgcatgc"))
