"""
Naive string search and naive prefix function.

Includes:
- naive_search: brute-force substring check at every position
- naive_prefix_function: O(n^3) prefix array built by comparing all
  prefix-suffix pairs explicitly (educational, not for production)
"""
from __future__ import annotations


def naive_search(text: str, pattern: str) -> str:
    """Return 'Yes' if pattern appears in text, 'No' otherwise.

    Checks every possible starting position by slicing.
    """
    for i in range(len(text) - len(pattern) + 1):
        if text[i:i + len(pattern)] == pattern:
            return "Yes"
    return "No"


def naive_prefix_function(s: str) -> list[int]:
    """Build the prefix (failure) array for s using a naive O(n^3) approach.

    For each position i, scan all proper prefixes of s[0:i] and find the
    longest one that is also a suffix.

    sp[i] is the length of that longest proper prefix-suffix.
    """
    sp = [0] * len(s)
    for i in range(2, len(s)):
        sub = s[0:i]
        for j in range(len(sub)):
            suffix = sub[-j:] if j > 0 else ""
            prefix = sub[:j]
            if suffix == prefix and j > 0:
                sp[i] = j
    return sp


if __name__ == "__main__":
    # Naive search on DNA sequences
    text = "AATGCATGCTAGCATGC"
    patterns = ["ATGC", "TTTT", "GCTA", "CATGC"]
    for p in patterns:
        print(f"naive_search('{p}'): {naive_search(text, p)}")

    # Naive prefix function
    print(f"\nnaive_prefix_function('ATGCATGC'): {naive_prefix_function('ATGCATGC')}")
    print(f"naive_prefix_function('AAAA'):     {naive_prefix_function('AAAA')}")
