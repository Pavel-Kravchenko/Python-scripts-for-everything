"""
Knuth-Morris-Pratt (KMP) pattern matching algorithm.

Includes:
- prefix_function: builds the failure/prefix array for a pattern
- kmp_search: finds all (overlapping) occurrences of a pattern in a text
- prefix_search: variant that concatenates pattern + separator + text and
  uses the prefix array directly to locate matches
"""
from __future__ import annotations


def prefix_function(p: str) -> list[int]:
    """Build the KMP prefix (failure) array for pattern p.

    sp[i] is the length of the longest proper prefix of p[0:i+1]
    that is also a suffix.
    """
    sp = [0] * len(p)
    j = 0
    for i in range(1, len(p)):
        while j >= 0 and p[j] != p[i]:
            if j - 1 >= 0:
                j = sp[j - 1]
            else:
                j = -1
        j += 1
        sp[i] = j
    return sp


def kmp_search(text: str, pattern: str) -> list[int]:
    """Return start positions of all (possibly overlapping) occurrences of
    pattern in text using the KMP algorithm.
    """
    matches = []
    f = prefix_function(pattern)
    n, m = len(text), len(pattern)
    j = 0
    for i in range(n):
        while j >= 0 and text[i] != pattern[j]:
            if j - 1 >= 0:
                j = f[j - 1]
            else:
                j = -1
        j += 1
        if j == m:
            # Full match: record start position
            matches.append(i - m + 1)
            # Use failure function to allow overlapping matches
            j = f[m - 1]
    return matches


def prefix_search(text: str, pattern: str) -> list[int]:
    """Find all occurrences of pattern in text by running the prefix function
    on the concatenated string  pattern + '@' + text + '$'.

    The sentinel characters '@' and '$' must not appear in the pattern or text.
    Returns a list of 0-based start positions in text.
    """
    combined = pattern + "@" + text + "$"
    sp = [0] * len(combined)
    j = 0
    i = 1
    # Stop at the sentinel '$' so we never read past it
    while combined[i] != "$":
        while j > 0 and combined[i] != combined[j]:
            j = sp[j - 1]
        if combined[i] == combined[j]:
            j += 1
        sp[i] = j
        i += 1

    m = len(pattern)
    out = []
    for idx in range(len(sp)):
        if sp[idx] == m:
            # Subtract pattern length twice: once for the pattern prefix,
            # once for the separator character
            out.append(idx - m * 2)
    return out


if __name__ == "__main__":
    # Example with DNA sequences
    text = "ATGCATGCATGC"
    pattern = "ATGC"

    positions = kmp_search(text, pattern)
    print(f"kmp_search positions of '{pattern}' in '{text}': {positions}")

    positions2 = prefix_search(text, pattern)
    print(f"prefix_search positions of '{pattern}' in '{text}': {positions2}")

    # Overlapping example
    text2 = "AAAA"
    pattern2 = "AA"
    print(f"kmp_search (overlapping) '{pattern2}' in '{text2}': {kmp_search(text2, pattern2)}")
