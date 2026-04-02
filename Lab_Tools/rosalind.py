"""
Rosalind Problem Collection Wrapper

A utility for fetching and working with problems from rosalind.info,
a platform for learning bioinformatics through problem solving.

Usage:
    from Lab_Tools.rosalind import Rosalind
    
    ros = Rosalind()
    
    # List available problems
    problems = ros.list_problems()
    
    # Get a specific problem
    problem = ros.get_problem("DNA")
    print(problem.title)
    print(problem.description)
    print(problem.sample_input)
    print(problem.sample_output)
    
    # Generate a notebook template
    ros.create_notebook("DNA", output_dir="./problems")
    
    # Test your solution
    def my_solution(data):
        # Your code here
        return result
    
    ros.test_solution("DNA", my_solution)
"""

import re
import json
import urllib.request
import urllib.error
from html.parser import HTMLParser
from dataclasses import dataclass, field
from typing import Optional, Callable, Any
from pathlib import Path


@dataclass
class RosalindProblem:
    """Represents a single Rosalind problem."""
    id: str
    title: str
    url: str
    description: str = ""
    problem_statement: str = ""
    sample_input: str = ""
    sample_output: str = ""
    solved_by: int = 0
    prerequisites: list = field(default_factory=list)
    glossary_terms: dict = field(default_factory=dict)
    
    def __str__(self) -> str:
        return f"{self.id}: {self.title}"
    
    def summary(self) -> str:
        """Return a formatted summary of the problem."""
        lines = [
            f"{'='*60}",
            f"{self.id}: {self.title}",
            f"{'='*60}",
            f"URL: {self.url}",
            f"Solved by: {self.solved_by:,} users",
            "",
        ]
        if self.prerequisites:
            lines.append(f"Prerequisites: {', '.join(self.prerequisites)}")
            lines.append("")
        if self.problem_statement:
            lines.append("PROBLEM:")
            lines.append(self.problem_statement)
            lines.append("")
        if self.sample_input:
            lines.append("SAMPLE INPUT:")
            lines.append(self.sample_input)
            lines.append("")
        if self.sample_output:
            lines.append("SAMPLE OUTPUT:")
            lines.append(self.sample_output)
        return "\n".join(lines)


class RosalindHTMLParser(HTMLParser):
    """Parse Rosalind problem pages to extract content."""
    
    def __init__(self):
        super().__init__()
        self.in_problem = False
        self.in_sample_dataset = False
        self.in_sample_output = False
        self.in_given_return = False
        self.in_blockquote = False
        self.in_h2 = False
        self.current_section = None
        self.content = {
            'title': '',
            'description': '',
            'problem': '',
            'sample_input': '',
            'sample_output': '',
            'solved_by': 0,
        }
        self.current_text = []
        
    def handle_starttag(self, tag, attrs):
        attrs_dict = dict(attrs)
        if tag == 'h2':
            self.in_h2 = True
        elif tag == 'blockquote':
            self.in_blockquote = True
        elif tag == 'code' or tag == 'pre':
            if self.in_sample_dataset:
                pass  # Will capture in handle_data
            elif self.in_sample_output:
                pass
                
    def handle_endtag(self, tag):
        if tag == 'h2':
            self.in_h2 = False
            text = ''.join(self.current_text).strip()
            self.current_text = []
            
            if 'Sample Dataset' in text:
                self.in_sample_dataset = True
                self.in_sample_output = False
                self.in_problem = False
            elif 'Sample Output' in text:
                self.in_sample_output = True
                self.in_sample_dataset = False
                self.in_problem = False
            elif 'Problem' in text:
                self.in_problem = True
                self.in_sample_dataset = False
                self.in_sample_output = False
            elif 'Given' in text or 'Return' in text:
                self.in_given_return = True
        elif tag == 'blockquote':
            self.in_blockquote = False
        elif tag == 'pre' or tag == 'code':
            text = ''.join(self.current_text).strip()
            if self.in_sample_dataset and text:
                self.content['sample_input'] = text
            elif self.in_sample_output and text:
                self.content['sample_output'] = text
            self.current_text = []
            
    def handle_data(self, data):
        if self.in_h2:
            self.current_text.append(data)
        elif self.in_sample_dataset or self.in_sample_output:
            self.current_text.append(data)


class Rosalind:
    """
    Interface for working with Rosalind bioinformatics problems.
    
    Examples:
        >>> ros = Rosalind()
        >>> problem = ros.get_problem("DNA")
        >>> print(problem.title)
        Counting DNA Nucleotides
        
        >>> # Test a solution
        >>> def count_nucleotides(dna):
        ...     return f"{dna.count('A')} {dna.count('C')} {dna.count('G')} {dna.count('T')}"
        >>> ros.test_solution("DNA", count_nucleotides)
        ✓ Solution passed for DNA!
    """
    
    BASE_URL = "https://rosalind.info"
    PROBLEMS_URL = f"{BASE_URL}/problems/list-view/"
    
    # Problem locations/tracks
    LOCATIONS = {
        'stronghold': '',  # Default bioinformatics problems
        'python-village': 'python-village',
        'armory': 'bioinformatics-armory', 
        'textbook': 'bioinformatics-textbook-track',
        'algorithmic': 'algorithmic-heights',
    }
    
    def __init__(self, cache_dir: Optional[Path] = None):
        """
        Initialize the Rosalind wrapper.
        
        Args:
            cache_dir: Optional directory to cache fetched problems.
        """
        self.cache_dir = Path(cache_dir) if cache_dir else None
        if self.cache_dir:
            self.cache_dir.mkdir(parents=True, exist_ok=True)
        self._problems_cache: dict[str, RosalindProblem] = {}
    
    def _fetch_url(self, url: str) -> str:
        """Fetch content from a URL."""
        headers = {
            'User-Agent': 'Mozilla/5.0 (compatible; RosalindWrapper/1.0; Python)'
        }
        req = urllib.request.Request(url, headers=headers)
        try:
            with urllib.request.urlopen(req, timeout=30) as response:
                return response.read().decode('utf-8')
        except urllib.error.HTTPError as e:
            raise ConnectionError(f"Failed to fetch {url}: HTTP {e.code}") from e
        except urllib.error.URLError as e:
            raise ConnectionError(f"Failed to fetch {url}: {e.reason}") from e
    
    def list_problems(self, location: str = 'stronghold') -> list[dict]:
        """
        List available problems from a specific track.
        
        Args:
            location: One of 'stronghold', 'python-village', 'armory', 
                     'textbook', or 'algorithmic'
        
        Returns:
            List of dicts with 'id', 'title', 'solved_by' keys.
        """
        import html
        
        loc = self.LOCATIONS.get(location, '')
        url = f"{self.PROBLEMS_URL}?location={loc}" if loc else self.PROBLEMS_URL
        
        content = self._fetch_url(url)
        problems = []
        
        # Parse problem list from HTML table
        # Pattern: href="/problems/ID/">Title</a>
        pattern = r'href="/problems/([a-z]+)/"[^>]*>([^<]+)</a>'
        matches = re.findall(pattern, content, re.IGNORECASE)
        
        # Find solved counts from "recent" links
        solved_pattern = r'href="/problems/([a-z]+)/recent/"[^>]*>(\d+)</a>'
        solved_matches = {pid.upper(): int(count) for pid, count in 
                        re.findall(solved_pattern, content, re.IGNORECASE)}
        
        seen = set()
        skip_titles = {'Topics', 'List View', 'Tree View'}
        
        for problem_id, title in matches:
            pid = problem_id.upper()
            # Clean up title: decode HTML entities, strip quotes and whitespace
            title = html.unescape(title).strip().strip('"')
            
            # Skip navigation links and duplicates
            if pid in seen or title in skip_titles or not title:
                continue
            if problem_id.lower() in ('topics', 'list', 'tree'):
                continue
                
            seen.add(pid)
            problems.append({
                'id': pid,
                'title': title,
                'solved_by': solved_matches.get(pid, 0),
                'url': f"{self.BASE_URL}/problems/{problem_id.lower()}/"
            })
        
        return problems
    
    def get_problem(self, problem_id: str) -> RosalindProblem:
        """
        Fetch detailed information about a specific problem.
        
        Args:
            problem_id: The problem ID (e.g., "DNA", "RNA", "GC")
        
        Returns:
            RosalindProblem with full details including sample I/O.
        """
        pid = problem_id.upper()
        
        # Check cache
        if pid in self._problems_cache:
            return self._problems_cache[pid]
        
        # Check file cache
        if self.cache_dir:
            cache_file = self.cache_dir / f"{pid}.json"
            if cache_file.exists():
                with open(cache_file) as f:
                    data = json.load(f)
                    problem = RosalindProblem(**data)
                    self._problems_cache[pid] = problem
                    return problem
        
        # Fetch from web
        url = f"{self.BASE_URL}/problems/{pid.lower()}/"
        content = self._fetch_url(url)
        
        problem = self._parse_problem_page(pid, url, content)
        self._problems_cache[pid] = problem
        
        # Save to file cache
        if self.cache_dir:
            cache_file = self.cache_dir / f"{pid}.json"
            with open(cache_file, 'w') as f:
                json.dump({
                    'id': problem.id,
                    'title': problem.title,
                    'url': problem.url,
                    'description': problem.description,
                    'problem_statement': problem.problem_statement,
                    'sample_input': problem.sample_input,
                    'sample_output': problem.sample_output,
                    'solved_by': problem.solved_by,
                    'prerequisites': problem.prerequisites,
                    'glossary_terms': problem.glossary_terms,
                }, f, indent=2)
        
        return problem
    
    def _parse_problem_page(self, problem_id: str, url: str, content: str) -> RosalindProblem:
        """Parse a problem page to extract all relevant information."""
        
        # Extract title from <h1> tag
        title_match = re.search(r'<h1>([^<]+)', content)
        title = title_match.group(1).strip() if title_match else problem_id
        
        # Extract solved count
        solved_match = re.search(r'solved by (\d+)', content)
        solved_by = int(solved_match.group(1)) if solved_match else 0
        
        # Extract description from blockquote
        desc_match = re.search(r'<blockquote>(.*?)</blockquote>', content, re.DOTALL)
        description = ""
        if desc_match:
            desc_html = desc_match.group(1)
            # Strip HTML tags but keep text
            desc_text = re.sub(r'<[^>]+>', ' ', desc_html)
            desc_text = re.sub(r'\s+', ' ', desc_text)
            description = desc_text.strip()[:1000]  # Limit length
        
        # Extract problem statement (between <h2 id="problem"> and <h2 id="sample-dataset">)
        problem_match = re.search(
            r'<h2 id="problem">.*?</h2>(.*?)(?=<h2 id="sample-dataset">)',
            content,
            re.DOTALL | re.IGNORECASE
        )
        problem_statement = ""
        if problem_match:
            ps_html = problem_match.group(1)
            # Clean up HTML
            ps_text = re.sub(r'<[^>]+>', ' ', ps_html)
            ps_text = re.sub(r'\s+', ' ', ps_text)
            problem_statement = ps_text.strip()
        
        # Extract sample input (in <div class="codehilite"><pre> after Sample Dataset)
        input_match = re.search(
            r'<h2 id="sample-dataset">.*?</h2>\s*<div class="codehilite"><pre>([^<]+)</pre>',
            content,
            re.DOTALL | re.IGNORECASE
        )
        sample_input = input_match.group(1).strip() if input_match else ""
        
        # Extract sample output (in <div class="codehilite"><pre> after Sample Output)
        output_match = re.search(
            r'<h2 id="sample-output">.*?</h2>\s*<div class="codehilite"><pre>([^<]+)</pre>',
            content,
            re.DOTALL | re.IGNORECASE
        )
        sample_output = output_match.group(1).strip() if output_match else ""
        
        # Extract glossary terms from tooltip titles
        glossary = {}
        glossary_matches = re.findall(
            r'href="/glossary/([^/]+)/"[^>]*title="[^:]*:\s*([^"]+)"',
            content
        )
        for slug, definition in glossary_matches:
            term = slug.replace('-', ' ').title()
            glossary[term] = definition.strip()
        
        # Extract prerequisites (not easily available on problem page)
        prereqs = []
        
        return RosalindProblem(
            id=problem_id,
            title=title,
            url=url,
            description=description,
            problem_statement=problem_statement,
            sample_input=sample_input,
            sample_output=sample_output,
            solved_by=solved_by,
            prerequisites=prereqs,
            glossary_terms=glossary,
        )
    
    def test_solution(
        self, 
        problem_id: str, 
        solution_func: Callable[[str], Any],
        verbose: bool = True
    ) -> bool:
        """
        Test a solution function against the sample input/output.
        
        Args:
            problem_id: The problem ID
            solution_func: Function that takes input string and returns result
            verbose: Whether to print results
        
        Returns:
            True if solution matches expected output, False otherwise.
        """
        problem = self.get_problem(problem_id)
        
        if not problem.sample_input or not problem.sample_output:
            if verbose:
                print(f"⚠ No sample data available for {problem_id}")
            return False
        
        try:
            result = solution_func(problem.sample_input)
            result_str = str(result).strip()
            expected = problem.sample_output.strip()
            
            # Normalize whitespace for comparison
            result_normalized = ' '.join(result_str.split())
            expected_normalized = ' '.join(expected.split())
            
            if result_normalized == expected_normalized:
                if verbose:
                    print(f"✓ Solution passed for {problem_id}!")
                return True
            else:
                if verbose:
                    print(f"✗ Solution failed for {problem_id}")
                    print(f"  Expected: {expected}")
                    print(f"  Got:      {result_str}")
                return False
                
        except Exception as e:
            if verbose:
                print(f"✗ Solution raised an error for {problem_id}: {e}")
            return False
    
    def create_notebook(
        self,
        problem_id: str,
        output_dir: str = ".",
        overwrite: bool = False
    ) -> Path:
        """
        Generate a Jupyter notebook template for a problem.
        
        Args:
            problem_id: The problem ID
            output_dir: Directory to save the notebook
            overwrite: Whether to overwrite existing notebooks
        
        Returns:
            Path to the created notebook.
        """
        problem = self.get_problem(problem_id)
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        
        notebook_path = output_path / f"rosalind_{problem_id.lower()}.ipynb"
        
        if notebook_path.exists() and not overwrite:
            raise FileExistsError(
                f"Notebook already exists: {notebook_path}. "
                "Use overwrite=True to replace."
            )
        
        # Build notebook structure
        cells = []
        
        # Title cell
        cells.append({
            "cell_type": "markdown",
            "metadata": {},
            "source": [
                f"# Rosalind: {problem.title}\n",
                "\n",
                f"**Problem ID:** {problem.id}  \n",
                f"**URL:** [{problem.url}]({problem.url})  \n",
                f"**Solved by:** {problem.solved_by:,} users\n",
            ]
        })
        
        # Problem description
        if problem.description:
            cells.append({
                "cell_type": "markdown", 
                "metadata": {},
                "source": [
                    "## Background\n",
                    "\n",
                    problem.description + "\n"
                ]
            })
        
        # Problem statement
        cells.append({
            "cell_type": "markdown",
            "metadata": {},
            "source": [
                "## Problem\n",
                "\n",
                problem.problem_statement + "\n"
            ]
        })
        
        # Sample data
        cells.append({
            "cell_type": "markdown",
            "metadata": {},
            "source": [
                "## Sample Data\n",
                "\n",
                "**Input:**\n",
                "```\n",
                problem.sample_input + "\n",
                "```\n",
                "\n",
                "**Output:**\n",
                "```\n",
                problem.sample_output + "\n",
                "```\n"
            ]
        })
        
        # Solution section
        cells.append({
            "cell_type": "markdown",
            "metadata": {},
            "source": ["## Solution\n"]
        })
        
        cells.append({
            "cell_type": "code",
            "metadata": {},
            "source": [
                f"def solve_{problem_id.lower()}(data: str) -> str:\n",
                '    """\n',
                f"    Solve Rosalind problem {problem_id}: {problem.title}\n",
                "    \n",
                "    Args:\n",
                "        data: Input string from Rosalind\n",
                "    \n",
                "    Returns:\n",
                "        Solution string to submit\n",
                '    """\n',
                "    # TODO: Implement your solution here\n",
                "    pass\n"
            ],
            "execution_count": None,
            "outputs": []
        })
        
        # Test with sample
        cells.append({
            "cell_type": "markdown",
            "metadata": {},
            "source": ["## Test with Sample Data\n"]
        })
        
        cells.append({
            "cell_type": "code",
            "metadata": {},
            "source": [
                "# Sample input\n",
                f'sample_input = """{problem.sample_input}"""\n',
                "\n",
                f'expected_output = """{problem.sample_output}"""\n',
                "\n",
                "# Test your solution\n",
                f"result = solve_{problem_id.lower()}(sample_input)\n",
                "print(f'Result: {result}')\n",
                "print(f'Expected: {expected_output}')\n",
                "print(f'Match: {str(result).strip() == expected_output.strip()}')\n"
            ],
            "execution_count": None,
            "outputs": []
        })
        
        # Run on actual data
        cells.append({
            "cell_type": "markdown",
            "metadata": {},
            "source": [
                "## Run on Rosalind Dataset\n",
                "\n",
                "Download your dataset from Rosalind and paste it below:\n"
            ]
        })
        
        cells.append({
            "cell_type": "code",
            "metadata": {},
            "source": [
                "# Paste your Rosalind dataset here\n",
                'rosalind_input = """PASTE_YOUR_DATA_HERE"""\n',
                "\n",
                f"answer = solve_{problem_id.lower()}(rosalind_input)\n",
                "print(answer)\n"
            ],
            "execution_count": None,
            "outputs": []
        })
        
        # Build notebook
        notebook = {
            "cells": cells,
            "metadata": {
                "kernelspec": {
                    "display_name": "Python 3",
                    "language": "python",
                    "name": "python3"
                },
                "language_info": {
                    "name": "python",
                    "version": "3.10.0"
                },
                "rosalind": {
                    "problem_id": problem.id,
                    "problem_title": problem.title,
                    "url": problem.url
                }
            },
            "nbformat": 4,
            "nbformat_minor": 4
        }
        
        with open(notebook_path, 'w') as f:
            json.dump(notebook, f, indent=2)
        
        return notebook_path
    
    def create_problem_set(
        self,
        problem_ids: list[str],
        output_dir: str = ".",
        readme: bool = True
    ) -> list[Path]:
        """
        Create notebooks for multiple problems.
        
        Args:
            problem_ids: List of problem IDs
            output_dir: Directory to save notebooks
            readme: Whether to create a README file
        
        Returns:
            List of paths to created notebooks.
        """
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        
        created = []
        problems_info = []
        
        for pid in problem_ids:
            try:
                path = self.create_notebook(pid, output_dir, overwrite=False)
                created.append(path)
                problem = self.get_problem(pid)
                problems_info.append(problem)
                print(f"✓ Created {path.name}")
            except FileExistsError:
                print(f"⊘ Skipped {pid} (already exists)")
            except Exception as e:
                print(f"✗ Failed {pid}: {e}")
        
        if readme and problems_info:
            readme_path = output_path / "README.md"
            with open(readme_path, 'w') as f:
                f.write("# Rosalind Problems\n\n")
                f.write("| # | ID | Problem | Solved By |\n")
                f.write("|---|----|---------|-----------|\n")
                for i, p in enumerate(problems_info, 1):
                    notebook_name = f"rosalind_{p.id.lower()}.ipynb"
                    f.write(f"| {i} | {p.id} | [{p.title}]({notebook_name}) | {p.solved_by:,} |\n")
                f.write(f"\n*Generated from [rosalind.info]({self.BASE_URL})*\n")
            created.append(readme_path)
            print(f"✓ Created README.md")
        
        return created


# Convenience functions for quick access
def fetch_problem(problem_id: str) -> RosalindProblem:
    """Quick fetch of a single problem."""
    return Rosalind().get_problem(problem_id)


def list_problems(location: str = 'stronghold') -> list[dict]:
    """Quick listing of problems."""
    return Rosalind().list_problems(location)


# Problem categories for course integration
COURSE_MAPPING = {
    # Tier 1 - Python Basics
    'python_basics': ['INI1', 'INI2', 'INI3', 'INI4', 'INI5', 'INI6'],
    
    # Tier 2 - Core Bioinformatics
    'dna_rna': ['DNA', 'RNA', 'REVC', 'GC', 'HAMM', 'SUBS'],
    'translation': ['PROT', 'MRNA', 'ORF', 'SPLC'],
    'sequence_analysis': ['CONS', 'LCSM', 'GRPH', 'LONG'],
    'motifs': ['MPRT', 'REVP', 'SSEQ', 'KMP'],
    
    # Tier 3 - Applied
    'phylogenetics': ['TREE', 'INOD', 'NWCK', 'CTBL'],
    'population_genetics': ['IPRB', 'IEV', 'LIA', 'AFRQ'],
    'assembly': ['LONG', 'CORR', 'DBRU'],
    
    # Tier 4 - Algorithms
    'dynamic_programming': ['EDIT', 'EDTA', 'LGIS', 'LCSQ'],
    'graph_algorithms': ['GRPH', 'LONG', 'DBRU', 'TREE'],
    'string_algorithms': ['KMP', 'TRIE', 'LREP', 'SUFF'],
}


# =============================================================================
# PROBLEM KNOWLEDGE BASE
# =============================================================================
# Detailed metadata for pedagogical guidance - concepts, hints, common mistakes

PROBLEM_KNOWLEDGE = {
    'DNA': {
        'concepts': ['string manipulation', 'counting', 'iteration'],
        'biology': ['nucleotides', 'DNA composition', 'base pairing'],
        'difficulty': 1,
        'hints': [
            "Think about what Python method counts occurrences of a character in a string.",
            "You need to count each of the four nucleotides: A, C, G, T.",
            "The output format matters - four integers separated by spaces.",
        ],
        'questions': [
            "What are the four nucleotides that make up DNA?",
            "How would you count how many times the letter 'A' appears in a string?",
            "What's the difference between a list and a string in Python?",
        ],
        'common_mistakes': [
            "Forgetting to handle newlines in input",
            "Wrong output format (commas instead of spaces)",
            "Case sensitivity issues",
        ],
        'next_problems': ['RNA', 'REVC'],
        'python_concepts': ['str.count()', 'string iteration', 'f-strings'],
    },
    'RNA': {
        'concepts': ['string manipulation', 'character replacement'],
        'biology': ['transcription', 'thymine vs uracil', 'DNA to RNA'],
        'difficulty': 1,
        'hints': [
            "Transcription replaces one nucleotide with another.",
            "Think about which base in DNA becomes which base in RNA.",
            "Python has a simple method for replacing characters in strings.",
        ],
        'questions': [
            "What is the biological process of creating RNA from DNA called?",
            "Which nucleotide in DNA is replaced during transcription, and with what?",
            "Does Python's replace() method modify the original string or return a new one?",
        ],
        'common_mistakes': [
            "Replacing the wrong nucleotide",
            "Forgetting that strings are immutable in Python",
        ],
        'next_problems': ['REVC', 'PROT'],
        'python_concepts': ['str.replace()', 'string immutability'],
    },
    'REVC': {
        'concepts': ['string reversal', 'complementation', 'mapping'],
        'biology': ['base pairing', 'complementary strands', 'antiparallel'],
        'difficulty': 1,
        'hints': [
            "This problem has two steps: complement AND reverse.",
            "A pairs with T, and G pairs with C.",
            "Think about the order - does it matter if you complement first or reverse first?",
        ],
        'questions': [
            "What does 'antiparallel' mean in the context of DNA strands?",
            "If you have 'ATGC', what is its complement?",
            "How do you reverse a string in Python?",
        ],
        'common_mistakes': [
            "Only complementing without reversing",
            "Only reversing without complementing",
            "Getting the complement pairs wrong",
        ],
        'next_problems': ['GC', 'HAMM'],
        'python_concepts': ['string slicing [::-1]', 'str.translate()', 'dictionary mapping'],
    },
    'GC': {
        'concepts': ['FASTA parsing', 'percentage calculation', 'max finding'],
        'biology': ['GC content', 'genome composition', 'FASTA format'],
        'difficulty': 2,
        'hints': [
            "First, you need to parse the FASTA format to separate sequences.",
            "GC content is (G + C) / total length * 100.",
            "You need to find the sequence with the HIGHEST GC content.",
        ],
        'questions': [
            "What does a FASTA header line start with?",
            "Why is GC content biologically significant?",
            "How would you split text by a delimiter in Python?",
        ],
        'common_mistakes': [
            "Not handling multi-line sequences in FASTA",
            "Integer division instead of float division",
            "Forgetting to multiply by 100 for percentage",
        ],
        'next_problems': ['HAMM', 'CONS'],
        'python_concepts': ['file parsing', 'str.split()', 'max() with key'],
    },
    'HAMM': {
        'concepts': ['string comparison', 'counting differences', 'zip'],
        'biology': ['point mutations', 'genetic distance', 'evolution'],
        'difficulty': 1,
        'hints': [
            "You need to compare two strings position by position.",
            "Count positions where the characters differ.",
            "Python's zip() function is useful for parallel iteration.",
        ],
        'questions': [
            "What is a point mutation in biology?",
            "How does Hamming distance relate to evolutionary distance?",
            "What does zip() do with two sequences of the same length?",
        ],
        'common_mistakes': [
            "Not handling strings of different lengths",
            "Off-by-one errors in indexing",
        ],
        'next_problems': ['SUBS', 'TRAN'],
        'python_concepts': ['zip()', 'list comprehension', 'sum()'],
    },
    'SUBS': {
        'concepts': ['substring search', 'multiple occurrences', '1-based indexing'],
        'biology': ['motifs', 'regulatory sequences', 'pattern recognition'],
        'difficulty': 2,
        'hints': [
            "You need to find ALL positions where the pattern occurs.",
            "Rosalind uses 1-based indexing, not 0-based.",
            "Overlapping occurrences count as separate matches.",
        ],
        'questions': [
            "What is a motif in bioinformatics?",
            "If 'ATA' appears at position 0 in 'ATATA', where else does it appear?",
            "What's the difference between 1-based and 0-based indexing?",
        ],
        'common_mistakes': [
            "Missing overlapping occurrences",
            "Forgetting to convert to 1-based indexing",
            "Using find() which only returns first occurrence",
        ],
        'next_problems': ['CONS', 'LCSM'],
        'python_concepts': ['str.find()', 'while loops', 'sliding window'],
    },
    'PROT': {
        'concepts': ['codon translation', 'genetic code', 'dictionary lookup'],
        'biology': ['translation', 'codons', 'amino acids', 'stop codons'],
        'difficulty': 2,
        'hints': [
            "RNA is read in groups of 3 nucleotides (codons).",
            "Each codon maps to an amino acid (or stop signal).",
            "You can use a dictionary to store the genetic code.",
        ],
        'questions': [
            "How many nucleotides make up a codon?",
            "What happens when a ribosome encounters a stop codon?",
            "How would you iterate through a string 3 characters at a time?",
        ],
        'common_mistakes': [
            "Including the stop codon in the output",
            "Using DNA codons instead of RNA codons",
            "Off-by-one errors in codon extraction",
        ],
        'next_problems': ['MRNA', 'ORF', 'SPLC'],
        'python_concepts': ['dictionary', 'range() with step', 'string slicing'],
    },
    'CONS': {
        'concepts': ['profile matrix', 'consensus sequence', '2D data'],
        'biology': ['multiple alignment', 'conservation', 'consensus'],
        'difficulty': 3,
        'hints': [
            "Build a count matrix: for each position, count each nucleotide.",
            "The consensus picks the most frequent nucleotide at each position.",
            "Make sure to parse FASTA format correctly first.",
        ],
        'questions': [
            "What does a consensus sequence represent biologically?",
            "How would you represent a matrix in Python?",
            "What if two nucleotides are tied for most frequent?",
        ],
        'common_mistakes': [
            "Incorrect matrix dimensions",
            "Output format not matching expected format",
            "Not handling ties consistently",
        ],
        'next_problems': ['PROB', 'LCSM'],
        'python_concepts': ['nested lists', 'collections.Counter', 'zip()'],
    },
    'IPRB': {
        'concepts': ['probability', 'Mendelian inheritance', 'combinatorics'],
        'biology': ['dominant/recessive', 'Punnett square', 'genotype'],
        'difficulty': 3,
        'hints': [
            "Calculate the probability of each possible mating.",
            "Use the law of total probability.",
            "Consider all possible pairs and their offspring probabilities.",
        ],
        'questions': [
            "What's the probability of a dominant phenotype from Aa × Aa?",
            "How do you calculate the probability of selecting two individuals?",
            "What does 'dominant' mean in Mendelian genetics?",
        ],
        'common_mistakes': [
            "Forgetting that selecting two individuals changes the pool",
            "Not considering all mating combinations",
            "Probability arithmetic errors",
        ],
        'next_problems': ['IEV', 'LIA'],
        'python_concepts': ['probability calculation', 'combinatorics'],
    },
    'GRPH': {
        'concepts': ['graph construction', 'suffix-prefix matching', 'output format'],
        'biology': ['overlap graphs', 'genome assembly', 'sequence similarity'],
        'difficulty': 3,
        'hints': [
            "An edge exists if the suffix of one equals the prefix of another.",
            "The overlap length k is given (usually 3).",
            "A sequence cannot have an edge to itself.",
        ],
        'questions': [
            "What is an overlap graph used for in bioinformatics?",
            "How do you get the last k characters of a string?",
            "Why can't a node connect to itself in this problem?",
        ],
        'common_mistakes': [
            "Including self-loops",
            "Wrong overlap length",
            "Output format issues",
        ],
        'next_problems': ['LONG', 'DBRU'],
        'python_concepts': ['string slicing', 'nested loops', 'set operations'],
    },
    'EDIT': {
        'concepts': ['dynamic programming', 'edit distance', '2D table'],
        'biology': ['sequence alignment', 'mutations', 'evolutionary distance'],
        'difficulty': 4,
        'hints': [
            "Build a 2D table where cell (i,j) is the edit distance of prefixes.",
            "Consider three operations: insert, delete, substitute.",
            "The recurrence relation is key - think about the last character.",
        ],
        'questions': [
            "What are the three edit operations?",
            "Why is dynamic programming more efficient than brute force here?",
            "What does cell (i,j) represent in the DP table?",
        ],
        'common_mistakes': [
            "Off-by-one errors in table indexing",
            "Incorrect base cases",
            "Wrong recurrence relation",
        ],
        'next_problems': ['EDTA', 'LCSQ'],
        'python_concepts': ['2D arrays', 'nested loops', 'min()'],
    },
}


# =============================================================================
# TUTOR CLASS - AI-ASSISTED LEARNING
# =============================================================================

class RosalindTutor:
    """
    AI-compatible tutor for Rosalind problems.
    
    Generates Socratic prompts that guide learning without giving answers.
    Designed to work with AI assistants like GitHub Copilot.
    
    Usage:
        tutor = RosalindTutor()
        
        # Get learning guidance (not answers!)
        guidance = tutor.get_guidance("DNA")
        print(guidance)
        
        # Get progressive hints
        hint = tutor.get_hint("DNA", level=1)
        
        # Generate AI prompt for assistance
        prompt = tutor.copilot_prompt("DNA", stuck_on="parsing input")
    """
    
    def __init__(self, rosalind: Optional[Rosalind] = None):
        self.rosalind = rosalind or Rosalind()
        self.knowledge = PROBLEM_KNOWLEDGE
    
    def get_knowledge(self, problem_id: str) -> dict:
        """Get pedagogical knowledge for a problem."""
        pid = problem_id.upper()
        return self.knowledge.get(pid, {})
    
    def get_guidance(self, problem_id: str) -> str:
        """
        Get learning guidance for a problem - concepts to understand,
        NOT the solution.
        """
        pid = problem_id.upper()
        problem = self.rosalind.get_problem(pid)
        knowledge = self.get_knowledge(pid)
        
        lines = [
            f"# Learning Guide: {problem.title} ({pid})",
            "",
            "## What You'll Learn",
            ""
        ]
        
        if knowledge.get('biology'):
            lines.append("**Biology concepts:**")
            for concept in knowledge['biology']:
                lines.append(f"  • {concept}")
            lines.append("")
        
        if knowledge.get('concepts'):
            lines.append("**Programming concepts:**")
            for concept in knowledge['concepts']:
                lines.append(f"  • {concept}")
            lines.append("")
        
        if knowledge.get('python_concepts'):
            lines.append("**Python tools you might use:**")
            for tool in knowledge['python_concepts']:
                lines.append(f"  • `{tool}`")
            lines.append("")
        
        if knowledge.get('questions'):
            lines.append("## Questions to Consider")
            lines.append("*Answer these before coding:*")
            lines.append("")
            for q in knowledge['questions']:
                lines.append(f"  1. {q}")
            lines.append("")
        
        if knowledge.get('common_mistakes'):
            lines.append("## Common Pitfalls")
            lines.append("*Watch out for:*")
            lines.append("")
            for mistake in knowledge['common_mistakes']:
                lines.append(f"  ⚠ {mistake}")
            lines.append("")
        
        if knowledge.get('next_problems'):
            lines.append("## After This Problem")
            lines.append(f"Try: {', '.join(knowledge['next_problems'])}")
        
        return '\n'.join(lines)
    
    def get_hint(self, problem_id: str, level: int = 1) -> str:
        """
        Get a progressive hint. Level 1 is vaguest, level 3 is most specific.
        Still doesn't give the answer!
        """
        pid = problem_id.upper()
        knowledge = self.get_knowledge(pid)
        hints = knowledge.get('hints', [])
        
        if not hints:
            return f"No hints available for {pid}. Try breaking down the problem step by step."
        
        # Clamp level to available hints
        level = max(1, min(level, len(hints)))
        
        return f"💡 Hint {level}/{len(hints)}: {hints[level - 1]}"
    
    def copilot_prompt(
        self, 
        problem_id: str, 
        stuck_on: Optional[str] = None,
        code_so_far: Optional[str] = None,
        error_message: Optional[str] = None
    ) -> str:
        """
        Generate a prompt for AI assistants that encourages learning
        without giving away the answer.
        
        Args:
            problem_id: The Rosalind problem ID
            stuck_on: What the student is stuck on
            code_so_far: Student's current code attempt
            error_message: Any error they're encountering
        
        Returns:
            A prompt string to paste to an AI assistant
        """
        pid = problem_id.upper()
        problem = self.rosalind.get_problem(pid)
        knowledge = self.get_knowledge(pid)
        
        prompt_parts = [
            "=" * 60,
            "ROSALIND LEARNING ASSISTANT - SOCRATIC MODE",
            "=" * 60,
            "",
            "I am working on a Rosalind bioinformatics problem and need",
            "GUIDANCE, not a complete solution. Please help me LEARN by:",
            "",
            "  ✓ Asking me clarifying questions",
            "  ✓ Explaining concepts I'm missing",
            "  ✓ Pointing out what to think about",
            "  ✓ Suggesting what to Google/research",
            "  ✓ Reviewing my approach (not writing code for me)",
            "",
            "  ✗ Do NOT write the complete solution",
            "  ✗ Do NOT give me code to copy-paste", 
            "  ✗ Do NOT solve the problem for me",
            "",
            "-" * 60,
            f"PROBLEM: {problem.title} ({pid})",
            "-" * 60,
            "",
            problem.problem_statement[:500] + "..." if len(problem.problem_statement) > 500 else problem.problem_statement,
            "",
        ]
        
        if knowledge.get('concepts'):
            prompt_parts.extend([
                "Key concepts involved:",
                f"  {', '.join(knowledge['concepts'])}",
                ""
            ])
        
        if stuck_on:
            prompt_parts.extend([
                f"WHAT I'M STUCK ON: {stuck_on}",
                ""
            ])
        
        if code_so_far:
            prompt_parts.extend([
                "MY CURRENT ATTEMPT:",
                "```python",
                code_so_far,
                "```",
                ""
            ])
        
        if error_message:
            prompt_parts.extend([
                f"ERROR I'M GETTING: {error_message}",
                ""
            ])
        
        prompt_parts.extend([
            "-" * 60,
            "Please help me understand what I'm missing, but let me",
            "write the solution myself. Ask me questions to guide my thinking.",
            "=" * 60,
        ])
        
        return '\n'.join(prompt_parts)
    
    def review_prompt(self, problem_id: str, code: str) -> str:
        """
        Generate a prompt for AI to review student code without
        giving away improvements unless asked.
        """
        pid = problem_id.upper()
        problem = self.rosalind.get_problem(pid)
        knowledge = self.get_knowledge(pid)
        
        prompt = f"""
{'=' * 60}
ROSALIND CODE REVIEW - LEARNING MODE
{'=' * 60}

I've written a solution for Rosalind problem {pid}: {problem.title}

Please review my code using the SOCRATIC METHOD:

  ✓ Ask questions about my design choices
  ✓ Point out potential issues as QUESTIONS, not fixes
  ✓ Help me discover improvements myself
  ✓ Check if I've considered edge cases (ask, don't tell)
  
  ✗ Do NOT rewrite my code
  ✗ Do NOT give me the "correct" solution
  ✗ Do NOT fix bugs directly - help me find them

{'-' * 60}
MY CODE:
{'-' * 60}

```python
{code}
```

{'-' * 60}
REVIEW QUESTIONS TO CONSIDER:
{'-' * 60}

1. Does my code handle all edge cases?
2. Is my approach efficient enough for large inputs?
3. Did I understand the problem correctly?
4. Are there any bugs I should look for?

Common mistakes for this problem:
{chr(10).join('  • ' + m for m in knowledge.get('common_mistakes', ['(no specific mistakes documented)']))}

{'=' * 60}
"""
        return prompt
    
    def explain_concept(self, concept: str) -> str:
        """
        Generate a prompt asking AI to explain a concept
        in the context of bioinformatics.
        """
        return f"""
{'=' * 60}
CONCEPT EXPLANATION REQUEST
{'=' * 60}

Please explain the concept of "{concept}" for a bioinformatics student.

Include:
1. What it is (simple definition)
2. Why it matters in biology/bioinformatics  
3. A simple example
4. How it's typically implemented in Python

Keep the explanation educational - I want to understand it,
not just use it mechanically.

{'=' * 60}
"""
    
    def debug_prompt(self, problem_id: str, code: str, error: str) -> str:
        """
        Generate a prompt for AI-assisted debugging that
        teaches debugging skills rather than fixing the bug.
        """
        pid = problem_id.upper()
        
        return f"""
{'=' * 60}
ROSALIND DEBUGGING ASSISTANT - TEACHING MODE
{'=' * 60}

I'm debugging my solution for Rosalind problem {pid}.

IMPORTANT: Help me LEARN to debug, don't just fix it!

  ✓ Ask what I've tried so far
  ✓ Suggest debugging strategies (print statements, etc.)
  ✓ Help me narrow down where the bug might be
  ✓ Explain how to read the error message
  ✓ Guide me to find the bug myself

  ✗ Do NOT tell me exactly what's wrong
  ✗ Do NOT write the fix for me
  ✗ Do NOT give away the answer

{'-' * 60}
MY CODE:
{'-' * 60}

```python
{code}
```

{'-' * 60}
THE ERROR:
{'-' * 60}

{error}

{'-' * 60}
Help me develop my debugging skills by guiding my investigation!
{'=' * 60}
"""
    
    def learning_path(self, start_problem: str = 'DNA') -> list[dict]:
        """
        Generate a recommended learning path starting from a problem.
        """
        path = []
        visited = set()
        queue = [start_problem.upper()]
        
        while queue and len(path) < 20:
            pid = queue.pop(0)
            if pid in visited:
                continue
            visited.add(pid)
            
            knowledge = self.get_knowledge(pid)
            if not knowledge:
                continue
            
            try:
                problem = self.rosalind.get_problem(pid)
                path.append({
                    'id': pid,
                    'title': problem.title,
                    'difficulty': knowledge.get('difficulty', 1),
                    'concepts': knowledge.get('concepts', []),
                    'biology': knowledge.get('biology', []),
                })
                
                # Add next problems to queue
                for next_pid in knowledge.get('next_problems', []):
                    if next_pid not in visited:
                        queue.append(next_pid)
            except Exception:
                continue
        
        # Sort by difficulty
        path.sort(key=lambda x: x['difficulty'])
        return path
    
    def create_tutor_notebook(
        self,
        problem_id: str,
        output_dir: str = ".",
        overwrite: bool = False
    ) -> Path:
        """
        Create a Jupyter notebook with AI-assisted learning scaffolds.
        Includes Copilot-compatible prompts and learning guidance.
        """
        pid = problem_id.upper()
        problem = self.rosalind.get_problem(pid)
        knowledge = self.get_knowledge(pid)
        
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        notebook_path = output_path / f"learn_{pid.lower()}.ipynb"
        
        if notebook_path.exists() and not overwrite:
            raise FileExistsError(f"Notebook exists: {notebook_path}")
        
        cells = []
        
        # Title with learning mode indicator
        cells.append({
            "cell_type": "markdown",
            "metadata": {},
            "source": [
                f"# 🎓 Rosalind Learning Mode: {problem.title}\n",
                "\n",
                f"**Problem ID:** {pid}  \n",
                f"**Difficulty:** {'⭐' * knowledge.get('difficulty', 1)}  \n",
                f"**URL:** [{problem.url}]({problem.url})\n",
                "\n",
                "---\n",
                "\n",
                "> 💡 **Learning Mode Active**  \n",
                "> This notebook is designed to help you LEARN, not just solve.  \n",
                "> Use the AI prompts to get guidance without spoilers!\n",
            ]
        })
        
        # Learning objectives
        cells.append({
            "cell_type": "markdown",
            "metadata": {},
            "source": [
                "## 📚 What You'll Learn\n",
                "\n",
                "**Biology:**\n",
            ] + [f"- {c}\n" for c in knowledge.get('biology', ['(explore the problem)'])] +
            [
                "\n**Programming:**\n",
            ] + [f"- {c}\n" for c in knowledge.get('concepts', ['(discover through solving)'])]
        })
        
        # Pre-flight questions
        if knowledge.get('questions'):
            cells.append({
                "cell_type": "markdown",
                "metadata": {},
                "source": [
                    "## ✋ Before You Code\n",
                    "\n",
                    "*Answer these questions first (write your answers below):*\n",
                    "\n",
                ] + [f"{i}. {q}\n" for i, q in enumerate(knowledge['questions'], 1)]
            })
            
            cells.append({
                "cell_type": "markdown",
                "metadata": {},
                "source": [
                    "**Your answers:**\n",
                    "\n",
                    "1. \n",
                    "2. \n",
                    "3. \n",
                ]
            })
        
        # Problem statement
        cells.append({
            "cell_type": "markdown",
            "metadata": {},
            "source": [
                "## 📋 Problem\n",
                "\n",
                problem.problem_statement + "\n",
            ]
        })
        
        # Sample data
        cells.append({
            "cell_type": "markdown",
            "metadata": {},
            "source": [
                "## 📊 Sample Data\n",
                "\n",
                "**Input:**\n",
                "```\n",
                problem.sample_input + "\n",
                "```\n",
                "\n",
                "**Expected Output:**\n",
                "```\n",
                problem.sample_output + "\n",
                "```\n",
            ]
        })
        
        # Hints section (collapsed)
        cells.append({
            "cell_type": "markdown",
            "metadata": {},
            "source": [
                "## 💡 Hints (Reveal Gradually)\n",
                "\n",
                "*Try to solve without hints first! Click to reveal only if stuck.*\n",
                "\n",
                "<details>\n",
                "<summary>Hint 1 (Gentlest)</summary>\n",
                "\n",
                knowledge.get('hints', ['Think about breaking the problem into steps.'])[0] + "\n",
                "\n",
                "</details>\n",
            ] + (
                [
                    "\n<details>\n",
                    "<summary>Hint 2 (More Specific)</summary>\n",
                    "\n",
                    knowledge.get('hints', ['', 'Consider the data structures you need.'])[1] + "\n" if len(knowledge.get('hints', [])) > 1 else "",
                    "\n</details>\n",
                ] if len(knowledge.get('hints', [])) > 1 else []
            ) + (
                [
                    "\n<details>\n",
                    "<summary>Hint 3 (Most Direct)</summary>\n",
                    "\n",
                    knowledge.get('hints', ['', '', 'Almost there...'])[2] + "\n" if len(knowledge.get('hints', [])) > 2 else "",
                    "\n</details>\n",
                ] if len(knowledge.get('hints', [])) > 2 else []
            )
        })
        
        # AI assistance prompt
        cells.append({
            "cell_type": "markdown",
            "metadata": {},
            "source": [
                "## 🤖 AI Learning Assistant\n",
                "\n",
                "If you're stuck, copy the prompt below to get GUIDANCE (not answers):\n",
                "\n",
                "```\n",
                self.copilot_prompt(pid, stuck_on="[describe what you're stuck on]").replace('```', '~~~') + "\n",
                "```\n",
            ]
        })
        
        # Solution workspace
        cells.append({
            "cell_type": "markdown",
            "metadata": {},
            "source": ["## ✏️ Your Solution\n"]
        })
        
        cells.append({
            "cell_type": "code",
            "metadata": {},
            "source": [
                "# Step 1: Understand the input format\n",
                "# What does the input look like? What do you need to extract?\n",
                "\n",
                f'sample_input = """{problem.sample_input}"""\n',
                "\n",
                "# Parse the input here:\n",
                "\n",
            ],
            "execution_count": None,
            "outputs": []
        })
        
        cells.append({
            "cell_type": "code",
            "metadata": {},
            "source": [
                "# Step 2: Implement your solution logic\n",
                "# What's your algorithm? Write it step by step.\n",
                "\n",
                f"def solve_{pid.lower()}(data: str) -> str:\n",
                '    """Your solution here."""\n',
                "    # TODO: Implement\n",
                "    pass\n",
            ],
            "execution_count": None,
            "outputs": []
        })
        
        cells.append({
            "cell_type": "code",
            "metadata": {},
            "source": [
                "# Step 3: Test with sample data\n",
                "\n",
                f"result = solve_{pid.lower()}(sample_input)\n",
                f'expected = """{problem.sample_output}"""\n',
                "\n",
                "print(f'Your result: {result}')\n",
                "print(f'Expected:    {expected}')\n",
                "print(f'Match: {str(result).strip() == expected.strip()}')\n",
            ],
            "execution_count": None,
            "outputs": []
        })
        
        # Reflection section
        cells.append({
            "cell_type": "markdown",
            "metadata": {},
            "source": [
                "## 🔍 Reflection\n",
                "\n",
                "*After solving, answer these:*\n",
                "\n",
                "1. What was the trickiest part?\n",
                "2. What would you do differently next time?\n",
                "3. What concepts do you need to review?\n",
                "\n",
                "**Your reflection:**\n",
                "\n",
            ]
        })
        
        # Next steps
        if knowledge.get('next_problems'):
            cells.append({
                "cell_type": "markdown",
                "metadata": {},
                "source": [
                    "## ➡️ Next Problems\n",
                    "\n",
                    "Ready for more? Try these:\n",
                ] + [f"- **{p}** - builds on what you learned here\n" 
                     for p in knowledge['next_problems']]
            })
        
        # Build notebook
        notebook = {
            "cells": cells,
            "metadata": {
                "kernelspec": {
                    "display_name": "Python 3",
                    "language": "python",
                    "name": "python3"
                },
                "rosalind_tutor": {
                    "problem_id": pid,
                    "mode": "learning",
                    "difficulty": knowledge.get('difficulty', 1),
                }
            },
            "nbformat": 4,
            "nbformat_minor": 4
        }
        
        with open(notebook_path, 'w') as f:
            json.dump(notebook, f, indent=2)
        
        return notebook_path


# =============================================================================
# COPILOT INTEGRATION HELPERS  
# =============================================================================

def get_tutor() -> RosalindTutor:
    """Get a tutor instance for quick access."""
    return RosalindTutor()


def copilot_help(problem_id: str, stuck_on: str = None) -> str:
    """
    Quick helper to generate a Copilot-compatible learning prompt.
    
    Usage:
        print(copilot_help("DNA", stuck_on="parsing the input"))
    """
    tutor = RosalindTutor()
    return tutor.copilot_prompt(problem_id, stuck_on=stuck_on)


def learning_mode(problem_id: str) -> str:
    """
    Get comprehensive learning guidance for a problem.
    
    Usage:
        print(learning_mode("GC"))
    """
    tutor = RosalindTutor()
    return tutor.get_guidance(problem_id)


if __name__ == "__main__":
    # Demo usage
    print("Rosalind Problem Wrapper Demo\n")
    
    ros = Rosalind()
    
    # List first 10 problems
    print("First 10 Bioinformatics Stronghold problems:")
    print("-" * 50)
    problems = ros.list_problems()[:10]
    for p in problems:
        print(f"  {p['id']:6s} - {p['title']} ({p['solved_by']:,} solved)")
    
    print("\n" + "=" * 50)
    
    # Fetch DNA problem details
    print("\nFetching 'DNA' problem details...")
    dna = ros.get_problem("DNA")
    print(dna.summary())
