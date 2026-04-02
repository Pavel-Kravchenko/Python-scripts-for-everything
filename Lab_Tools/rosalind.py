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
