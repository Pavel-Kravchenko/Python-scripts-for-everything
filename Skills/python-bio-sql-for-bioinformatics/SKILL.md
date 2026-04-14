---
name: python-bio-sql-for-bioinformatics
description: Ensembl, UCSC Genome Browser, NCBI, and dbSNP are all backed by relational databases. Even locally, SQLite is a practical way to store and query gene annotations, variant tables, and expression result
tool_type: python
primary_tool: Pandas
---

# SQL for Bioinformatics

## Complicated moments explained

**1. JOIN type determines which rows survive**
- `INNER JOIN`: only rows with matches in both tables.
- `LEFT JOIN`: all rows from the left table; NULL for unmatched right-side columns.

**2. `HAVING` is not the same as `WHERE`**
`WHERE` filters rows *before* grouping. `HAVING` filters *after* grouping. To filter groups by aggregate values (e.g., "genes with more than 3 variants"), use `HAVING COUNT(*) > 3`.

**3. Subqueries vs JOINs**
Both can answer the same question. A JOIN is generally more readable and often faster. Subqueries with `IN (SELECT ...)` are clearer when the inner set is small.

**4. Always use parameterized queries**
Never build SQL strings with f-strings containing user input — this allows SQL injection. Use `cursor.execute("WHERE gene = ?", (gene_name,))` instead.

```python
import sqlite3
import pandas as pd
import numpy as np

# Create an in-memory SQLite database for all examples
conn = sqlite3.connect(":memory:")
cursor = conn.cursor()

# Schema ----
cursor.executescript("""
CREATE TABLE genes (
    gene_id   INTEGER PRIMARY KEY,
    symbol    TEXT NOT NULL,
    chromosome TEXT,
    start_pos  INTEGER,
    end_pos    INTEGER,
    biotype   TEXT
);

CREATE TABLE variants (
    variant_id  INTEGER PRIMARY KEY,
    gene_id     INTEGER REFERENCES genes(gene_id),
    position    INTEGER,
    ref_allele  TEXT,
    alt_allele  TEXT,
    clinical_significance TEXT
);

CREATE TABLE expression (
    expr_id   INTEGER PRIMARY KEY,
    gene_id   INTEGER REFERENCES genes(gene_id),
    tissue    TEXT,
    condition TEXT,
    tpm       REAL
);

CREATE TABLE pathways (
    pathway_id   INTEGER PRIMARY KEY,
    pathway_name TEXT
);

CREATE TABLE gene_pathway (
    gene_id    INTEGER REFERENCES genes(gene_id),
    pathway_id INTEGER REFERENCES pathways(pathway_id)
);
""")

# Seed data ----
genes_data = [
    (1, 'BRCA1', 'chr17', 43044295, 43125483, 'protein_coding'),
    (2, 'TP53',  'chr17',  7661779,  7687538, 'protein_coding'),
    (3, 'EGFR',  'chr7',  55019017, 55207337, 'protein_coding'),
    (4, 'MYC',   'chr8', 127735434,127742951, 'protein_coding'),
    (5, 'KRAS',  'chr12', 25204789, 25250936, 'protein_coding'),
    (6, 'PTEN',  'chr10', 89692905, 89728532, 'protein_coding'),
    (7, 'RB1',   'chr13', 47775885, 47954065, 'protein_coding'),
]
cursor.executemany("INSERT INTO genes VALUES (?,?,?,?,?,?)", genes_data)

variants_data = [
    (1, 1, 43045629, 'A', 'T', 'pathogenic'),
    (2, 2, 7674220,  'C', 'T', 'pathogenic'),
    (3, 3, 55181320, 'G', 'A', 'likely_pathogenic'),
    (4, 5, 25245347, 'G', 'T', 'pathogenic'),
    (5, 6, 89711933, 'T', 'A', 'benign'),
]
cursor.executemany("INSERT INTO variants VALUES (?,?,?,?,?,?)", variants_data)

np.random.seed(42)
expr_rows = []
eid = 1
for gid in range(1, 8):
    for tissue in ['liver', 'kidney', 'brain']:
        for cond in ['normal', 'tumor']:
            base = np.random.uniform(10, 200)
            tpm = round(base * (1.8 if cond == 'tumor' else 1.0) + np.random.normal(0, 5), 2)
            expr_rows.append((eid, gid, tissue, cond, max(tpm, 0.1)))
            eid += 1
cursor.executemany("INSERT INTO expression VALUES (?,?,?,?,?)", expr_rows)

pathways_data = [(1,'DNA repair'), (2,'Cell cycle'), (3,'Apoptosis')]
cursor.executemany("INSERT INTO pathways VALUES (?,?)", pathways_data)

gp_data = [(1,1),(2,1),(2,2),(4,2),(3,3),(2,3)]
cursor.executemany("INSERT INTO gene_pathway VALUES (?,?)", gp_data)

conn.commit()
print("Database ready:",
      len(pd.read_sql_query("SELECT * FROM genes", conn)), "genes,",
      len(pd.read_sql_query("SELECT * FROM variants", conn)), "variants,",
      len(pd.read_sql_query("SELECT * FROM expression", conn)), "expression rows")
```python

## Basic Queries: SELECT, WHERE, ORDER BY

The core of SQL is `SELECT ... FROM ... WHERE`. Use `ORDER BY` to sort results and `LIMIT` to cap the number of rows returned.

```python
# All genes on chromosome 17
df = pd.read_sql_query(
    "SELECT symbol, chromosome, start_pos, end_pos FROM genes WHERE chromosome = 'chr17'",
    conn
)
print(df)

# Genes longer than 100 kb, ordered by length
df = pd.read_sql_query("""
    SELECT symbol, chromosome, (end_pos - start_pos) AS length
    FROM genes
    WHERE (end_pos - start_pos) > 100000
    ORDER BY length DESC
""", conn)
print(df)
```python

## Aggregate Functions and GROUP BY

`COUNT`, `AVG`, `MAX`, `MIN`, and `SUM` summarize groups of rows. `GROUP BY` splits the table into groups before aggregation. `HAVING` filters groups after aggregation (equivalent to `WHERE` on the aggregated result).

```python
# Average expression per tissue and condition
df = pd.read_sql_query("""
    SELECT tissue, condition, ROUND(AVG(tpm), 2) AS avg_tpm, COUNT(*) AS n_samples
    FROM expression
    GROUP BY tissue, condition
    ORDER BY tissue, condition
""", conn)
print(df)

# Genes with highest average tumor expression (HAVING filters after grouping)
df = pd.read_sql_query("""
    SELECT g.symbol, ROUND(AVG(e.tpm), 2) AS avg_tumor_tpm
    FROM genes g
    JOIN expression e ON g.gene_id = e.gene_id
    WHERE e.condition = 'tumor'
    GROUP BY g.symbol
    HAVING AVG(e.tpm) > 50
    ORDER BY avg_tumor_tpm DESC
""", conn)
print(df)
```python

## JOIN Operations

JOINs combine rows from two tables based on a related column.

- **INNER JOIN** — only rows that have a match in both tables.
- **LEFT JOIN** — all rows from the left table; NULLs where no match exists in the right table.
- **Self-join** — a table joined to itself, useful for comparing rows within the same table.

```python
# INNER JOIN: pathogenic variants with gene names
df = pd.read_sql_query("""
    SELECT g.symbol, g.chromosome, v.position, v.ref_allele, v.alt_allele, v.clinical_significance
    FROM variants v
    INNER JOIN genes g ON v.gene_id = g.gene_id
    WHERE v.clinical_significance = 'pathogenic'
    ORDER BY g.symbol
""", conn)
print(df)

# LEFT JOIN: all genes with variant count (including genes with 0 variants)
df = pd.read_sql_query("""
    SELECT g.symbol, COUNT(v.variant_id) AS n_variants
    FROM genes g
    LEFT JOIN variants v ON g.gene_id = v.gene_id
    GROUP BY g.symbol
    ORDER BY n_variants DESC
""", conn)
print(df)
```python

## Subqueries

```python
# Genes that are both highly expressed in tumors AND carry pathogenic variants
df = pd.read_sql_query("""
    SELECT symbol FROM genes
    WHERE gene_id IN (
        SELECT gene_id FROM expression
        WHERE condition = 'tumor'
        GROUP BY gene_id
        HAVING AVG(tpm) > 50
    )
    AND gene_id IN (
        SELECT gene_id FROM variants
        WHERE clinical_significance = 'pathogenic'
    )
""", conn)
print('Highly expressed tumor genes with pathogenic variants:')
print(df)
```python

## CREATE, INSERT, UPDATE, DELETE

SQL is not read-only. You can create new tables to store analysis results, insert rows, update values, and delete records. This is useful for storing intermediate results from a pipeline directly in a database.

```python
# Add a new analysis results table
cursor.execute('''
    CREATE TABLE IF NOT EXISTS de_results (
        gene_id INTEGER REFERENCES genes(gene_id),
        log2fc REAL,
        pvalue REAL,
        padj REAL,
        significant INTEGER
    )
''')

# Insert differential expression results
de_data = [
    (1, 2.3,  0.001,  0.01,  1),
    (2, 1.8,  0.005,  0.03,  1),
    (3, 0.2,  0.45,   0.78,  0),
    (5, 3.1,  0.0001, 0.002, 1),
]
cursor.executemany('INSERT INTO de_results VALUES (?,?,?,?,?)', de_data)
conn.commit()

# Query the results joined back to gene names
df = pd.read_sql_query("""
    SELECT g.symbol, d.log2fc, d.padj,
           CASE WHEN d.significant = 1 THEN 'Yes' ELSE 'No' END AS significant
    FROM de_results d
    JOIN genes g ON d.gene_id = g.gene_id
    ORDER BY d.padj
""", conn)
print(df)
```python

# Your query here
df = pd.read_sql_query("""
    -- fill query
""", conn)
print(df)
```python

**Exercise 2** (★★) — Calculate the fold change (tumor AVG TPM / normal AVG TPM) for each gene in breast tissue. Return gene symbol and fold change, ordered by fold change descending.

```python
# Your query here
df = pd.read_sql_query("""
    -- fill query
""", conn)
print(df)
```python

**Exercise 3** (★★) — Find genes that have more than one pathogenic variant. Return gene symbol and the count of pathogenic variants, ordered by count descending.

```python
# Your query here
df = pd.read_sql_query("""
    -- fill query
""", conn)
print(df)

conn.close()
```python

## Pitfalls

- **Mutable default arguments**: Never use `def f(x=[])` — use `def f(x=None)` and set inside the function
- **Off-by-one errors**: Python ranges are half-open `[start, stop)` — bioinformatics coordinates are often 1-based
- **Deep vs shallow copy**: Nested data structures require `copy.deepcopy()` — `list.copy()` only copies the top level
