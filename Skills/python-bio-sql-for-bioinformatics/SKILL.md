---
name: python-bio-sql-for-bioinformatics
description: Ensembl, UCSC Genome Browser, NCBI, and dbSNP are all backed by relational databases. Even locally, SQLite is a practical way to store and query gene annotations, variant tables, and expression result
tool_type: python
primary_tool: Pandas
---

# SQL for Bioinformatics

## Key Concepts

**JOIN type determines which rows survive:** `INNER JOIN` — only matched rows; `LEFT JOIN` — all left rows, NULL for unmatched right side.

**`HAVING` vs `WHERE`:** `WHERE` filters before grouping; `HAVING` filters after. Use `HAVING COUNT(*) > 3` to filter aggregate groups.

**Subqueries vs JOINs:** JOIN is more readable and usually faster. Subqueries with `IN (SELECT ...)` are clearer for small inner sets.

**Always use parameterized queries:** Never build SQL with f-strings containing user input (SQL injection). Use `cursor.execute("WHERE gene = ?", (gene_name,))`.

```python
import sqlite3
import pandas as pd
import numpy as np

conn = sqlite3.connect(":memory:")
cursor = conn.cursor()

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
```

## SELECT, WHERE, ORDER BY

```python
# All genes on chromosome 17
df = pd.read_sql_query(
    "SELECT symbol, chromosome, start_pos, end_pos FROM genes WHERE chromosome = 'chr17'",
    conn
)

# Genes longer than 100 kb, ordered by length
df = pd.read_sql_query("""
    SELECT symbol, chromosome, (end_pos - start_pos) AS length
    FROM genes
    WHERE (end_pos - start_pos) > 100000
    ORDER BY length DESC
""", conn)
```

## Aggregate Functions and GROUP BY

```python
# Average expression per tissue and condition
df = pd.read_sql_query("""
    SELECT tissue, condition, ROUND(AVG(tpm), 2) AS avg_tpm, COUNT(*) AS n_samples
    FROM expression
    GROUP BY tissue, condition
    ORDER BY tissue, condition
""", conn)

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
```

## JOIN Operations

- **INNER JOIN** — only rows with a match in both tables.
- **LEFT JOIN** — all rows from left table; NULLs where no match on right.

```python
# INNER JOIN: pathogenic variants with gene names
df = pd.read_sql_query("""
    SELECT g.symbol, g.chromosome, v.position, v.ref_allele, v.alt_allele, v.clinical_significance
    FROM variants v
    INNER JOIN genes g ON v.gene_id = g.gene_id
    WHERE v.clinical_significance = 'pathogenic'
    ORDER BY g.symbol
""", conn)

# LEFT JOIN: all genes with variant count (including genes with 0 variants)
df = pd.read_sql_query("""
    SELECT g.symbol, COUNT(v.variant_id) AS n_variants
    FROM genes g
    LEFT JOIN variants v ON g.gene_id = v.gene_id
    GROUP BY g.symbol
    ORDER BY n_variants DESC
""", conn)
```

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
```

## CREATE, INSERT, UPDATE, DELETE

```python
cursor.execute('''
    CREATE TABLE IF NOT EXISTS de_results (
        gene_id INTEGER REFERENCES genes(gene_id),
        log2fc REAL,
        pvalue REAL,
        padj REAL,
        significant INTEGER
    )
''')

de_data = [
    (1, 2.3,  0.001,  0.01,  1),
    (2, 1.8,  0.005,  0.03,  1),
    (3, 0.2,  0.45,   0.78,  0),
    (5, 3.1,  0.0001, 0.002, 1),
]
cursor.executemany('INSERT INTO de_results VALUES (?,?,?,?,?)', de_data)
conn.commit()

df = pd.read_sql_query("""
    SELECT g.symbol, d.log2fc, d.padj,
           CASE WHEN d.significant = 1 THEN 'Yes' ELSE 'No' END AS significant
    FROM de_results d
    JOIN genes g ON d.gene_id = g.gene_id
    ORDER BY d.padj
""", conn)

conn.close()
```

## Pitfalls

- **JOIN type**: INNER loses unmatched rows; LEFT preserves them — choose deliberately.
- **HAVING vs WHERE**: Using WHERE on an aggregate silently fails or raises an error; use HAVING.
- **Parameterized queries**: f-string SQL with user input = SQL injection risk.
- **Off-by-one**: Python ranges are half-open `[start, stop)`; bioinformatics coordinates are often 1-based.
