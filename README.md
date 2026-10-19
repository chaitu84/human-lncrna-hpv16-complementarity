Theoretical Therapeutic Discovery: HPV RNA Analysis
Overview

This is a computational biology project looking at whether human non-coding RNA sequences could point toward future RNA-based ways of targeting high-risk HPV.

Right now we are focused on HPV16 and HPV18, specifically the E6 and E7 oncogenes, which involves with human tumor-suppressor pathways.

The project is still in progress. We haven’t locked in the final approach, tools, or implementation details yet.

Problem

High-risk HPV types express E6 and E7, and those genes play a role in how cervical and other HPV-driven cancers develop. The question we want to look at is whether comparing public human non-coding RNA data against HPV reference sequences turns up any sequence relationships worth following up on.

Goal

Build a Python tool that compares public HPV reference sequences to human non-coding RNA datasets and summarizes the matches and predicted RNA interactions that come out of it.

What the Final Product Could Look Like

A tool that can let someone explore and compare possible RNA-based strategies for targeting HPV oncogene transcripts.

Features that we are considering:

Comparing HPV16 and HPV18 E6/E7 reference sequences
Pulling in human non-coding RNA datasets from public databases
Finding and summarizing sequence matches that look relevant
Running RNA structure or interaction predictions, depending on which tools I go with
Putting the results into a readable report or a simple interactive visualization

Data Sources
NCBI Nucleotide or PaVE for HPV reference sequences
GENCODE or Ensembl for human non-coding RNA annotations

Tools Required

Python is the main language I’m planning on. Beyond that:

Biopython for handling sequences
pandas for organizing and analyzing the data
BLAST+ or something similar for sequence comparison
ViennaRNA for RNA structure and interaction predictions
Some kind of visualization tool or lightweight dashboard if it helps show the results

Scope and Limitations

A sequence match or a predicted RNA interaction doesn’t mean a molecule will actually suppress HPV gene expression, or that it’d be safe, or that it’d work as a treatment.

The project only looks at computational evidence. Experimental validation, clinical effectiveness, and anything involving delivery through an engineered viral vector are all outside what this software can show.

Current Status

Midterm — proposal and planning.

Team

Jacky - Computer Science - Coder for the project
Arnav - Data Science - Researcher for the project
Cleo - Computer Science - Coder
Ekijah - Mathematics & Economics - Researcher for the project
