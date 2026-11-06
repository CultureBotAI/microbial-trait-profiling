"""Cross-disease comparison of decreased cellulose-annotated taxa (#550, #552, #557, #1127).

:mod:`mtp.cellulose_decr.analysis` writes the tables, summary JSON and provenance
in ``data/analysis/cellulose_decr/`` (``just cellulose-decr``).
:mod:`mtp.cellulose_decr.figures` draws the current figures from those tables
(``just cellulose-decr-figures``).

The October 16 historical Tufte plates (``figures_tufte/``) and 16-slide deck
(``cellulose_decr_slide_deck.html``) are retired because they contain unverified
numbers (#1572). They are excluded from default generation.
:mod:`mtp.cellulose_decr.tufte_plates` remains available for audit only via
``just cellulose-decr-tufte --audit-retired``. Every historical plate, standalone
SVG and HTML slide carries ``RETIRED: historical, unverified numbers (#1572)``.
A replacement built from current rerun data is pending under #1670.
methods/cellulose_decr.md describes the analysis and figure pipeline.
"""
