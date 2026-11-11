import marimo

__generated_with = "0.25.2"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo

    return (mo,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # Schema Evaluation

    ## Measurment metrics
    - Loci: Total number of loci that were evaluated.
    - Alleles: Total number of alleles.

    ## List of evaluation criteria:
    - Valid Alleles: Total number of valid alleles. An allele is considered valid if its sequence size is a multiple of 3, if it has a single start and stop codon, if has no in-frame stop codons, and if it contains no ambiguous bases.
    """)
    return


@app.cell
def unzip():
    print("Hello world")
    return


if __name__ == "__main__":
    app.run()
