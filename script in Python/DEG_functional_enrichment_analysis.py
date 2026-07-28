#!/usr/bin/env python3

"""
RNA-seq DEG Functional Enrichment Analysis
==========================================

Script:
    DEG_functional_enrichment_analysis.py

Description:
    This script reads the complete DESeq2 result table, identifies significant
    differentially expressed genes, separates upregulated and downregulated
    genes, performs GO and KEGG enrichment through the official g:Profiler API,
    and generates enrichment tables, bar plots and bubble plots.

Input:
    data/DESeq2_All_Genes_Results.csv

Required columns:
    Gene_ID
    log2FoldChange
    padj

DEG thresholds:
    padj < 0.05
    absolute log2FoldChange > 1

Enrichment groups:
    1. All significant DEGs
    2. Upregulated DEGs
    3. Downregulated DEGs

Enrichment sources:
    GO:BP
    GO:CC
    GO:MF
    KEGG

Output:
    results/functional_enrichment/
        figures/
        tables/
        Session_Info.txt
"""

from __future__ import annotations

import json
import math
import platform
import sys
import time
from pathlib import Path
from typing import Iterable, Optional

import Bio
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests


# =============================================================================
# 1. USER SETTINGS
# =============================================================================

INPUT_FILE = Path("data") / "DESeq2_All_Genes_Results.csv"

OUTPUT_DIRECTORY = Path("results") / "functional_enrichment"
FIGURE_DIRECTORY = OUTPUT_DIRECTORY / "figures"
TABLE_DIRECTORY = OUTPUT_DIRECTORY / "tables"

# g:Profiler organism code for Oryza sativa
GPROFILER_ORGANISM = "osativa"

# DEG thresholds
PADJ_THRESHOLD = 0.05
ABS_LOG2FC_THRESHOLD = 1.0

# Enrichment threshold
ENRICHMENT_FDR_THRESHOLD = 0.05

# Number of enriched terms shown in each plot
TOP_TERMS_TO_PLOT = 10

# Figure settings
PLOT_WIDTH = 13
PLOT_HEIGHT = 8
PLOT_DPI = 300
LABEL_WRAP_WIDTH = 48

# Official g:Profiler enrichment endpoint
GPROFILER_URL = (
    "https://biit.cs.ut.ee/gprofiler/api/gost/profile/"
)

ENRICHMENT_SOURCES = [
    "GO:BP",
    "GO:CC",
    "GO:MF",
    "KEGG",
]

SOURCE_FILE_LABELS = {
    "GO:BP": "GO_BP",
    "GO:CC": "GO_CC",
    "GO:MF": "GO_MF",
    "KEGG": "KEGG",
}


# =============================================================================
# 2. CREATE OUTPUT DIRECTORIES
# =============================================================================

for directory in [
    OUTPUT_DIRECTORY,
    FIGURE_DIRECTORY,
    TABLE_DIRECTORY,
]:
    directory.mkdir(
        parents=True,
        exist_ok=True,
    )


# =============================================================================
# 3. HELPER FUNCTIONS
# =============================================================================

def normalise_column_name(name: str) -> str:
    """
    Standardise a column name for flexible matching.
    """

    return (
        str(name)
        .strip()
        .lower()
        .replace(" ", "")
        .replace("_", "")
        .replace(".", "")
        .replace("-", "")
    )


def detect_column(
    dataframe: pd.DataFrame,
    possible_names: Iterable[str],
    description: str,
) -> str:
    """
    Identify a required column using case-insensitive matching.
    """

    existing_columns = {
        normalise_column_name(column): column
        for column in dataframe.columns
    }

    for candidate in possible_names:

        normalised_candidate = normalise_column_name(
            candidate
        )

        if normalised_candidate in existing_columns:

            return existing_columns[
                normalised_candidate
            ]

    available_columns = ", ".join(
        map(str, dataframe.columns)
    )

    accepted_columns = ", ".join(
        possible_names
    )

    raise ValueError(
        f"Unable to identify the {description} column.\n\n"
        f"Available columns:\n{available_columns}\n\n"
        f"Accepted names include:\n{accepted_columns}"
    )


def clean_gene_list(
    values: Iterable[object],
) -> list[str]:
    """
    Remove missing, empty and duplicate gene identifiers.
    """

    genes: list[str] = []
    observed: set[str] = set()

    for value in values:

        if pd.isna(value):
            continue

        gene = str(value).strip()

        if not gene:
            continue

        if gene.lower() == "nan":
            continue

        if gene not in observed:

            genes.append(gene)
            observed.add(gene)

    return genes


def wrap_label(
    text: str,
    width: int = LABEL_WRAP_WIDTH,
) -> str:
    """
    Wrap long GO or KEGG labels across multiple lines.
    """

    words = str(text).split()

    lines: list[str] = []
    current_line: list[str] = []
    current_length = 0

    for word in words:

        additional_length = (
            len(word) + 1
            if current_line
            else len(word)
        )

        if (
            current_line
            and current_length + additional_length > width
        ):

            lines.append(
                " ".join(current_line)
            )

            current_line = [word]
            current_length = len(word)

        else:

            current_line.append(word)
            current_length += additional_length

    if current_line:

        lines.append(
            " ".join(current_line)
        )

    return "\n".join(lines)


def calculate_minus_log10(
    value: float,
) -> float:
    """
    Calculate -log10 while protecting against zero p-values.
    """

    minimum_float = np.finfo(float).tiny

    return -math.log10(
        max(
            float(value),
            minimum_float,
        )
    )


def save_dataframe(
    dataframe: pd.DataFrame,
    filename: str,
) -> None:
    """
    Save a DataFrame as a CSV file.
    """

    output_path = TABLE_DIRECTORY / filename

    dataframe.to_csv(
        output_path,
        index=False,
    )

    print(
        f"Saved table: {output_path}"
    )


def save_figure(
    figure: plt.Figure,
    base_filename: str,
) -> None:
    """
    Save a figure in PNG, PDF and SVG formats.
    """

    for extension in [
        "png",
        "pdf",
        "svg",
    ]:

        output_path = (
            FIGURE_DIRECTORY
            / f"{base_filename}.{extension}"
        )

        save_options = {
            "bbox_inches": "tight",
            "facecolor": "white",
        }

        if extension == "png":

            save_options["dpi"] = PLOT_DPI

        figure.savefig(
            output_path,
            **save_options,
        )

        print(
            f"Saved figure: {output_path}"
        )

    plt.close(figure)


# =============================================================================
# 4. READ DESEQ2 RESULTS
# =============================================================================

if not INPUT_FILE.exists():

    raise FileNotFoundError(
        f"Input file not found:\n"
        f"{INPUT_FILE.resolve()}\n\n"
        "Place DESeq2_All_Genes_Results.csv "
        "inside the data directory."
    )


print(
    f"Reading DESeq2 results: {INPUT_FILE}"
)


deseq2_results = pd.read_csv(
    INPUT_FILE
)


if deseq2_results.empty:

    raise ValueError(
        "The DESeq2 input file contains no rows."
    )


# =============================================================================
# 5. DETECT REQUIRED COLUMNS
# =============================================================================

gene_id_column = detect_column(
    dataframe=deseq2_results,
    possible_names=[
        "Gene_ID",
        "GeneID",
        "Gene",
        "gene_id",
        "geneid",
        "gene",
        "Feature",
        "Gene_Feature",
        "Row.names",
        "X",
    ],
    description="gene identifier",
)


log2fc_column = detect_column(
    dataframe=deseq2_results,
    possible_names=[
        "log2FoldChange",
        "log2FC",
        "logFC",
        "LFC",
        "log2_fold_change",
    ],
    description="log2 fold-change",
)


padj_column = detect_column(
    dataframe=deseq2_results,
    possible_names=[
        "padj",
        "FDR",
        "Adjusted_P_Value",
        "Adjusted.P.Value",
        "adj.P.Val",
        "qvalue",
        "q_value",
    ],
    description="adjusted p-value",
)


print(
    f"Detected gene ID column: {gene_id_column}"
)

print(
    f"Detected log2FC column: {log2fc_column}"
)

print(
    f"Detected padj column: {padj_column}"
)


# =============================================================================
# 6. STANDARDISE DESEQ2 COLUMNS
# =============================================================================

deseq2_results = deseq2_results.rename(
    columns={
        gene_id_column: "Gene_ID",
        log2fc_column: "log2FoldChange",
        padj_column: "padj",
    }
).copy()


deseq2_results["Gene_ID"] = (
    deseq2_results["Gene_ID"]
    .astype("string")
    .str.strip()
)


deseq2_results["log2FoldChange"] = (
    pd.to_numeric(
        deseq2_results["log2FoldChange"],
        errors="coerce",
    )
)


deseq2_results["padj"] = (
    pd.to_numeric(
        deseq2_results["padj"],
        errors="coerce",
    )
)


deseq2_results = deseq2_results[
    deseq2_results["Gene_ID"].notna()
    & deseq2_results["Gene_ID"].ne("")
].copy()


# Retain one result for each gene.
# If duplicate gene IDs occur, retain the row with the smallest padj.

deseq2_results = (
    deseq2_results
    .sort_values(
        ["Gene_ID", "padj"],
        na_position="last",
    )
    .drop_duplicates(
        subset="Gene_ID",
        keep="first",
    )
    .reset_index(drop=True)
)


# =============================================================================
# 7. CREATE BACKGROUND GENE LIST
# =============================================================================

background_genes = clean_gene_list(
    deseq2_results["Gene_ID"]
)


if not background_genes:

    raise ValueError(
        "No valid genes were available for the "
        "enrichment background."
    )


# =============================================================================
# 8. IDENTIFY SIGNIFICANT DEGS
# =============================================================================

significant_mask = (

    deseq2_results["padj"].notna()

    & deseq2_results[
        "log2FoldChange"
    ].notna()

    & (
        deseq2_results["padj"]
        < PADJ_THRESHOLD
    )

    & (
        deseq2_results[
            "log2FoldChange"
        ].abs()
        > ABS_LOG2FC_THRESHOLD
    )
)


significant_degs = (
    deseq2_results.loc[
        significant_mask
    ].copy()
)


significant_degs["Regulation"] = np.where(

    significant_degs[
        "log2FoldChange"
    ] > ABS_LOG2FC_THRESHOLD,

    "Upregulated",

    "Downregulated",
)


upregulated_degs = significant_degs[

    significant_degs[
        "log2FoldChange"
    ] > ABS_LOG2FC_THRESHOLD

].copy()


downregulated_degs = significant_degs[

    significant_degs[
        "log2FoldChange"
    ] < -ABS_LOG2FC_THRESHOLD

].copy()


# =============================================================================
# 9. EXPORT DEG TABLES
# =============================================================================

deg_summary = pd.DataFrame(
    {
        "Category": [
            "All DESeq2-tested genes",
            "Significant DEGs",
            "Upregulated significant DEGs",
            "Downregulated significant DEGs",
        ],

        "Number_of_Genes": [
            len(background_genes),
            len(significant_degs),
            len(upregulated_degs),
            len(downregulated_degs),
        ],

        "Criteria": [
            "All genes in complete DESeq2 results",

            (
                f"padj < {PADJ_THRESHOLD} and "
                f"|log2FoldChange| > "
                f"{ABS_LOG2FC_THRESHOLD}"
            ),

            (
                f"padj < {PADJ_THRESHOLD} and "
                f"log2FoldChange > "
                f"{ABS_LOG2FC_THRESHOLD}"
            ),

            (
                f"padj < {PADJ_THRESHOLD} and "
                f"log2FoldChange < "
                f"-{ABS_LOG2FC_THRESHOLD}"
            ),
        ],
    }
)


save_dataframe(
    deg_summary,
    "DEG_Summary.csv",
)


save_dataframe(
    deseq2_results,
    "All_DESeq2_Tested_Genes.csv",
)


save_dataframe(
    significant_degs,
    "All_Significant_DEGs.csv",
)


save_dataframe(
    upregulated_degs,
    "Upregulated_Significant_DEGs.csv",
)


save_dataframe(
    downregulated_degs,
    "Downregulated_Significant_DEGs.csv",
)


print(
    f"Background genes: {len(background_genes)}"
)

print(
    f"Significant DEGs: {len(significant_degs)}"
)

print(
    f"Upregulated DEGs: {len(upregulated_degs)}"
)

print(
    f"Downregulated DEGs: {len(downregulated_degs)}"
)


# =============================================================================
# 10. RUN G:PROFILER ENRICHMENT
# =============================================================================

def run_gprofiler_enrichment(
    genes: Iterable[object],
    analysis_name: str,
    background: Iterable[object],
    maximum_retries: int = 3,
) -> Optional[pd.DataFrame]:
    """
    Perform GO and KEGG enrichment using the official g:Profiler API.
    """

    query_genes = clean_gene_list(
        genes
    )

    background_gene_list = clean_gene_list(
        background
    )


    if not query_genes:

        print(
            f"No genes available for {analysis_name}. "
            "Enrichment analysis was skipped."
        )

        return None


    request_payload = {

        "organism": GPROFILER_ORGANISM,

        "query": query_genes,

        "sources": ENRICHMENT_SOURCES,

        "user_threshold": (
            ENRICHMENT_FDR_THRESHOLD
        ),

        "significance_threshold_method": "fdr",

        "all_results": False,

        "ordered": False,

        "combined": False,

        "measure_underrepresentation": False,

        "no_iea": False,

        "no_evidences": False,

        "domain_scope": "custom_annotated",

        "background": background_gene_list,
    }


    request_headers = {

        "Content-Type": "application/json",

        "Accept": "application/json",

        "User-Agent": (
            "RNAseq-Functional-Enrichment-Analysis/1.0"
        ),
    }


    print(
        f"\nRunning enrichment for {analysis_name} "
        f"using {len(query_genes)} genes..."
    )


    response: Optional[
        requests.Response
    ] = None


    for attempt in range(
        1,
        maximum_retries + 1,
    ):

        try:

            response = requests.post(

                GPROFILER_URL,

                json=request_payload,

                headers=request_headers,

                timeout=180,
            )


            response.raise_for_status()

            break


        except requests.RequestException as error:

            if attempt == maximum_retries:

                print(
                    f"g:Profiler failed for "
                    f"{analysis_name}: {error}"
                )

                return None


            waiting_time = attempt * 3

            print(
                f"Attempt {attempt} failed. "
                f"Retrying in {waiting_time} seconds..."
            )

            time.sleep(
                waiting_time
            )


    if response is None:

        return None


    response_json = response.json()

    result_records = response_json.get(
        "result",
        [],
    )


    if not result_records:

        print(
            f"No significant enriched terms "
            f"were found for {analysis_name}."
        )

        return None


    enrichment = pd.DataFrame(
        result_records
    )


    required_columns = {

        "source",

        "native",

        "name",

        "p_value",

        "term_size",

        "query_size",

        "intersection_size",

        "effective_domain_size",
    }


    missing_columns = required_columns.difference(
        enrichment.columns
    )


    if missing_columns:

        raise ValueError(
            f"g:Profiler output is missing columns: "
            f"{', '.join(sorted(missing_columns))}"
        )


    enrichment["Analysis"] = analysis_name

    enrichment["GO_KEGG_ID"] = (
        enrichment["native"].astype(str)
    )

    enrichment["Term_Name"] = (
        enrichment["name"].astype(str)
    )

    enrichment["Source"] = (
        enrichment["source"].astype(str)
    )


    enrichment["Adjusted_P_Value"] = (
        pd.to_numeric(
            enrichment["p_value"],
            errors="coerce",
        )
    )


    enrichment[
        "Minus_Log10_Adjusted_P"
    ] = enrichment[
        "Adjusted_P_Value"
    ].map(
        calculate_minus_log10
    )


    enrichment["Gene_Count"] = (
        pd.to_numeric(
            enrichment["intersection_size"],
            errors="coerce",
        )
    )


    enrichment["Query_Size"] = (
        pd.to_numeric(
            enrichment["query_size"],
            errors="coerce",
        )
    )


    enrichment["Term_Size"] = (
        pd.to_numeric(
            enrichment["term_size"],
            errors="coerce",
        )
    )


    enrichment[
        "Effective_Domain_Size"
    ] = pd.to_numeric(

        enrichment[
            "effective_domain_size"
        ],

        errors="coerce",
    )


    enrichment["Gene_Ratio"] = np.where(

        enrichment["Query_Size"] > 0,

        enrichment["Gene_Count"]
        / enrichment["Query_Size"],

        np.nan,
    )


    enrichment["Background_Ratio"] = np.where(

        enrichment[
            "Effective_Domain_Size"
        ] > 0,

        enrichment["Term_Size"]
        / enrichment[
            "Effective_Domain_Size"
        ],

        np.nan,
    )


    enrichment["Plot_Label"] = (

        enrichment["GO_KEGG_ID"]

        + " - "

        + enrichment["Term_Name"]
    )


    if "intersections" in enrichment.columns:

        enrichment["Genes_in_Term"] = (

            enrichment["intersections"]
            .apply(
                lambda value:
                "; ".join(map(str, value))

                if isinstance(value, list)

                else str(value)
            )
        )

    elif "intersection" in enrichment.columns:

        enrichment["Genes_in_Term"] = (

            enrichment["intersection"]
            .apply(
                lambda value:
                "; ".join(map(str, value))

                if isinstance(value, list)

                else str(value)
            )
        )

    else:

        enrichment["Genes_in_Term"] = pd.NA


    preferred_columns = [

        "Analysis",

        "GO_KEGG_ID",

        "Term_Name",

        "Source",

        "Gene_Count",

        "Query_Size",

        "Gene_Ratio",

        "Term_Size",

        "Background_Ratio",

        "Adjusted_P_Value",

        "Minus_Log10_Adjusted_P",

        "Genes_in_Term",
    ]


    remaining_columns = [

        column

        for column in enrichment.columns

        if column not in preferred_columns
    ]


    enrichment = enrichment[

        preferred_columns
        + remaining_columns

    ].sort_values(

        [
            "Source",
            "Adjusted_P_Value",
            "Gene_Count",
        ],

        ascending=[
            True,
            True,
            False,
        ],

    ).reset_index(
        drop=True
    )


    save_dataframe(

        enrichment,

        (
            f"{analysis_name}_"
            "Complete_Functional_Enrichment.csv"
        ),
    )


    for source in ENRICHMENT_SOURCES:

        source_table = enrichment[

            enrichment["Source"] == source

        ].copy()


        if not source_table.empty:

            save_dataframe(

                source_table,

                (
                    f"{analysis_name}_"
                    f"{SOURCE_FILE_LABELS[source]}_"
                    "Enrichment.csv"
                ),
            )


    metadata = response_json.get(
        "meta",
        {},
    )


    metadata_file = (

        TABLE_DIRECTORY

        / (
            f"{analysis_name}_"
            "gProfiler_Metadata.json"
        )
    )


    with metadata_file.open(
        "w",
        encoding="utf-8",
    ) as output_handle:

        json.dump(
            metadata,
            output_handle,
            indent=2,
            default=str,
        )


    print(
        f"Saved metadata: {metadata_file}"
    )


    return enrichment


# =============================================================================
# 11. RUN ENRICHMENT FOR THREE DEG GROUPS
# =============================================================================

all_significant_enrichment = (
    run_gprofiler_enrichment(

        genes=significant_degs["Gene_ID"],

        analysis_name=(
            "All_Significant_DEGs"
        ),

        background=background_genes,
    )
)


upregulated_enrichment = (
    run_gprofiler_enrichment(

        genes=upregulated_degs["Gene_ID"],

        analysis_name=(
            "Upregulated_DEGs"
        ),

        background=background_genes,
    )
)


downregulated_enrichment = (
    run_gprofiler_enrichment(

        genes=downregulated_degs["Gene_ID"],

        analysis_name=(
            "Downregulated_DEGs"
        ),

        background=background_genes,
    )
)


# =============================================================================
# 12. SELECT TOP ENRICHED TERMS
# =============================================================================

def select_top_terms(
    enrichment: Optional[pd.DataFrame],
    source: str,
    top_n: int = TOP_TERMS_TO_PLOT,
) -> Optional[pd.DataFrame]:
    """
    Select the most significant enriched terms from one database.
    """

    if enrichment is None:

        return None


    if enrichment.empty:

        return None


    selected = enrichment[

        (enrichment["Source"] == source)

        & enrichment[
            "Adjusted_P_Value"
        ].notna()

        & (
            enrichment[
                "Adjusted_P_Value"
            ]
            <= ENRICHMENT_FDR_THRESHOLD
        )

    ].copy()


    if selected.empty:

        return None


    selected = (

        selected

        .sort_values(

            [
                "Adjusted_P_Value",
                "Gene_Count",
            ],

            ascending=[
                True,
                False,
            ],
        )

        .head(top_n)

        .copy()
    )


    selected["Wrapped_Label"] = (

        selected["Plot_Label"]
        .map(wrap_label)
    )


    # Reverse order so the most significant term is displayed at the top.

    selected = (

        selected

        .iloc[::-1]

        .reset_index(
            drop=True
        )
    )


    return selected


# =============================================================================
# 13. BAR PLOT FUNCTION
# =============================================================================

def create_barplot(
    enrichment: Optional[pd.DataFrame],
    analysis_name: str,
    source: str,
    top_n: int = TOP_TERMS_TO_PLOT,
) -> None:
    """
    Generate a horizontal enrichment bar plot.
    """

    plot_data = select_top_terms(

        enrichment=enrichment,

        source=source,

        top_n=top_n,
    )


    if plot_data is None:

        print(
            f"No {source} terms available "
            f"for the {analysis_name} bar plot."
        )

        return


    figure, axis = plt.subplots(

        figsize=(
            PLOT_WIDTH,
            PLOT_HEIGHT,
        )
    )


    significance_values = plot_data[

        "Minus_Log10_Adjusted_P"

    ].to_numpy()


    adjusted_p_values = plot_data[

        "Adjusted_P_Value"

    ].to_numpy()


    y_positions = np.arange(
        len(plot_data)
    )


    colour_normalisation = (

        matplotlib.colors.LogNorm(

            vmin=max(

                adjusted_p_values.min(),

                np.finfo(float).tiny,
            ),

            vmax=max(

                adjusted_p_values.max(),

                adjusted_p_values.min()
                * 1.001,
            ),
        )
    )


    colour_map = plt.get_cmap(
        "viridis_r"
    )


    axis.barh(

        y_positions,

        significance_values,

        color=colour_map(

            colour_normalisation(
                adjusted_p_values
            )
        ),

        height=0.72,
    )


    axis.set_yticks(
        y_positions
    )


    axis.set_yticklabels(

        plot_data["Wrapped_Label"],

        fontsize=10,
    )


    axis.set_xlabel(

        r"$-\log_{10}$(adjusted p-value)",

        fontsize=13,

        fontweight="bold",
    )


    axis.set_ylabel(

        "Enriched terms",

        fontsize=13,

        fontweight="bold",
    )


    axis.set_title(

        f"Top {len(plot_data)} Enriched Terms\n"
        f"{source} ({analysis_name})",

        fontsize=16,

        fontweight="bold",

        pad=14,
    )


    axis.spines["top"].set_visible(
        False
    )

    axis.spines["right"].set_visible(
        False
    )


    colour_object = (

        matplotlib.cm.ScalarMappable(

            norm=colour_normalisation,

            cmap=colour_map,
        )
    )


    colour_object.set_array([])


    colour_bar = figure.colorbar(

        colour_object,

        ax=axis,

        pad=0.02,
    )


    colour_bar.set_label(

        "Adjusted p-value",

        fontsize=11,

        fontweight="bold",
    )


    figure.tight_layout()


    base_filename = (

        f"{analysis_name}_"

        f"{SOURCE_FILE_LABELS[source]}_"

        f"Top{top_n}_Barplot"
    )


    save_figure(

        figure,

        base_filename,
    )


# =============================================================================
# 14. BUBBLE PLOT FUNCTION
# =============================================================================

def create_bubbleplot(
    enrichment: Optional[pd.DataFrame],
    analysis_name: str,
    source: str,
    top_n: int = TOP_TERMS_TO_PLOT,
) -> None:
    """
    Generate an enrichment bubble plot.

    X-axis:
        -log10 adjusted p-value

    Y-axis:
        GO or KEGG terms

    Bubble size:
        Gene count

    Bubble colour:
        Adjusted p-value
    """

    plot_data = select_top_terms(

        enrichment=enrichment,

        source=source,

        top_n=top_n,
    )


    if plot_data is None:

        print(
            f"No {source} terms available "
            f"for the {analysis_name} bubble plot."
        )

        return


    figure, axis = plt.subplots(

        figsize=(
            PLOT_WIDTH,
            PLOT_HEIGHT,
        )
    )


    x_values = plot_data[

        "Minus_Log10_Adjusted_P"

    ].to_numpy()


    adjusted_p_values = plot_data[

        "Adjusted_P_Value"

    ].to_numpy()


    gene_counts = plot_data[

        "Gene_Count"

    ].to_numpy()


    y_positions = np.arange(
        len(plot_data)
    )


    minimum_bubble_size = 70
    maximum_bubble_size = 500


    if gene_counts.max() == gene_counts.min():

        bubble_sizes = np.full_like(

            gene_counts,

            fill_value=(

                minimum_bubble_size
                + maximum_bubble_size

            ) / 2,

            dtype=float,
        )

    else:

        bubble_sizes = (

            minimum_bubble_size

            + (

                (
                    gene_counts
                    - gene_counts.min()
                )

                / (
                    gene_counts.max()
                    - gene_counts.min()
                )
            )

            * (
                maximum_bubble_size
                - minimum_bubble_size
            )
        )


    colour_normalisation = (

        matplotlib.colors.LogNorm(

            vmin=max(

                adjusted_p_values.min(),

                np.finfo(float).tiny,
            ),

            vmax=max(

                adjusted_p_values.max(),

                adjusted_p_values.min()
                * 1.001,
            ),
        )
    )


    scatter = axis.scatter(

        x_values,

        y_positions,

        s=bubble_sizes,

        c=adjusted_p_values,

        cmap="viridis_r",

        norm=colour_normalisation,

        alpha=0.9,

        edgecolors="black",

        linewidths=0.4,
    )


    axis.set_yticks(
        y_positions
    )


    axis.set_yticklabels(

        plot_data["Wrapped_Label"],

        fontsize=10,
    )


    axis.set_xlabel(

        r"$-\log_{10}$(adjusted p-value)",

        fontsize=13,

        fontweight="bold",
    )


    axis.set_ylabel(

        "Enriched terms",

        fontsize=13,

        fontweight="bold",
    )


    axis.set_title(

        f"Top {len(plot_data)} Enriched Terms\n"
        f"{source} ({analysis_name})",

        fontsize=16,

        fontweight="bold",

        pad=14,
    )


    axis.spines["top"].set_visible(
        False
    )

    axis.spines["right"].set_visible(
        False
    )


    colour_bar = figure.colorbar(

        scatter,

        ax=axis,

        pad=0.02,
    )


    colour_bar.set_label(

        "Adjusted p-value",

        fontsize=11,

        fontweight="bold",
    )


    legend_counts = np.unique(

        np.round(

            np.linspace(

                gene_counts.min(),

                gene_counts.max(),

                4,
            )

        ).astype(int)
    )


    legend_handles = []


    for count in legend_counts:

        if gene_counts.max() == gene_counts.min():

            size = (

                minimum_bubble_size
                + maximum_bubble_size

            ) / 2

        else:

            size = (

                minimum_bubble_size

                + (

                    (
                        count
                        - gene_counts.min()
                    )

                    / (
                        gene_counts.max()
                        - gene_counts.min()
                    )
                )

                * (
                    maximum_bubble_size
                    - minimum_bubble_size
                )
            )


        legend_handle = axis.scatter(

            [],

            [],

            s=size,

            facecolors="none",

            edgecolors="black",

            label=str(int(count)),
        )


        legend_handles.append(
            legend_handle
        )


    axis.legend(

        handles=legend_handles,

        title="Gene count",

        loc="lower right",

        frameon=False,

        fontsize=9,

        title_fontsize=10,
    )


    figure.tight_layout()


    base_filename = (

        f"{analysis_name}_"

        f"{SOURCE_FILE_LABELS[source]}_"

        f"Top{top_n}_Bubbleplot"
    )


    save_figure(

        figure,

        base_filename,
    )


# =============================================================================
# 15. GENERATE TABLES AND FIGURES
# =============================================================================

def process_enrichment_group(
    enrichment: Optional[pd.DataFrame],
    analysis_name: str,
) -> None:
    """
    Generate all GO and KEGG plots for one DEG group.
    """

    if enrichment is None:

        print(
            f"No enrichment results available "
            f"for {analysis_name}."
        )

        return


    if enrichment.empty:

        print(
            f"No enrichment results available "
            f"for {analysis_name}."
        )

        return


    for source in ENRICHMENT_SOURCES:

        top_terms = select_top_terms(

            enrichment=enrichment,

            source=source,

            top_n=TOP_TERMS_TO_PLOT,
        )


        if top_terms is not None:

            save_dataframe(

                top_terms,

                (
                    f"{analysis_name}_"

                    f"{SOURCE_FILE_LABELS[source]}_"

                    f"Top{TOP_TERMS_TO_PLOT}_Terms.csv"
                ),
            )


        create_barplot(

            enrichment=enrichment,

            analysis_name=analysis_name,

            source=source,

            top_n=TOP_TERMS_TO_PLOT,
        )


        create_bubbleplot(

            enrichment=enrichment,

            analysis_name=analysis_name,

            source=source,

            top_n=TOP_TERMS_TO_PLOT,
        )


process_enrichment_group(

    enrichment=all_significant_enrichment,

    analysis_name="All_Significant_DEGs",
)


process_enrichment_group(

    enrichment=upregulated_enrichment,

    analysis_name="Upregulated_DEGs",
)


process_enrichment_group(

    enrichment=downregulated_enrichment,

    analysis_name="Downregulated_DEGs",
)


# =============================================================================
# 16. SAVE ANALYSIS SETTINGS
# =============================================================================

analysis_settings = pd.DataFrame(
    {
        "Setting": [

            "Input file",

            "Organism",

            "g:Profiler organism code",

            "DEG adjusted p-value threshold",

            "DEG absolute log2FC threshold",

            "Enrichment FDR threshold",

            "Number of terms per plot",

            "Background definition",

            "Enrichment databases",

            "Plot formats",

            "Plot resolution",
        ],

        "Value": [

            str(
                INPUT_FILE.resolve()
            ),

            "Oryza sativa",

            GPROFILER_ORGANISM,

            PADJ_THRESHOLD,

            ABS_LOG2FC_THRESHOLD,

            ENRICHMENT_FDR_THRESHOLD,

            TOP_TERMS_TO_PLOT,

            (
                "All genes tested in the "
                "complete DESeq2 results"
            ),

            "; ".join(
                ENRICHMENT_SOURCES
            ),

            "PNG; PDF; SVG",

            f"{PLOT_DPI} dpi",
        ],
    }
)


save_dataframe(

    analysis_settings,

    "Analysis_Settings.csv",
)


# =============================================================================
# 17. SAVE SESSION INFORMATION
# =============================================================================

session_information_file = (

    OUTPUT_DIRECTORY

    / "Session_Info.txt"
)


with session_information_file.open(

    "w",

    encoding="utf-8",

) as output_handle:

    output_handle.write(
        "RNA-seq Functional Enrichment Analysis\n"
    )

    output_handle.write(
        "======================================\n\n"
    )

    output_handle.write(
        f"Python version: {sys.version}\n"
    )

    output_handle.write(
        f"Platform: {platform.platform()}\n"
    )

    output_handle.write(
        f"Biopython version: {Bio.__version__}\n"
    )

    output_handle.write(
        f"pandas version: {pd.__version__}\n"
    )

    output_handle.write(
        f"NumPy version: {np.__version__}\n"
    )

    output_handle.write(
        f"matplotlib version: "
        f"{matplotlib.__version__}\n"
    )

    output_handle.write(
        f"requests version: "
        f"{requests.__version__}\n"
    )

    output_handle.write(
        f"g:Profiler endpoint: "
        f"{GPROFILER_URL}\n"
    )


print(
    f"Saved session information: "
    f"{session_information_file}"
)


# =============================================================================
# 18. FINAL SUMMARY
# =============================================================================

print(
    "\n"
    + "=" * 78
)

print(
    "FUNCTIONAL ENRICHMENT ANALYSIS COMPLETED"
)

print(
    "=" * 78
)


print(
    f"All DESeq2-tested genes: "
    f"{len(background_genes)}"
)

print(
    f"Significant DEGs: "
    f"{len(significant_degs)}"
)

print(
    f"Upregulated significant DEGs: "
    f"{len(upregulated_degs)}"
)

print(
    f"Downregulated significant DEGs: "
    f"{len(downregulated_degs)}"
)


print(
    "\nDEG thresholds:"
)

print(
    f"  padj < {PADJ_THRESHOLD}"
)

print(
    f"  |log2FoldChange| > "
    f"{ABS_LOG2FC_THRESHOLD}"
)


print(
    "\nEnrichment databases:"
)


for source in ENRICHMENT_SOURCES:

    print(
        f"  {source}"
    )


print(
    f"\nFigures saved to:\n"
    f"  {FIGURE_DIRECTORY.resolve()}"
)

print(
    f"\nTables saved to:\n"
    f"  {TABLE_DIRECTORY.resolve()}"
)

print(
    f"\nSession information saved to:\n"
    f"  {session_information_file.resolve()}"
)

print(
    "=" * 78
)
