"""TE data processing and cleanup."""

from __future__ import annotations

import logging
from typing import Optional

import pandas as pd

logger = logging.getLogger(__name__)

# Common TE family -> class mapping for filling missing te_class values
_FAMILY_TO_CLASS = {
    "Alu": "SINE", "AluJ": "SINE", "AluS": "SINE", "AluY": "SINE",
    "MIR": "SINE", "tRNA": "SINE",
    "L1": "LINE", "L1HS": "LINE", "L1PA2": "LINE", "L1PA3": "LINE",
    "L2": "LINE", "CR1": "LINE", "RTE": "LINE",
    "ERV1": "LTR", "ERVL": "LTR", "ERVL-MaLR": "LTR", "ERVK": "LTR",
    "HERV-K": "LTR", "Gypsy": "LTR", "Copia": "LTR",
    "SVA": "Retroposon", "SVA_A": "Retroposon", "SVA_B": "Retroposon",
    "SVA_C": "Retroposon", "SVA_D": "Retroposon", "SVA_E": "Retroposon",
    "SVA_F": "Retroposon",
    "hAT": "DNA", "hAT-Charlie": "DNA", "hAT-Tip100": "DNA",
    "TcMar": "DNA", "TcMar-Tigger": "DNA", "TcMar-Mariner": "DNA",
    "MULE": "DNA", "Helitron": "DNA", "CMC": "DNA",
    "ALR/Alpha": "Satellite", "HSATII": "Satellite", "centr": "Satellite",
}


def process_te_data(
    input_file: str,
    reference: Optional[str],
    output_file: str,
) -> None:
    """Process TE data: fill missing values, recalculate percentages, validate.

    Cleans a TE annotation CSV by:
    1. Merging with reference table to obtain seq_length if missing.
    2. Recalculating te_pct from te_bp / seq_length when possible.
    3. Filling missing te_class from te_family using a built-in mapping.
    4. Validating: no negative values, percentages clamped to 0-100.

    Args:
        input_file: TE annotation CSV file.
        reference: Optional reference table (CSV/TSV) with columns
            containing sequence name and length for merging.
        output_file: Output cleaned CSV file path.
    """
    df = pd.read_csv(input_file)
    logger.info(f"Loaded {len(df)} rows from {input_file}")
    initial_cols = list(df.columns)

    # --- Step 1: Merge with reference for seq_length ---
    if reference:
        ref_df = _load_reference(reference)
        if ref_df is not None and not ref_df.empty:
            # Find merge key
            merge_key = _find_merge_key(df, ref_df)
            if merge_key:
                df_key, ref_key = merge_key
                before = len(df)
                df = df.merge(
                    ref_df[[ref_key, "seq_length"]].drop_duplicates(),
                    left_on=df_key, right_on=ref_key, how="left",
                    suffixes=("", "_ref"),
                )
                # Use reference seq_length to fill missing
                if "eccDNA_len" in df.columns:
                    df["eccDNA_len"] = df["eccDNA_len"].fillna(df.get("seq_length"))
                elif "seq_length" not in initial_cols:
                    df.rename(columns={"seq_length": "eccDNA_len"}, inplace=True)
                # Clean up extra columns
                for col in ["seq_length_ref", ref_key + "_ref"]:
                    if col in df.columns and col not in initial_cols:
                        df.drop(columns=[col], inplace=True, errors="ignore")
                logger.info(f"Merged with reference ({merge_key}): {before} -> {len(df)} rows")
            else:
                logger.warning("Could not find matching merge key between input and reference")

    # --- Step 2: Recalculate te_pct ---
    len_col = _find_column(df, ["eccDNA_len", "seq_length", "length", "len"])
    bp_col = _find_column(df, ["te_bp", "overlap_bp", "te_coverage_bp"])

    if len_col and bp_col:
        df[bp_col] = pd.to_numeric(df[bp_col], errors="coerce").fillna(0)
        df[len_col] = pd.to_numeric(df[len_col], errors="coerce")
        valid = df[len_col] > 0
        df.loc[valid, "te_pct"] = (df.loc[valid, bp_col] / df.loc[valid, len_col] * 100)
        n_filled = valid.sum()
        logger.info(f"Recalculated te_pct for {n_filled} rows from {bp_col}/{len_col}")
    elif "te_pct" not in df.columns:
        logger.warning("Cannot compute te_pct: need (te_bp + eccDNA_len) or te_pct column")

    # --- Step 3: Fill missing te_class from te_family ---
    if "te_family" in df.columns and "te_class" in df.columns:
        missing_class = df["te_class"].isna() | (df["te_class"] == "") | (df["te_class"] == "Unknown")
        if missing_class.any():
            filled = 0
            for idx in df.index[missing_class]:
                family = str(df.at[idx, "te_family"])
                # Try exact match, then prefix match
                mapped = _FAMILY_TO_CLASS.get(family)
                if not mapped:
                    for prefix, cls in _FAMILY_TO_CLASS.items():
                        if family.startswith(prefix):
                            mapped = cls
                            break
                if mapped:
                    df.at[idx, "te_class"] = mapped
                    filled += 1
            logger.info(f"Filled {filled}/{missing_class.sum()} missing te_class values")

    # --- Step 4: Validate ---
    n_fixes = 0
    # No negative bp values
    for col in ["te_bp", "overlap_bp", "te_coverage_bp", "eccDNA_len", "seq_length"]:
        if col in df.columns:
            neg = df[col] < 0
            if neg.any():
                df.loc[neg, col] = 0
                n_fixes += neg.sum()

    # Clamp percentages to 0-100
    for col in ["te_pct", "class_pct", "family_pct"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
            out_of_range = (df[col] < 0) | (df[col] > 100)
            if out_of_range.any():
                df[col] = df[col].clip(0, 100)
                n_fixes += out_of_range.sum()

    if n_fixes > 0:
        logger.info(f"Fixed {n_fixes} out-of-range values")

    # Round percentage columns
    for col in ["te_pct", "class_pct", "family_pct"]:
        if col in df.columns:
            df[col] = df[col].round(2)

    # Save
    df.to_csv(output_file, index=False)
    logger.info(f"Saved processed data ({len(df)} rows) to {output_file}")


def _load_reference(filepath: str) -> Optional[pd.DataFrame]:
    """Load reference table, auto-detecting delimiter and length column."""
    try:
        # Try tab first, then comma
        for sep in ["\t", ","]:
            df = pd.read_csv(filepath, sep=sep, comment="#")
            if df.shape[1] >= 2:
                break

        # Find length column
        len_col = _find_column(df, [
            "seq_length", "length", "len", "size", "chromEnd",
            "eccDNA_len", "eEnd",
        ])
        if len_col is None:
            # If only 2 columns, assume name + length
            if df.shape[1] == 2:
                df.columns = ["name", "seq_length"]
                return df
            logger.warning(f"Cannot find length column in reference: {list(df.columns)}")
            return None

        if len_col != "seq_length":
            df = df.rename(columns={len_col: "seq_length"})
        df["seq_length"] = pd.to_numeric(df["seq_length"], errors="coerce")
        return df
    except Exception as e:
        logger.warning(f"Failed to load reference {filepath}: {e}")
        return None


def _find_column(df: pd.DataFrame, candidates: list[str]) -> Optional[str]:
    """Find first matching column from candidates."""
    for col in candidates:
        if col in df.columns:
            return col
    return None


def _find_merge_key(
    df: pd.DataFrame, ref_df: pd.DataFrame
) -> Optional[tuple[str, str]]:
    """Find a matching key column between two DataFrames."""
    key_candidates = [
        ("eccDNA_name", "eccDNA_name"),
        ("eccDNA_name", "name"),
        ("name", "name"),
        ("eccDNA_name", "eccDNA_id"),
        ("eccDNA_id", "name"),
    ]
    for df_key, ref_key in key_candidates:
        if df_key in df.columns and ref_key in ref_df.columns:
            return (df_key, ref_key)
    # Fall back to first string column overlap
    for dc in df.columns:
        for rc in ref_df.columns:
            if dc == rc and df[dc].dtype == object:
                return (dc, rc)
    return None
