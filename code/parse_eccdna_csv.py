import pandas as pd
import re
import argparse


def parse_seqname(s):
    try:
        parts = s.split("|")
        if len(parts) != 2:
            return pd.Series([None, None, None, None])
        chrom_info = parts[1]
        chrom, start, end = chrom_info.split("-")
        start = int(float(start))
        end = int(float(end))
        length = end - start + 1
        return pd.Series([chrom, start, end, length])
    except Exception:
        return pd.Series([None, None, None, None])


def parse_attribute(attr):
    try:
        motif_match = re.search(r'Target "Motif:([^"]+)"\s+(\d+)\s+(\d+)', str(attr))
        if motif_match:
            motif = motif_match.group(1)
            motif_start = int(float(motif_match.group(2)))
            motif_end = int(float(motif_match.group(3)))
            return pd.Series([motif, motif_start, motif_end])
        else:
            return pd.Series([None, None, None])
    except Exception:
        return pd.Series([None, None, None])


def main(input_file, output_file):
    df = pd.read_csv(input_file)
    df[["chr", "seq_start", "seq_end", "seq_length"]] = df["seqname"].apply(
        parse_seqname
    )
    df[["motif", "motif_start", "motif_end"]] = df["attribute"].apply(parse_attribute)

    # 计算 motif_length
    df["motif_length"] = df.apply(
        lambda row: (
            int(row["motif_end"]) - int(row["motif_start"]) + 1
            if pd.notnull(row["motif_start"]) and pd.notnull(row["motif_end"])
            else None
        ),
        axis=1,
    )

    # motif_length / seq_length * 100, 保留两位小数
    df["motif_percent"] = df.apply(
        lambda row: (
            round(row["motif_length"] / row["seq_length"] * 100, 2)
            if pd.notnull(row["motif_length"])
            and pd.notnull(row["seq_length"])
            and row["seq_length"] != 0
            else None
        ),
        axis=1,
    )

    # 所有整数列处理为 int（忽略 None），如果不为 None
    for col in [
        "seq_start",
        "seq_end",
        "seq_length",
        "motif_start",
        "motif_end",
        "motif_length",
    ]:
        df[col] = df[col].apply(lambda x: int(x) if pd.notnull(x) else None)

    # 保存
    df.to_csv(output_file, index=False)
    print(f"处理完成！已输出文件: {output_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Parse eccDNA csv and expand columns")
    parser.add_argument(
        "--input", "-i", type=str, required=True, help="输入文件名（csv）"
    )
    parser.add_argument(
        "--output", "-o", type=str, required=True, help="输出文件名（csv）"
    )
    args = parser.parse_args()
    main(args.input, args.output)
