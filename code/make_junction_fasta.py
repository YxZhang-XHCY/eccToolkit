#!/usr/bin/env python3
"""
make_junction_fasta.py
从环状DNA序列创建连接点FASTA文件
"""

from Bio import SeqIO
import argparse
import os


def parse_arguments():
    """解析命令行参数"""
    parser = argparse.ArgumentParser(
        description="从环状DNA序列创建连接点FASTA文件",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
使用示例:
    python make_junction_fasta.py -i all_cecc_renamed.fa -o all_cecc_junc.fa
    python make_junction_fasta.py -i input.fa -o output.fa -k 500
        """,
    )

    parser.add_argument("-i", "--input", required=True, help="输入FASTA文件路径")
    parser.add_argument("-o", "--output", required=True, help="输出FASTA文件路径")
    parser.add_argument(
        "-k",
        "--junction-length",
        type=int,
        default=300,
        help="连接点长度，从末端和起始各取的碱基数 (default: 300)",
    )

    return parser.parse_args()


def main():
    """主函数"""
    args = parse_arguments()

    # 检查输入文件是否存在
    if not os.path.exists(args.input):
        print(f"错误: 输入文件不存在: {args.input}")
        return

    # 创建输出目录
    output_dir = os.path.dirname(args.output)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)

    processed_count = 0

    try:
        with open(args.input) as fa_in, open(args.output, "w") as fa_out:
            for rec in SeqIO.parse(fa_in, "fasta"):
                # 检查序列长度是否足够
                if len(rec.seq) < 2 * args.junction_length:
                    print(
                        f"警告: 序列 {rec.id} 长度({len(rec.seq)})小于2倍连接点长度({2*args.junction_length})，跳过"
                    )
                    continue

                # 创建连接点序列 (末端 + 起始)
                junc_seq = (
                    rec.seq[-args.junction_length :] + rec.seq[: args.junction_length]
                )

                # 更新序列信息
                new_rec = rec
                new_rec.id += "_JUNC"
                new_rec.description = (
                    f"junction_of_{rec.id}_length_{args.junction_length}"
                )
                new_rec.seq = junc_seq

                # 写入输出文件
                SeqIO.write(new_rec, fa_out, "fasta")
                processed_count += 1

        print(f"成功处理 {processed_count} 个序列")
        print(f"连接点FASTA文件已保存到: {args.output}")

    except Exception as e:
        print(f"错误: {e}")


if __name__ == "__main__":
    main()
