"""展示资料摘要，或执行合成对照与空间 placebo 诊断。"""
import argparse
import json
from pathlib import Path

from .domain import summarize
from .io import load
from .synthetic import analyze


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("command", choices=["inspect", "analyze"])
    parser.add_argument("file")
    parser.add_argument("--output", help="把结果写入该路径；不会覆盖输入文件")
    args = parser.parse_args()

    value, digest = load(args.file)
    if args.command == "inspect":
        result = {"source_sha256": digest, "summary": summarize(value)}
    else:
        result = analyze(value, source_sha256=digest)

    text = json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False)
    if args.output:
        target = Path(args.output)
        if target.resolve() == Path(args.file).resolve():
            raise ValueError("输出路径不能与输入资料相同")
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(text + "\n", encoding="utf8")
        print(json.dumps({"written": str(target), "source_sha256": digest}, ensure_ascii=False))
    else:
        print(text)


if __name__ == "__main__":
    main()
