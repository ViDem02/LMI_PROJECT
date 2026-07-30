#!/usr/bin/env python3
"""
Generate a Markdown report of this codebase: a directory tree followed by
the contents of every text file, each under a "# <path>" heading in a
fenced code block.

Usage:
    python3 generate_report.py

Edit the CONFIG section below to control what gets included.
"""

import os
from pathlib import Path

# ----------------------------- CONFIG -----------------------------------

ROOT_DIR = Path(__file__).resolve().parent
OUTPUT_FILE = ROOT_DIR / "CODEBASE_REPORT.md"

# Folders to exclude entirely (from both the tree and the content dump).
# Paths are relative to ROOT_DIR, use "/" as separator. Matches the folder
# itself and everything under it.
EXCLUDE_DIRS = {
    ".git",
    ".vscode",
    "results",
    "class_code",
    "class_code/slprj",
}

# Individual files to exclude entirely (relative to ROOT_DIR).
EXCLUDE_FILES = {
    ".DS_Store",
    "CODEBASE_REPORT.md",
    "generate_report.py",
}

# File extensions to list in the tree but SKIP in the content dump
# (binary / generated data that doesn't render usefully as text).
BINARY_EXTENSIONS = {
    ".mat", ".png", ".jpg", ".jpeg", ".fig", ".slx", ".slxc",
    ".pdf", ".DS_Store", ".mex", ".mexmaci64", ".mexa64", ".mexw64",
}

# Map file extensions to Markdown code-fence language tags.
LANGUAGE_MAP = {
    ".m": "matlab",
    ".py": "python",
    ".json": "json",
    ".md": "markdown",
    ".txt": "text",
    ".xml": "xml",
    ".yml": "yaml",
    ".yaml": "yaml",
    ".sh": "bash",
    ".gitattributes": "text",
}

# --------------------------------------------------------------------------


def is_excluded_dir(rel_path: Path) -> bool:
    rel_str = rel_path.as_posix()
    return any(rel_str == d or rel_str.startswith(d + "/") for d in EXCLUDE_DIRS)


def is_excluded_file(rel_path: Path) -> bool:
    rel_str = rel_path.as_posix()
    if rel_str in EXCLUDE_FILES or rel_path.name in EXCLUDE_FILES:
        return True
    return any(rel_str.startswith(d + "/") for d in EXCLUDE_DIRS) or rel_str in EXCLUDE_DIRS


def collect_entries():
    """Return sorted list of (Path relative to ROOT_DIR) for dirs/files to show in the tree."""
    entries = []
    for dirpath, dirnames, filenames in os.walk(ROOT_DIR):
        rel_dir = Path(dirpath).relative_to(ROOT_DIR)

        if rel_dir != Path(".") and is_excluded_dir(rel_dir):
            dirnames[:] = []
            continue

        dirnames[:] = sorted(
            d for d in dirnames
            if not is_excluded_dir((rel_dir / d) if rel_dir != Path(".") else Path(d))
        )

        for f in sorted(filenames):
            rel_file = (rel_dir / f) if rel_dir != Path(".") else Path(f)
            if is_excluded_file(rel_file):
                continue
            entries.append(rel_file)

    return sorted(entries, key=lambda p: p.as_posix())


def build_tree_lines():
    """Build an ASCII tree string of ROOT_DIR honoring the exclude rules."""
    lines = [ROOT_DIR.name + "/"]

    def walk(dir_path: Path, rel_path: Path, prefix: str):
        try:
            children = sorted(os.listdir(dir_path), key=str.lower)
        except PermissionError:
            return

        visible = []
        for name in children:
            child_rel = (rel_path / name) if rel_path != Path(".") else Path(name)
            full = dir_path / name
            if full.is_dir():
                if is_excluded_dir(child_rel):
                    continue
                visible.append((name, True, child_rel))
            else:
                if is_excluded_file(child_rel):
                    continue
                visible.append((name, False, child_rel))

        for i, (name, is_dir, child_rel) in enumerate(visible):
            connector = "└── " if i == len(visible) - 1 else "├── "
            lines.append(f"{prefix}{connector}{name}{'/' if is_dir else ''}")
            if is_dir:
                extension = "    " if i == len(visible) - 1 else "│   "
                walk(dir_path / name, child_rel, prefix + extension)

    walk(ROOT_DIR, Path("."), "")
    return lines


def read_text_safely(path: Path):
    try:
        return path.read_text(encoding="utf-8")
    except UnicodeDecodeError:
        return path.read_text(encoding="latin-1")


def main():
    tree_lines = build_tree_lines()
    files = collect_entries()

    out = []
    out.append(f"# Codebase Report: {ROOT_DIR.name}\n")
    out.append("## Directory Tree\n")
    out.append("```")
    out.extend(tree_lines)
    out.append("```\n")

    out.append("## File Contents\n")
    for rel_file in files:
        full_path = ROOT_DIR / rel_file
        ext = full_path.suffix.lower()
        out.append(f"# {rel_file.as_posix()}\n")

        if ext in BINARY_EXTENSIONS:
            out.append("*(binary file, content omitted)*\n")
            continue

        lang = LANGUAGE_MAP.get(ext, "")
        try:
            content = read_text_safely(full_path)
        except Exception as e:
            out.append(f"*(could not read file: {e})*\n")
            continue

        out.append(f"```{lang}")
        out.append(content.rstrip("\n"))
        out.append("```\n")

    OUTPUT_FILE.write_text("\n".join(out), encoding="utf-8")
    print(f"Wrote {OUTPUT_FILE} ({len(files)} files listed)")


if __name__ == "__main__":
    main()
