#!/usr/bin/env python3
"""
Installation verification script for DAH7PS project.

This script checks that all required tools and dependencies are installed
and accessible in the current environment.
"""

import sys
import subprocess
import importlib
from pathlib import Path
from typing import Tuple, List


def check_command(cmd: str, version_flag: str = "--version") -> Tuple[bool, str]:
    """Check if a command-line tool is available."""
    try:
        result = subprocess.run(
            [cmd, version_flag],
            capture_output=True,
            text=True,
            timeout=5
        )
        # Check if command produced any output (some tools return non-zero for help)
        output = (result.stdout + result.stderr).strip()
        if output or result.returncode == 0:
            # Get first line from either stdout or stderr
            first_line = output.split('\n')[0] if output else "OK"
            return True, first_line
        else:
            return False, "Command failed"
    except FileNotFoundError:
        return False, "Not found"
    except subprocess.TimeoutExpired:
        return False, "Timeout"
    except Exception as e:
        return False, str(e)


def check_python_package(package: str, import_name: str = None) -> Tuple[bool, str]:
    """Check if a Python package is installed."""
    if import_name is None:
        import_name = package

    try:
        module = importlib.import_module(import_name)
        version = getattr(module, '__version__', 'unknown')
        return True, version
    except ImportError:
        return False, "Not installed"
    except Exception as e:
        return False, str(e)


def check_r_package(package: str) -> Tuple[bool, str]:
    """Check if an R package is installed."""
    try:
        result = subprocess.run(
            ["Rscript", "-e", f"library({package}); cat(as.character(packageVersion('{package}')))"],
            capture_output=True,
            text=True,
            timeout=10
        )
        if result.returncode == 0:
            return True, result.stdout.strip()
        else:
            return False, "Not installed"
    except Exception as e:
        return False, str(e)


def main():
    print("=" * 70)
    print("DAH7PS Project - Installation Check")
    print("=" * 70)
    print()

    all_ok = True

    # Core command-line tools
    print("Checking core bioinformatics tools...")
    print("-" * 70)

    tools = [
        ("hmmsearch", "-h"),
        ("hmmscan", "-h"),
        ("cd-hit", "-h"),
        ("mmseqs", "version"),
        ("mafft", "--version"),
        ("muscle", "-version"),
        ("trimal", "--version"),
        ("iqtree", "--version"),
        ("raxml-ng", "--version"),
    ]

    for tool, flag in tools:
        status, info = check_command(tool, flag)
        status_symbol = "✓" if status else "✗"
        print(f"  {status_symbol} {tool:20s} {info}")
        if not status:
            all_ok = False

    print()

    # Optional tools
    print("Checking optional tools...")
    print("-" * 70)

    optional_tools = [
        ("signalp6", "--version"),
        ("targetp", "--version"),
        ("interproscan.sh", "--version"),
    ]

    for tool, flag in optional_tools:
        status, info = check_command(tool, flag)
        status_symbol = "✓" if status else "○"  # Circle for optional
        print(f"  {status_symbol} {tool:20s} {info}")

    print()

    # Python packages
    print("Checking Python packages...")
    print("-" * 70)

    python_packages = [
        ("biopython", "Bio"),
        ("pandas", "pandas"),
        ("numpy", "numpy"),
        ("scipy", "scipy"),
        ("matplotlib", "matplotlib"),
        ("seaborn", "seaborn"),
        ("ete3", "ete3"),
        ("dendropy", "dendropy"),
        ("snakemake", "snakemake"),
    ]

    for package, import_name in python_packages:
        status, info = check_python_package(package, import_name)
        status_symbol = "✓" if status else "✗"
        print(f"  {status_symbol} {package:20s} v{info}")
        if not status:
            all_ok = False

    print()

    # R packages
    print("Checking R packages...")
    print("-" * 70)

    # First check if R is available
    r_status, r_version = check_command("Rscript", "--version")
    if r_status:
        print(f"  ✓ R                  {r_version}")

        # Core R packages (required)
        r_packages_core = ["ape", "phytools", "ggplot2"]
        # Optional R packages
        r_packages_optional = ["ggtree"]

        for package in r_packages_core:
            status, info = check_r_package(package)
            status_symbol = "✓" if status else "✗"
            print(f"  {status_symbol} {package:20s} v{info}")
            if not status:
                all_ok = False

        # Check optional packages (don't fail if missing)
        for package in r_packages_optional:
            status, info = check_r_package(package)
            status_symbol = "✓" if status else "○"
            print(f"  {status_symbol} {package:20s} v{info} (optional)")
            # Don't set all_ok to False for optional packages
    else:
        print(f"  ✗ R                  Not found")
        all_ok = False

    print()

    # File structure check
    print("Checking project structure...")
    print("-" * 70)

    required_dirs = [
        "data/raw",
        "data/processed",
        "seqs",
        "msa",
        "trees",
        "asr",
        "traits",
        "structures",
        "figs",
        "workflow",
        "scripts",
        "docs",
        "logs",
    ]

    project_root = Path(__file__).parent.parent

    for dir_path in required_dirs:
        full_path = project_root / dir_path
        status = full_path.exists() and full_path.is_dir()
        status_symbol = "✓" if status else "✗"
        print(f"  {status_symbol} {dir_path}")
        if not status:
            all_ok = False

    print()

    # Configuration files
    print("Checking configuration files...")
    print("-" * 70)

    config_files = [
        "workflow/config.yaml",
        "env/dah7ps.yaml",
        "docs/plan.md",
    ]

    for file_path in config_files:
        full_path = project_root / file_path
        status = full_path.exists() and full_path.is_file()
        status_symbol = "✓" if status else "✗"
        print(f"  {status_symbol} {file_path}")
        if not status:
            all_ok = False

    print()
    print("=" * 70)

    if all_ok:
        print("✓ All required components are installed and configured!")
        print()
        print("Next steps:")
        print("  1. Download Pfam database (see env/README.md)")
        print("  2. Configure workflow/config.yaml with your paths")
        print("  3. Run: snakemake --cores 8 --use-conda all")
        return 0
    else:
        print("✗ Some required components are missing.")
        print()
        print("Please install missing dependencies:")
        print("  conda env update -f env/dah7ps.yaml")
        print()
        print("For optional tools, see env/README.md for manual installation instructions.")
        return 1


if __name__ == "__main__":
    sys.exit(main())
