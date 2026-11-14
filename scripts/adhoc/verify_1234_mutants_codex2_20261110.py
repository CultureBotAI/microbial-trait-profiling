"""Verify #1234 policy and palette mutants on scratch files; print verdicts, 2026-11-10."""

import importlib.util
import sys
from pathlib import Path
from tempfile import TemporaryDirectory

ROOT = Path(__file__).resolve().parents[2]
TEST = Path("tests/test_legacy_facultative_labels_1234.py")


def main() -> None:
    spec = importlib.util.spec_from_file_location("checks_1234", ROOT / TEST)
    if spec is None or spec.loader is None:
        raise RuntimeError("Cannot load label checks")
    checks = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(checks)
    paths = [Path("scripts") / f"{name}.py" for name in checks.SCRIPTS]
    originals = {path: (ROOT / path).read_bytes() for path in [TEST, *paths]}
    mutants = [
        (path, "return FACULTATIVE_UNION_LABEL", "return FACULTATIVE_RAW_LABEL") for path in paths
    ]
    mutants += [
        (paths[1], 'FACULTATIVE_UNION_LABEL: "#0072b2"', 'FACULTATIVE_UNION_LABEL: "#000000"'),
        (paths[3], "return result", "return frame"),
        (paths[0], "label=resolve_facultative_display_label(),", "label=FACULTATIVE_RAW_LABEL,"),
    ]
    sys.path.insert(0, str(ROOT / "scripts"))
    try:
        with TemporaryDirectory() as directory:
            for index in range(len(mutants) + 1):
                root = Path(directory) / str(index)
                for path, content in originals.items():
                    target = root / path
                    target.parent.mkdir(parents=True, exist_ok=True)
                    target.write_bytes(content)
                if index:
                    path, before, after = mutants[index - 1]
                    target = root / path
                    source = target.read_text()
                    if before not in source:
                        raise RuntimeError(f"Missing anchor {index}")
                    target.write_text(source.replace(before, after, 1))
                try:
                    for name in checks.SCRIPTS:
                        module = checks.load_script(name, root)
                        case = root / name
                        case.mkdir()
                        checks.check_module(module, case)
                except AssertionError:
                    if not index:
                        raise
                    print(f"mutant {index}: killed by AssertionError")
                else:
                    if index:
                        raise RuntimeError(f"mutant {index} survived")
    finally:
        sys.path.pop(0)
        for path, content in originals.items():
            if (ROOT / path).read_bytes() != content:
                raise RuntimeError(f"Changed worktree: {path}")
    print("8/8 killed; worktree targets byte-identical")


if __name__ == "__main__":
    main()
