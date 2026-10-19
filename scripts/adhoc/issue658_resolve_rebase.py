"""Resolve additive PR #660 rebase conflicts against main.

Inputs: git versions at 83f23aa and 99160fe. Outputs: combined schema, Justfile,
and owner wrappers. Serves #658 / #660; 2026-10-19.
"""

from pathlib import Path

from git import Repo

ROOT = Path(__file__).resolve().parents[2]


def source(revision: str, path: str) -> str:
    return Repo(ROOT).git.show(f"{revision}:{path}") + "\n"


def main() -> None:
    base = "83f23aa3440ee78f291df3227c02b43f470f3aad"
    feature = "99160fee6338ddb5e8f30ed1da06739fefff0599"
    schema = source(base, "schema/workflow.yaml")
    additions = source(feature, "schema/workflow.yaml").split("  RecordedPath:", 1)[1]
    (ROOT / "schema/workflow.yaml").write_text(schema.rstrip() + "\n  RecordedPath:" + additions)
    justfile = source(base, "Justfile")
    reruns = (
        source(feature, "Justfile")
        .split("# Build one shared raw-to-core execution", 1)[1]
        .split("\nsync:\n", 1)[0]
    )
    (ROOT / "Justfile").write_text(
        justfile.replace(
            "\nsync:\n",
            "\n# Build one shared raw-to-core execution" + reruns + "\nsync:\n",
            1,
        )
    )
    for name in ("composition_arm_jed", "directional_consistency_arm_noel"):
        path = f"scripts/rerun_{name}.sh"
        wrapper = source(feature, path).replace(
            "set -euo pipefail\n",
            "set -euo pipefail\n"
            "# Keep rendering headless independently of the caller's session (#719).\n"
            "export MPLBACKEND=Agg\n",
            1,
        )
        (ROOT / path).write_text(wrapper)


if __name__ == "__main__":
    main()
