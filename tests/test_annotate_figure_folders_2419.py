"""Every figure folder under trait-profile-investigation/ is labelled or exempt (#2419).

`scripts/annotate_trait_profile_investigation_figures.py` refuses any PNG under
trait-profile-investigation/ that has neither a guild-specificity label nor an exempt
folder. `just taxa-coverage-548` (#2386) writes coverage histograms into
taxa-coverage-tier1-2-subject-controls/, which was neither, so the last step of
`just trait-profile-investigation` failed once that recipe had run.

Which folders are checked: the Justfile names output folders under
trait-profile-investigation/ on a command line only for that one recipe; every other
script sets its folder in code. The folder test below therefore covers the folders the
Justfile names, and the annotate script's own check covers the rest when the recipe runs.
"""

import importlib.util
import re
import sys
from pathlib import Path
from types import ModuleType

import pytest
from PIL import Image

REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = REPO_ROOT / "scripts" / "annotate_trait_profile_investigation_figures.py"
JUSTFILE_FOLDER = re.compile(
    r"trait-visualization/trait-profile-investigation/([A-Za-z0-9_-]+)/?(?=[\s\"']|$)"
)


def load_annotate() -> ModuleType:
    spec = importlib.util.spec_from_file_location("annotate_2419", SCRIPT)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


annotate = load_annotate()


def covered_folders() -> set[str]:
    labelled = {path.parts[0] for path in annotate.FIGURE_LABELS if len(path.parts) > 1}
    exempt: set[str] = annotate.UNLABELED_FOLDERS | annotate.TOP_LABEL_FOLDERS
    return labelled | exempt


def test_justfile_names_the_taxa_coverage_548_folder() -> None:
    folders = set(JUSTFILE_FOLDER.findall((REPO_ROOT / "Justfile").read_text(encoding="utf-8")))
    assert "taxa-coverage-tier1-2-subject-controls" in folders


def test_every_justfile_output_folder_is_labelled_or_exempt() -> None:
    folders = set(JUSTFILE_FOLDER.findall((REPO_ROOT / "Justfile").read_text(encoding="utf-8")))
    assert folders
    assert sorted(folders - covered_folders()) == []


@pytest.mark.parametrize("folder", ["taxa-coverage", "taxa-coverage-tier1-2-subject-controls"])
def test_annotate_accepts_coverage_histograms(
    folder: str, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    (tmp_path / folder).mkdir()
    Image.new("RGB", (4, 4), "white").save(
        tmp_path / folder / "taxa_coverage_occurrences_vs_unique.png"
    )
    monkeypatch.setattr(annotate, "FIGURE_ROOT", tmp_path)
    monkeypatch.setattr(annotate, "FIGURE_LABELS", {})

    annotate.main()

    assert capsys.readouterr().out.strip() == "Updated 0 figure(s)."


def test_annotate_still_refuses_an_unknown_folder(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    (tmp_path / "new-folder").mkdir()
    Image.new("RGB", (4, 4), "white").save(tmp_path / "new-folder" / "figure.png")
    monkeypatch.setattr(annotate, "FIGURE_ROOT", tmp_path)
    monkeypatch.setattr(annotate, "FIGURE_LABELS", {})

    with pytest.raises(
        ValueError, match=re.escape("Missing guild-specificity labels for: new-folder/figure.png")
    ):
        annotate.main()
