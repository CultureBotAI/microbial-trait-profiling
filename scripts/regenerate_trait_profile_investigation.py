"""Regenerate the current trait-profile investigation outputs."""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
import tempfile
import zipfile
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
INVESTIGATION_DIR = REPO_ROOT / "trait-visualization" / "trait-profile-investigation"
PPTX_PATH = INVESTIGATION_DIR / "trait_profile_investigation_figures.pptx"

FACULTATIVE_DIR = INVESTIGATION_DIR / "facultative-aerobe-v2"
ANAEROBE_DIR = INVESTIGATION_DIR / "anaerobe-depletion-v2"

SLIDE_IMAGE_SOURCES = {
    "ppt/media/image.png": FACULTATIVE_DIR / "candidate_taxa_increased_profile_tendency.png",
    "ppt/media/image2.png": FACULTATIVE_DIR / "condition_increased_guild_signal.png",
    "ppt/media/image3.png": FACULTATIVE_DIR / "guild_trait_direction_shift.png",
    "ppt/media/image4.png": FACULTATIVE_DIR / "anatomical_sites_diagram.png",
    "ppt/media/image5.png": FACULTATIVE_DIR / "target_taxa_condition_breakdown_bar.png",
    "ppt/media/image6.png": ANAEROBE_DIR / "candidate_taxa_depletion_bias.png",
    "ppt/media/image7.png": ANAEROBE_DIR / "condition_decreased_guild_signal.png",
    "ppt/media/image8.png": ANAEROBE_DIR / "guild_trait_direction_shift.png",
    "ppt/media/image9.png": ANAEROBE_DIR / "anatomical_sites_diagram.png",
    "ppt/media/image10.png": ANAEROBE_DIR / "target_taxa_condition_breakdown_bar.png",
    "ppt/media/image11.png": INVESTIGATION_DIR / "profile_enrichment_guild_assignment_summary.png",
}

PPTX_TEXT_REPLACEMENTS = {
    "ppt/slides/slide6.xml": {
        "GUT ANAEROBE V2 DEPLETION": "ANAEROBE V2 DEPLETION SIGNAL",
    },
    "ppt/slides/slide7.xml": {
        "GUT ANAEROBE V2 DEPLETION": "ANAEROBE V2 DEPLETION SIGNAL",
        "Conditions with high decreased-profile gut anaerobe v2 guild signal": (
            "Conditions with high decreased-profile anaerobe v2 guild signal"
        ),
    },
    "ppt/slides/slide8.xml": {
        "GUT ANAEROBE V2 DEPLETION": "ANAEROBE V2 DEPLETION SIGNAL",
        "V2 guild traits are higher in decreased gut-context profiles": (
            "V2 guild traits are higher in decreased profiles"
        ),
    },
    "ppt/slides/slide9.xml": {
        "GUT ANAEROBE V2 DEPLETION": "ANAEROBE V2 DEPLETION SIGNAL",
        "Anatomical sites for decreased-profile gut anaerobe v2 signal": (
            "Anatomical sites for decreased-profile anaerobe v2 signal"
        ),
    },
    "ppt/slides/slide10.xml": {
        "GUT ANAEROBE V2 DEPLETION": "ANAEROBE V2 DEPLETION SIGNAL",
    },
}


def run_script(script_name: str) -> None:
    """Run one repo script using this interpreter."""
    script = REPO_ROOT / "scripts" / script_name
    print(f"==> {script.relative_to(REPO_ROOT)}")
    subprocess.run([sys.executable, str(script)], cwd=REPO_ROOT, check=True)


def patch_pptx_entry(name: str, data: bytes) -> bytes:
    """Return replacement bytes for a PPTX zip entry."""
    if name in SLIDE_IMAGE_SOURCES:
        return SLIDE_IMAGE_SOURCES[name].read_bytes()
    if name in PPTX_TEXT_REPLACEMENTS:
        text = data.decode("utf-8")
        for old, new in PPTX_TEXT_REPLACEMENTS[name].items():
            text = text.replace(old, new)
        return text.encode("utf-8")
    return data


def refresh_powerpoint() -> None:
    """Replace slide media in the existing PowerPoint deck with regenerated PNGs."""
    if not PPTX_PATH.exists():
        raise FileNotFoundError(f"Missing PowerPoint deck: {PPTX_PATH}")

    missing_sources = [path for path in SLIDE_IMAGE_SOURCES.values() if not path.exists()]
    if missing_sources:
        missing = "\n".join(str(path) for path in missing_sources)
        raise FileNotFoundError(f"Missing slide image source(s):\n{missing}")

    fd, temp_name = tempfile.mkstemp(suffix=".pptx", dir=str(PPTX_PATH.parent))
    os.close(fd)
    temp_path = Path(temp_name)

    try:
        with (
            zipfile.ZipFile(PPTX_PATH, "r") as source,
            zipfile.ZipFile(
                temp_path,
                "w",
                compression=zipfile.ZIP_DEFLATED,
            ) as target,
        ):
            for item in source.infolist():
                target.writestr(item, patch_pptx_entry(item.filename, source.read(item.filename)))
        shutil.move(str(temp_path), str(PPTX_PATH))
    finally:
        if temp_path.exists():
            temp_path.unlink()

    print(f"==> refreshed {PPTX_PATH.relative_to(REPO_ROOT)}")


def main() -> None:
    run_script("investigate_facultative_aerobe_expansion.py")
    run_script("investigate_gut_anaerobe_depletion.py")
    run_script("create_profile_enrichment_guild_assignment_diagram.py")
    run_script("annotate_trait_profile_investigation_figures.py")
    refresh_powerpoint()
    print("Regenerated current trait-profile investigation outputs.")


if __name__ == "__main__":
    main()
