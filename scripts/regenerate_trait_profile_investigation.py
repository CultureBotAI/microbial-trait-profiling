"""Regenerate the current trait-profile investigation outputs."""

from __future__ import annotations

import os
import re
import runpy
import shutil
import sys
import tempfile
import zipfile
from pathlib import Path

import defusedxml.ElementTree as ET

REPO_ROOT = Path(__file__).resolve().parents[1]
INVESTIGATION_DIR = REPO_ROOT / "trait-visualization" / "trait-profile-investigation"
PPTX_PATH = INVESTIGATION_DIR / "trait_profile_investigation_figures.pptx"

FACULTATIVE_DIR = INVESTIGATION_DIR / "aerobe-enrichment-v2"
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
    "ppt/media/image12.png": FACULTATIVE_DIR / "condition_target_taxa_stacked_bar.png",
    "ppt/media/image13.png": ANAEROBE_DIR / "condition_target_taxa_stacked_bar.png",
}

STACKED_SLIDE_SPECS = [
    {
        "source_slide": 5,
        "slide": 12,
        "media": "image12.png",
        "source_title": "Top conditions for increased target-taxon occurrences",
        "title": "Stacked focal taxa by increased-profile condition",
        "presentation_rel_id": "RfacultativeStackedTaxaSlide",
    },
    {
        "source_slide": 10,
        "slide": 13,
        "media": "image13.png",
        "source_title": "Top conditions for decreased target-taxon occurrences",
        "title": "Stacked focal taxa by decreased-profile condition",
        "presentation_rel_id": "RanaerobeStackedTaxaSlide",
    },
]
SLIDE_ORDER = [1, 2, 3, 4, 5, 12, 6, 7, 8, 9, 10, 13, 11]
TOTAL_SLIDES = len(SLIDE_ORDER)
P_NS = "http://schemas.openxmlformats.org/presentationml/2006/main"
R_NS = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
REL_NS = "http://schemas.openxmlformats.org/package/2006/relationships"
CT_NS = "http://schemas.openxmlformats.org/package/2006/content-types"
SLIDE_REL_TYPE = "http://schemas.openxmlformats.org/officeDocument/2006/relationships/slide"
IMAGE_REL_TYPE = "http://schemas.openxmlformats.org/officeDocument/2006/relationships/image"
NOTES_REL_TYPE = "http://schemas.openxmlformats.org/officeDocument/2006/relationships/notesSlide"
SLIDE_CONTENT_TYPE = "application/vnd.openxmlformats-officedocument.presentationml.slide+xml"

PPTX_TEXT_REPLACEMENTS = {
    "ppt/slides/slide1.xml": {
        "FACULTATIVE-AEROBE V2 INCREASED-PROFILE SIGNAL": (
            "AEROBE V2 TRAIT ENRICHMENT IN INCREASED PROFILES"
        ),
        "Target taxa most often reported in increased profiles": (
            "Target taxa most often reported in increased profiles"
        ),
    },
    "ppt/slides/slide2.xml": {
        "FACULTATIVE-AEROBE V2 INCREASED-PROFILE SIGNAL": (
            "AEROBE V2 TRAIT ENRICHMENT IN INCREASED PROFILES"
        ),
        "Conditions with high increased-profile v2 guild signal": (
            "Conditions with aerobe v2 trait enrichment in increased profiles"
        ),
        "Conditions with facultative-aerobe v2 trait enrichment in increased profiles": (
            "Conditions with aerobe v2 trait enrichment in increased profiles"
        ),
    },
    "ppt/slides/slide3.xml": {
        "FACULTATIVE-AEROBE V2 INCREASED-PROFILE SIGNAL": (
            "AEROBE V2 TRAIT ENRICHMENT IN INCREASED PROFILES"
        ),
        "V2 guild traits are higher in increased profiles": (
            "V2 traits enriched in increased profiles"
        ),
    },
    "ppt/slides/slide4.xml": {
        "FACULTATIVE-AEROBE V2 INCREASED-PROFILE SIGNAL": (
            "AEROBE V2 TRAIT ENRICHMENT IN INCREASED PROFILES"
        ),
        "Anatomical sites for increased-profile facultative-aerobe v2 signal": (
            "Anatomical sites for aerobe v2 trait enrichment in increased profiles"
        ),
        "Anatomical sites for facultative-aerobe v2 trait enrichment in increased profiles": (
            "Anatomical sites for aerobe v2 trait enrichment in increased profiles"
        ),
        "Anatomical sites for increased-profile aerobe v2 signal": (
            "Anatomical sites for aerobe v2 trait enrichment in increased profiles"
        ),
    },
    "ppt/slides/slide5.xml": {
        "FACULTATIVE-AEROBE V2 INCREASED-PROFILE SIGNAL": (
            "AEROBE V2 TRAIT ENRICHMENT IN INCREASED PROFILES"
        ),
        "Top conditions for increased target-taxon occurrences": (
            "Top conditions for increased-profile target-taxon occurrences"
        ),
    },
    "ppt/slides/slide12.xml": {
        "FACULTATIVE-AEROBE V2 INCREASED-PROFILE SIGNAL": (
            "AEROBE V2 TRAIT ENRICHMENT IN INCREASED PROFILES"
        ),
    },
    "ppt/slides/slide6.xml": {
        "GUT ANAEROBE V2 DEPLETION": "ANAEROBE V2 DEPLETION SIGNAL",
        "Target taxa most often reported in decreased profiles": (
            "Target taxa most often reported in decreased profiles"
        ),
    },
    "ppt/slides/slide7.xml": {
        "GUT ANAEROBE V2 DEPLETION": "ANAEROBE V2 DEPLETION SIGNAL",
        "Conditions with high decreased-profile gut anaerobe v2 guild signal": (
            "Conditions with anaerobe v2 trait enrichment in decreased profiles"
        ),
        "Conditions with high decreased-profile anaerobe v2 guild signal": (
            "Conditions with anaerobe v2 trait enrichment in decreased profiles"
        ),
    },
    "ppt/slides/slide8.xml": {
        "GUT ANAEROBE V2 DEPLETION": "ANAEROBE V2 DEPLETION SIGNAL",
        "V2 guild traits are higher in decreased gut-context profiles": (
            "V2 traits enriched in decreased profiles"
        ),
        "V2 guild traits are higher in decreased profiles": (
            "V2 traits enriched in decreased profiles"
        ),
    },
    "ppt/slides/slide9.xml": {
        "GUT ANAEROBE V2 DEPLETION": "ANAEROBE V2 DEPLETION SIGNAL",
        "Anatomical sites for decreased-profile gut anaerobe v2 signal": (
            "Anatomical sites for anaerobe v2 trait enrichment in decreased profiles"
        ),
        "Anatomical sites for decreased-profile anaerobe v2 signal": (
            "Anatomical sites for anaerobe v2 trait enrichment in decreased profiles"
        ),
    },
    "ppt/slides/slide10.xml": {
        "GUT ANAEROBE V2 DEPLETION": "ANAEROBE V2 DEPLETION SIGNAL",
        "Top conditions for decreased target-taxon occurrences": (
            "Top conditions for decreased-profile target-taxon occurrences"
        ),
    },
}


def run_script(script_name: str) -> None:
    """Run one repo script using this interpreter."""
    script = REPO_ROOT / "scripts" / script_name
    print(f"==> {script.relative_to(REPO_ROOT)}")
    previous_cwd = Path.cwd()
    previous_argv = sys.argv[:]
    try:
        os.chdir(REPO_ROOT)
        sys.argv = [str(script)]
        runpy.run_path(str(script), run_name="__main__")
    finally:
        sys.argv = previous_argv
        os.chdir(previous_cwd)


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


def decode_xml(data: bytes) -> str:
    return data.decode("utf-8-sig")


def encode_xml(text: str) -> bytes:
    return text.encode("utf-8")


def replace_page_label(slide_xml: str, page_label: str) -> str:
    """Replace a slide footer page label such as 5/11."""
    return re.sub(r"(<a:t>)\d+/\d+(</a:t>)", rf"\g<1>{page_label}\2", slide_xml, count=1)


def remove_notes_relationships(rels_xml: str) -> str:
    root = ET.fromstring(rels_xml)
    for relationship in list(root):
        if relationship.attrib.get("Type") == NOTES_REL_TYPE:
            root.remove(relationship)
    return ET.tostring(root, encoding="unicode", xml_declaration=True)


def set_first_image_relationship(rels_xml: str, media_name: str) -> str:
    root = ET.fromstring(rels_xml)
    for relationship in root:
        if relationship.attrib.get("Type") == IMAGE_REL_TYPE:
            relationship.set("Target", f"/ppt/media/{media_name}")
            break
    return ET.tostring(root, encoding="unicode", xml_declaration=True)


def ensure_slide_content_type(entries: dict[str, bytes], slide_number: int) -> None:
    root = ET.fromstring(decode_xml(entries["[Content_Types].xml"]))
    part_name = f"/ppt/slides/slide{slide_number}.xml"
    for override in root:
        if override.attrib.get("PartName") == part_name:
            return
    root.append(
        root.makeelement(
            f"{{{CT_NS}}}Override",
            {
                "PartName": part_name,
                "ContentType": SLIDE_CONTENT_TYPE,
            },
        )
    )
    entries["[Content_Types].xml"] = encode_xml(
        ET.tostring(root, encoding="unicode", xml_declaration=True)
    )


def normalize_slide_target(target: str) -> str:
    target = target.lstrip("/")
    if target.startswith("slides/"):
        return f"ppt/{target}"
    return target


def ensure_presentation_slide_relationships(entries: dict[str, bytes]) -> dict[int, str]:
    root = ET.fromstring(decode_xml(entries["ppt/_rels/presentation.xml.rels"]))
    target_to_id = {
        normalize_slide_target(relationship.attrib.get("Target", "")): relationship.attrib["Id"]
        for relationship in root
        if relationship.attrib.get("Type") == SLIDE_REL_TYPE
    }
    slide_rel_ids: dict[int, str] = {}
    for slide_number in SLIDE_ORDER:
        target = f"ppt/slides/slide{slide_number}.xml"
        relationship_id = target_to_id.get(target)
        if relationship_id is None:
            spec = next(item for item in STACKED_SLIDE_SPECS if item["slide"] == slide_number)
            relationship_id = str(spec["presentation_rel_id"])
            root.append(
                root.makeelement(
                    f"{{{REL_NS}}}Relationship",
                    {
                        "Type": SLIDE_REL_TYPE,
                        "Target": f"/{target}",
                        "Id": relationship_id,
                    },
                )
            )
        slide_rel_ids[slide_number] = relationship_id

    entries["ppt/_rels/presentation.xml.rels"] = encode_xml(
        ET.tostring(root, encoding="unicode", xml_declaration=True)
    )
    return slide_rel_ids


def update_presentation_slide_order(
    entries: dict[str, bytes],
    slide_rel_ids: dict[int, str],
) -> None:
    root = ET.fromstring(decode_xml(entries["ppt/presentation.xml"]))
    slide_id_list = root.find(f"{{{P_NS}}}sldIdLst")
    if slide_id_list is None:
        msg = "Missing slide ID list in ppt/presentation.xml"
        raise ValueError(msg)

    existing_ids = {
        slide_id.attrib.get(f"{{{R_NS}}}id"): int(slide_id.attrib["id"])
        for slide_id in slide_id_list
        if slide_id.attrib.get(f"{{{R_NS}}}id") and slide_id.attrib.get("id")
    }
    next_id = max(existing_ids.values(), default=255) + 1
    for child in list(slide_id_list):
        slide_id_list.remove(child)
    for slide_number in SLIDE_ORDER:
        relationship_id = slide_rel_ids[slide_number]
        slide_id_value = existing_ids.get(relationship_id)
        if slide_id_value is None:
            slide_id_value = next_id
            next_id += 1
        slide_id_list.append(
            slide_id_list.makeelement(
                f"{{{P_NS}}}sldId",
                {
                    "id": str(slide_id_value),
                    f"{{{R_NS}}}id": relationship_id,
                },
            )
        )

    entries["ppt/presentation.xml"] = encode_xml(
        ET.tostring(root, encoding="unicode", xml_declaration=True)
    )


def ensure_stacked_slides(entries: dict[str, bytes]) -> None:
    """Add or refresh the two stacked target-taxon slides in the PPTX package."""
    for spec in STACKED_SLIDE_SPECS:
        source_slide = int(spec["source_slide"])
        slide_number = int(spec["slide"])
        source_slide_name = f"ppt/slides/slide{source_slide}.xml"
        source_rels_name = f"ppt/slides/_rels/slide{source_slide}.xml.rels"
        slide_name = f"ppt/slides/slide{slide_number}.xml"
        rels_name = f"ppt/slides/_rels/slide{slide_number}.xml.rels"
        if source_slide_name not in entries or source_rels_name not in entries:
            msg = f"Missing source slide template for slide {slide_number}"
            raise FileNotFoundError(msg)

        slide_xml = decode_xml(entries.get(slide_name, entries[source_slide_name]))
        for source_title in (
            str(spec["source_title"]),
            "Top conditions for increased target-taxon occurrences",
            "Top conditions for decreased target-taxon occurrences",
            "Top conditions for increased-profile target-taxon occurrences",
            "Top conditions for decreased-profile target-taxon occurrences",
        ):
            slide_xml = slide_xml.replace(source_title, str(spec["title"]))
        entries[slide_name] = encode_xml(slide_xml)

        rels_xml = decode_xml(entries.get(rels_name, entries[source_rels_name]))
        rels_xml = remove_notes_relationships(rels_xml)
        rels_xml = set_first_image_relationship(rels_xml, str(spec["media"]))
        entries[rels_name] = encode_xml(rels_xml)
        media_name = f"ppt/media/{spec['media']}"
        entries[media_name] = SLIDE_IMAGE_SOURCES[media_name].read_bytes()
        ensure_slide_content_type(entries, slide_number)

    for position, slide_number in enumerate(SLIDE_ORDER, start=1):
        slide_name = f"ppt/slides/slide{slide_number}.xml"
        if slide_name in entries:
            slide_xml = decode_xml(entries[slide_name])
            entries[slide_name] = encode_xml(
                replace_page_label(slide_xml, f"{position}/{TOTAL_SLIDES}")
            )

    slide_rel_ids = ensure_presentation_slide_relationships(entries)
    update_presentation_slide_order(entries, slide_rel_ids)


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
        with zipfile.ZipFile(PPTX_PATH, "r") as source:
            entries = {
                item.filename: patch_pptx_entry(item.filename, source.read(item.filename))
                for item in source.infolist()
            }
        ensure_stacked_slides(entries)
        with zipfile.ZipFile(
            temp_path,
            "w",
            compression=zipfile.ZIP_DEFLATED,
        ) as target:
            for name, data in entries.items():
                target.writestr(name, data)
        shutil.move(str(temp_path), str(PPTX_PATH))
    finally:
        if temp_path.exists():
            temp_path.unlink()

    print(f"==> refreshed {PPTX_PATH.relative_to(REPO_ROOT)}")


def main() -> None:
    run_script("investigate_facultative_aerobe_expansion.py")
    run_script("investigate_facultative_aerobe_condition_followup.py")
    run_script("investigate_gut_anaerobe_depletion.py")
    run_script("create_condition_body_site_guild_maps.py")
    run_script("create_all_bugsigdb_condition_body_site_maps.py")
    run_script("create_profile_enrichment_guild_assignment_diagram.py")
    run_script("create_trait_profile_infographic.py")
    run_script("analyze_parkinsons_disease_trait_profiles.py")
    run_script("annotate_trait_profile_investigation_figures.py")
    refresh_powerpoint()
    print("Regenerated current trait-profile investigation outputs.")


if __name__ == "__main__":
    main()
