"""Ensure trait-profile investigation figures have non-redundant guild labels."""

from __future__ import annotations

from pathlib import Path

from PIL import Image, ImageDraw, ImageFont, PngImagePlugin

REPO_ROOT = Path(__file__).resolve().parents[1]
FIGURE_ROOT = REPO_ROOT / "trait-visualization" / "trait-profile-investigation"
FOOTER_LINE_COLOR = (203, 213, 225)
TOP_LABEL_FOLDERS = {
    "facultative-aerobe-trait-expansion",
    "gut-anaerobe-fermentation-depletion",
}
UNLABELED_FIGURES = {
    Path("profile_enrichment_guild_assignment_pre_direction_grid_summary.png"),
    Path("profile_enrichment_guild_assignment_summary.png"),
    Path("trait_profile_paired_contrast_infographic.png"),
}

FACULTATIVE_BODY_SITE = (
    "Guild-specific: yes - color encodes the facultative-aerobe guild score by body-site group."
)
FACULTATIVE_V2_BODY_SITE = (
    "Guild-specific: yes - color encodes the facultative-aerobe v2 guild score by body-site group."
)
GUT_BODY_SITE = (
    "Guild-specific: yes - color encodes the gut anaerobe/fermentation guild score "
    "by body-site group."
)
GUT_V2_BODY_SITE = (
    "Guild-specific: yes - color encodes the anaerobe v2 depletion guild score by body-site group."
)
NON_GUILD_DIRECTION_COUNT = (
    "Guild-specific: no - fixed target-taxon direction-count view; "
    "not recomputed from the guild score."
)
NON_GUILD_DIRECTION_COUNT_V2 = (
    "Guild-specific: no - fixed target-taxon direction-count view; "
    "not recomputed from the v2 guild score."
)
NON_GUILD_OCCURRENCE_COUNT = (
    "Guild-specific: no - fixed target-taxon occurrence-count view; "
    "not recomputed from the guild score."
)
NON_GUILD_OCCURRENCE_COUNT_V2 = (
    "Guild-specific: no - fixed target-taxon occurrence-count view; "
    "not recomputed from the v2 guild score."
)
NON_GUILD_TAXA_HEATMAP = (
    "Guild-specific: no - taxon direction-count heatmap; guild traits are annotation context."
)
FACULTATIVE_GUILD_SHIFT = (
    "Guild-specific: yes - shows selected guild traits and their increased-minus-decreased shifts."
)
FACULTATIVE_V2_GUILD_SHIFT = (
    "Guild-specific: yes - shows selected v2 guild traits and their "
    "increased-minus-decreased shifts."
)
GUT_GUILD_SHIFT = (
    "Guild-specific: yes - shows selected guild traits and their decreased-minus-increased shifts."
)
GUT_V2_GUILD_SHIFT = (
    "Guild-specific: yes - shows selected v2 guild traits and their "
    "decreased-minus-increased shifts."
)

FIGURE_LABELS = {
    Path("facultative-aerobe-trait-expansion/anatomical_sites_diagram.png"): (
        FACULTATIVE_BODY_SITE
    ),
    Path("facultative-aerobe-trait-expansion/autism_autoimmune_taxa_signed_heatmap.png"): (
        NON_GUILD_TAXA_HEATMAP
    ),
    Path("facultative-aerobe-trait-expansion/candidate_taxa_increased_profile_tendency.png"): (
        NON_GUILD_DIRECTION_COUNT
    ),
    Path(
        "facultative-aerobe-trait-expansion/condition_increased_guild_signal.png"
    ): "Guild-specific: yes - uses the facultative-aerobe guild score by condition.",
    Path("facultative-aerobe-trait-expansion/guild_trait_direction_shift.png"): (
        FACULTATIVE_GUILD_SHIFT
    ),
    Path("facultative-aerobe-trait-expansion/target_taxa_condition_breakdown_bar.png"): (
        NON_GUILD_OCCURRENCE_COUNT
    ),
    Path("facultative-aerobe-v2/anatomical_sites_diagram.png"): FACULTATIVE_V2_BODY_SITE,
    Path("facultative-aerobe-v2/autism_autoimmune_taxa_signed_heatmap.png"): (
        NON_GUILD_TAXA_HEATMAP
    ),
    Path("facultative-aerobe-v2/candidate_taxa_increased_profile_tendency.png"): (
        NON_GUILD_DIRECTION_COUNT_V2
    ),
    Path(
        "facultative-aerobe-v2/condition_increased_guild_signal.png"
    ): "Guild-specific: yes - uses the facultative-aerobe v2 guild score by condition.",
    Path("facultative-aerobe-v2/condition_target_taxa_stacked_bar.png"): (
        NON_GUILD_OCCURRENCE_COUNT_V2
    ),
    Path("facultative-aerobe-v2/guild_trait_direction_shift.png"): (FACULTATIVE_V2_GUILD_SHIFT),
    Path("facultative-aerobe-v2/target_taxa_condition_breakdown_bar.png"): (
        NON_GUILD_OCCURRENCE_COUNT_V2
    ),
    Path("gut-anaerobe-fermentation-depletion/anatomical_sites_diagram.png"): GUT_BODY_SITE,
    Path("gut-anaerobe-fermentation-depletion/candidate_taxa_depletion_bias.png"): (
        NON_GUILD_DIRECTION_COUNT
    ),
    Path(
        "gut-anaerobe-fermentation-depletion/condition_decreased_guild_signal.png"
    ): "Guild-specific: yes - uses the gut anaerobe/fermentation guild score by condition.",
    Path("gut-anaerobe-fermentation-depletion/guild_trait_direction_shift.png"): (GUT_GUILD_SHIFT),
    Path("gut-anaerobe-fermentation-depletion/target_taxa_condition_breakdown_bar.png"): (
        NON_GUILD_OCCURRENCE_COUNT
    ),
    Path("anaerobe-depletion-v2/anatomical_sites_diagram.png"): GUT_V2_BODY_SITE,
    Path("anaerobe-depletion-v2/anatomical_sites_enriched_profile_diagram.png"): (
        "Guild-specific: yes - counts anatomical-site increased/decreased profiles "
        "with at least one enriched anaerobe v2 guild trait."
    ),
    Path("anaerobe-depletion-v2/candidate_taxa_depletion_bias.png"): (NON_GUILD_DIRECTION_COUNT_V2),
    Path(
        "anaerobe-depletion-v2/condition_decreased_guild_signal.png"
    ): "Guild-specific: yes - uses the anaerobe v2 depletion guild score by condition.",
    Path("anaerobe-depletion-v2/condition_enriched_profile_histogram.png"): (
        "Guild-specific: yes - counts condition-level increased/decreased profiles "
        "with at least one enriched anaerobe v2 guild trait."
    ),
    Path("anaerobe-depletion-v2/condition_target_taxa_stacked_bar.png"): (
        NON_GUILD_OCCURRENCE_COUNT_V2
    ),
    Path("anaerobe-depletion-v2/guild_trait_direction_shift.png"): GUT_V2_GUILD_SHIFT,
    Path("anaerobe-depletion-v2/target_taxa_condition_breakdown_bar.png"): (
        NON_GUILD_OCCURRENCE_COUNT_V2
    ),
}


def load_font(label: str, max_width: int) -> ImageFont.ImageFont:
    """Load a font sized to keep label text within the image width."""
    scratch = ImageDraw.Draw(Image.new("RGB", (1, 1)))
    for size in range(28, 13, -1):
        try:
            font = ImageFont.truetype("arial.ttf", size)
        except OSError:
            font = ImageFont.load_default()
        bbox = scratch.textbbox((0, 0), label, font=font)
        if bbox[2] <= max_width:
            return font
    return ImageFont.load_default()


def is_footer_separator_row(row: list[tuple[int, int, int]]) -> bool:
    """Return whether a sampled row matches the old appended footer separator."""
    matching = 0
    for pixel in row:
        if all(abs(pixel[index] - FOOTER_LINE_COLOR[index]) <= 2 for index in range(3)):
            matching += 1
    return matching / len(row) > 0.8


def crop_footer(image: Image.Image) -> tuple[Image.Image, bool]:
    """Remove the old bottom guild-specificity footer, if present."""
    rgb = image.convert("RGB")
    width, height = rgb.size
    step = max(1, width // 500)
    for y in range(max(0, height - 250), height):
        row = [rgb.getpixel((x, y)) for x in range(0, width, step)]
        if is_footer_separator_row(row):
            return rgb.crop((0, 0, width, y)), True
    return rgb, False


def add_top_label(image: Image.Image, label: str) -> Image.Image:
    """Add a top guild-specificity band for legacy figures lacking title labels."""
    width, height = image.size
    text_x = max(18, int(width * 0.018))
    font = load_font(label, width - (text_x * 2))
    scratch = ImageDraw.Draw(Image.new("RGB", (1, 1)))
    text_bbox = scratch.textbbox((0, 0), label, font=font)
    text_height = text_bbox[3] - text_bbox[1]
    band_height = max(56, text_height + 30)

    output = Image.new("RGB", (width, height + band_height), "white")
    output.paste(image, (0, band_height))
    draw = ImageDraw.Draw(output)
    text_y = max(12, (band_height - text_height) // 2)
    draw.text((text_x, text_y), label, fill=(71, 85, 105), font=font)
    draw.rectangle((0, band_height - 1, width, band_height), fill=FOOTER_LINE_COLOR)
    return output


def process_figure(relative_path: Path, label: str) -> str:
    """Remove footer labels and add top labels only for legacy figure folders."""
    path = FIGURE_ROOT / relative_path
    needs_top_label = relative_path.parts[0] in TOP_LABEL_FOLDERS
    with Image.open(path) as source:
        dpi = source.info.get("dpi")
        existing_top_label = source.info.get("guild_specific_top_label")
        cropped, removed_footer = crop_footer(source)

    changed = removed_footer
    output = cropped
    if needs_top_label and existing_top_label != label:
        output = add_top_label(output, label)
        changed = True

    if not changed:
        return "unchanged"

    metadata = PngImagePlugin.PngInfo()
    if needs_top_label:
        metadata.add_text("guild_specific_top_label", label)
    save_kwargs: dict[str, object] = {"pnginfo": metadata}
    if dpi:
        save_kwargs["dpi"] = dpi
    output.save(path, **save_kwargs)
    if removed_footer and needs_top_label:
        return "removed footer and added top label"
    if removed_footer:
        return "removed footer"
    return "added top label"


def main() -> None:
    existing_figures = {path.relative_to(FIGURE_ROOT) for path in FIGURE_ROOT.rglob("*.png")}
    labelable_figures = existing_figures - UNLABELED_FIGURES
    missing_labels = sorted(labelable_figures - set(FIGURE_LABELS))
    stale_labels = sorted(set(FIGURE_LABELS) - labelable_figures)
    if missing_labels or stale_labels:
        messages = []
        if missing_labels:
            messages.append(
                "Missing guild-specificity labels for: "
                + ", ".join(path.as_posix() for path in missing_labels)
            )
        if stale_labels:
            messages.append(
                "Labels refer to missing figures: "
                + ", ".join(path.as_posix() for path in stale_labels)
            )
        raise ValueError("; ".join(messages))

    updated = 0
    for relative_path, label in sorted(FIGURE_LABELS.items()):
        status = process_figure(relative_path, label)
        if status != "unchanged":
            updated += 1
        print(f"{status}: {relative_path.as_posix()}")
    print(f"Updated {updated} figure(s).")


if __name__ == "__main__":
    main()
