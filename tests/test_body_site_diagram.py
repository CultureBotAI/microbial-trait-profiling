import pytest

from mtp.body_site_diagram import OTHER_ANCHOR, UNKNOWN_ANCHOR, anatomical_site_base

# Top body-site groups of study_trait_profile_study_metadata.tsv with per-experiment
# labels (#555), plus Nasopharynx, which they replaced.
TOP_GROUPS = [
    "Feces",
    "Mouth",
    "Nasopharynx",
    "Oral cavity",
    "Saliva",
    "Subgingival dental plaque",
    "Throat",
    "Uterine cervix",
    "Vagina",
]


def test_throat_and_pharynx_sit_on_the_neck() -> None:
    neck = anatomical_site_base("Throat")
    assert anatomical_site_base("Oropharynx") == neck
    assert anatomical_site_base("Nasopharynx") != neck
    assert neck != OTHER_ANCHOR


def test_cervix_and_uterine_sites_sit_in_the_pelvis() -> None:
    pelvis = anatomical_site_base("Uterine cervix")
    assert anatomical_site_base("Uterus") == pelvis
    assert pelvis != anatomical_site_base("Vagina")
    assert pelvis != OTHER_ANCHOR


@pytest.mark.parametrize("group", ["Throat", "Uterine cervix"])
def test_new_top_groups_do_not_share_another_groups_anchor(group: str) -> None:
    anchor = anatomical_site_base(group)
    others = [other for other in TOP_GROUPS if other != group]
    assert [other for other in others if anatomical_site_base(other) == anchor] == []


def test_other_and_unknown_stay_off_the_body() -> None:
    assert anatomical_site_base("Other") == OTHER_ANCHOR
    for value in ["unknown", None, " NA ", "nan"]:
        assert anatomical_site_base(value) == UNKNOWN_ANCHOR


def test_unmapped_body_site_group_raises() -> None:
    with pytest.raises(ValueError, match="'Skin'"):
        anatomical_site_base("Skin")
