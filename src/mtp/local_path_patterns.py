"""Shared absolute local-path spellings for provenance and the repository guard.

Callers explicitly choose case sensitivity and whether matches name their
family. Slashless-prefix, account-name and filesystem-root checks belong to
provenance validation and are intentionally separate from these text patterns.
"""

from __future__ import annotations

import re

# Spell alternatives so their source does not match the repository's own scan.
USER_HOME_SOURCE = r"/(?:Users|home)/[A-Za-z0-9._-]+|[A-Za-z]:[\\/]+Users[\\/]+"
# Keep the final separator as a lookahead: it may also start a user-home hit.
TEMP_ROOT_SOURCE = r"(?:/private)?/var/(?:folders)(?=/)|/private/(?:tmp)(?=/)"


def compile_absolute_local_paths(*, flags: int, classified: bool) -> re.Pattern[str]:
    """Compile the shared alternatives with an explicit consumer matching policy."""
    source = (
        f"(?P<user_home>{USER_HOME_SOURCE})|(?P<temp_root>{TEMP_ROOT_SOURCE})"
        if classified
        else f"{USER_HOME_SOURCE}|{TEMP_ROOT_SOURCE}"
    )
    return re.compile(source, flags)
