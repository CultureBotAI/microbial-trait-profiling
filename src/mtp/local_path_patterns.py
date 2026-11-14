"""Shared local-path spellings for provenance and the repository guard.

Callers explicitly choose case sensitivity and whether matches name their
family. Additional POSIX slashless-prefix, account-name and filesystem-root checks belong to
provenance validation and are intentionally separate from these text patterns.
"""

from __future__ import annotations

import re

# Spell alternatives so their source does not match the repository's own scan.
# Newly recognized Windows forms need an account after the separators. A single
# backslash followed by a JSON escape introducer is ambiguous raw text: accept it.
# Actual backslashes in serialized JSON are doubled. Escape letters stay case-sensitive
# even when the consumer matches Users case-insensitively.
_WINDOWS_HOME_TAIL = (
    r'(?:[\\]{2}[\\/]*|[\\](?!(?-i:[\\ntrbf"/u]))|/[\\/]*)'
    r"(?=[\w.-])"
)
USER_HOME_SOURCE = (
    r"/(?:Users|home)/(?:[\w.-]|\\u[0-9A-Fa-f]{4})+|[A-Za-z]:[\\/]+Users[\\/]+"
    rf"|(?:[A-Za-z]:|(?<![\\])[\\]+)Users{_WINDOWS_HOME_TAIL}"
    rf"|(?<![\w])(?<!\w-)(?:Users)(?=[\\]){_WINDOWS_HOME_TAIL}"
)
# Windows also includes drive-relative and rooted/UNC components; bare relative
# backslash forms retain the ordinary-word/internal-hyphen exemption (#1830).
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
