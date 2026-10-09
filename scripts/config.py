"""
Backward-compatible imports.

Legacy analysis scripts can keep using:

    from config import *

The canonical workflow package is xtt_project/workflow.
"""

from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[1]
XTT_PROJECT = PROJECT_ROOT  # alias legado; camada xtt_project/ removida
SCRIPTS_METADATA = Path(__file__).resolve().parent / "metadata"

# Insert order matters: sys.path.insert(0, ...) pushes the newest entry to
# the front, so the LAST path inserted here wins for "import workflow" and
# similar top-level package names. XTT_PROJECT must win over PROJECT_ROOT,
# otherwise "from workflow.config.paths import *" silently resolves to the
# stale duplicate at PROJECT_ROOT/workflow instead of the canonical
# xtt_project/workflow.
for path in (SCRIPTS_METADATA, PROJECT_ROOT):
    path_str = str(path)
    if path_str not in sys.path:
        sys.path.insert(0, path_str)

from workflow.config.paths import *
