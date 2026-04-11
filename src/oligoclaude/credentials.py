"""Secure storage for API keys.

Reads the AlphaGenome API key from (in priority order):
    1. $ALPHAGENOME_API_KEY environment variable
    2. ~/.oligoclaude/credentials.json (file mode 0600)
    3. `dna_api_key` field in the config JSON (legacy, with warning)

Never commit the key to the config file — the legacy path is only a
fallback for backward compatibility and emits a DeprecationWarning.
"""
from __future__ import annotations

import json
import os
import stat
import warnings
from pathlib import Path
from typing import Optional

CRED_DIR = Path.home() / ".oligoclaude"
CRED_PATH = CRED_DIR / "credentials.json"
ENV_VAR = "ALPHAGENOME_API_KEY"

_PLACEHOLDERS = {"", "REPLACE_WITH_YOUR_KEY", "YOUR_ALPHAGENOME_API_KEY"}


def _read_credentials_file() -> dict:
    if not CRED_PATH.exists():
        return {}
    try:
        return json.loads(CRED_PATH.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return {}


def get_alphagenome_api_key(cfg_value: Optional[str] = None) -> Optional[str]:
    """Resolve the AlphaGenome API key from env, credentials file, or config.

    Returns None if no key is available. A legacy config-embedded key
    triggers a DeprecationWarning but is still honored.
    """
    env_key = os.environ.get(ENV_VAR, "").strip()
    if env_key and env_key not in _PLACEHOLDERS:
        return env_key

    file_key = _read_credentials_file().get("alphagenome_api_key", "")
    if isinstance(file_key, str):
        file_key = file_key.strip()
        if file_key and file_key not in _PLACEHOLDERS:
            return file_key

    if cfg_value:
        cfg_value = cfg_value.strip()
        if cfg_value and cfg_value not in _PLACEHOLDERS:
            warnings.warn(
                "dna_api_key is embedded in the config file. Move it to the "
                f"{ENV_VAR} environment variable or run `oligoclaude set-api-key` "
                "(saves to ~/.oligoclaude/credentials.json, mode 0600). "
                "Embedding keys in shared config files is a security risk.",
                DeprecationWarning,
                stacklevel=2,
            )
            return cfg_value

    return None


def save_alphagenome_api_key(key: str) -> Path:
    """Persist the key to ~/.oligoclaude/credentials.json with mode 0600."""
    key = key.strip()
    if not key or key in _PLACEHOLDERS:
        raise ValueError("Refusing to save an empty or placeholder key.")

    CRED_DIR.mkdir(parents=True, exist_ok=True)
    data = _read_credentials_file()
    data["alphagenome_api_key"] = key
    CRED_PATH.write_text(json.dumps(data, indent=2), encoding="utf-8")
    try:
        os.chmod(CRED_PATH, stat.S_IRUSR | stat.S_IWUSR)
    except OSError:
        pass
    return CRED_PATH


def clear_alphagenome_api_key() -> bool:
    """Remove the stored key. Returns True if a key was present."""
    data = _read_credentials_file()
    if "alphagenome_api_key" not in data:
        return False
    del data["alphagenome_api_key"]
    if data:
        CRED_PATH.write_text(json.dumps(data, indent=2), encoding="utf-8")
    else:
        try:
            CRED_PATH.unlink()
        except FileNotFoundError:
            pass
    return True


def require_alphagenome_api_key(cfg_value: Optional[str] = None) -> str:
    """Same as get_alphagenome_api_key, but raise if missing."""
    key = get_alphagenome_api_key(cfg_value)
    if not key:
        raise RuntimeError(
            "No AlphaGenome API key found. Set one of:\n"
            f"  - {ENV_VAR} environment variable\n"
            f"  - Run: oligoclaude set-api-key <KEY>\n"
            "  - (Legacy, discouraged) `dna_api_key` field in the config JSON\n"
            "Or pass --skip-alphagenome to run without it."
        )
    return key
