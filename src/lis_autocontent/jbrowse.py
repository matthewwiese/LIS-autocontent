"""Configuration of LIS's JBrowse 2 instances, read from jbrowse.yml."""

import functools
import os

import yaml

CONFIG_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "jbrowse.yml")


class JBrowseError(Exception):
    """Raised when jbrowse.yml can't be used."""


@functools.lru_cache(maxsize=None)
def load_config(path=CONFIG_PATH):
    """jbrowse.yml's instances and Jekyll link instance; JBrowseError if unusable."""
    try:
        with open(path, encoding="utf-8") as handle:
            config = yaml.safe_load(handle)
    except (OSError, yaml.YAMLError) as err:
        raise JBrowseError(f"cannot read {path}: {err}") from err
    config = config if isinstance(config, dict) else {}
    instances = config.get("instances")
    if not (
        isinstance(instances, dict)
        and instances
        and all(isinstance(k, str) and isinstance(v, str) for k, v in instances.items())
    ):
        raise JBrowseError(f"{path}: instances must map each id to a URL")
    if config.get("jekyll_instance") not in instances:
        raise JBrowseError(f"{path}: jekyll_instance must name one of the instances")
    return config


def jekyll_url():
    """Base URL of the instance the Jekyll site's resource links open."""
    config = load_config()
    return config["instances"][config["jekyll_instance"]]
