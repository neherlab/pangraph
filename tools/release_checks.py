#!/usr/bin/env python3

import argparse
import re


def parse_args():
    """Create argument parser and return parsed arguments"""
    parser = argparse.ArgumentParser(
        description="""
    Script to perform preliminary checks on release versioning:
    - check that release tag has the correct fom X.Y.Z
    - check that changelog contains the corresponding entry.
    - check that the Project.toml file contains the right version.
    """
    )
    parser.add_argument("--release", help="desired release version: X.X.X", type=str)
    parser.add_argument("--changelog", help="changelog file", type=str)
    parser.add_argument("--project", help="pangraph's Project.toml file", type=str)
    return parser.parse_args()


def capture_first_changelog_entries(changelog_file):
    """Captures and returns the first two lines of the changelog that start with '## '"""

    with open(changelog_file, "r") as f:
        changelog = f.readlines()

    # capture first two entries
    versions = [l.strip() for l in changelog if l.startswith("## ")]
    return versions[0], versions[1]


def check_newer_version(new_release, old_release):
    """Check the correct ordering of versions"""
    x, y, z = map(int, new_release.split("."))
    xo, yo, zo = map(int, old_release.split("."))

    order = x > xo
    order |= (x == xo) and (y > yo)
    order |= (x == xo) and (y == yo) and (z > zo)

    assert (
        order
    ), f"new release version v{new_release} is not newer than previous release v{old_release}"


def capture_toml_version(project_file):
    """Capture line `version = "X.Y.Z"` in toml file and returns the "X.Y.Z" string"""

    with open(project_file, "r") as f:
        project = f.readlines()

    versions = [l for l in project if l.startswith("version =")]
    assert len(versions) == 1, "Missing or extra version line in Project.toml"
    version = versions[0]

    pattern = r'^version = "([\d]+\.[\d]+.[\d]+)"$'
    m = re.match(pattern, version)
    assert bool(
        m
    ), 'version in Project.toml does have the expected pattern: `version = "X.Y.Z"`'
    return m.group(1)


if __name__ == "__main__":

    args = parse_args()
    release = args.release

    # 1) check that the version has the correct pattern
    pattern = r"^[\d]+\.[\d]+\.[\d]+$"
    match = re.match(pattern, release)
    assert (
        match
    ), f"specified release version {release} does not match expected pattern X.Y.Z"

    # 2) check that changelog has the correct corresponding entry
    # and that it is newer than the old one
    next_version, prev_version = capture_first_changelog_entries(args.changelog)
    # checks that the changelog has the release `## vX.Y.Z` entry
    assert (
        next_version == f"## v{release}"
    ), f"""
    First entry in changelog: '{next_version}'
    does not correspond to expected: '## v{release}'
    """
    # checks that this entry is newer than the previous one.
    check_newer_version(release, prev_version[len("## v") :])

    # 3) check that the Project.toml file has the correct version
    version = capture_toml_version(args.project)
    assert (
        version == release
    ), f"""
    version in Project.toml does not correspond to the expected `version = "{release}"`.
    Please modify it.
    """
