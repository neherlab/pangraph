#!/bin/bash

# Dependencies:
# - python3
# - gh (see this page for installation instructions on linux: https://github.com/cli/cli/blob/trunk/docs/install_linux.md )
# 
# Usage example:
# bash release.sh 0.6.1

set -euo pipefail

# takes as argument the desired version
version=$1

# Directory where this script resides
THIS_DIR="$(
  cd "$(dirname "${BASH_SOURCE[0]}")"
  pwd
)"
CHANGELOG="$THIS_DIR/../CHANGELOG.md"
PROJECT="$THIS_DIR/../Project.toml"
REPO="neherlab/pangraph"

# check authentication
gh auth status >/dev/null

# check for uncommitted changes
if [ -n "$(git status --porcelain)" ]; then
  echo "ERROR: there are uncommitted changes in the repo" >/dev/stderr
  echo "please commit before releasing a new version." >/dev/stderr
  exit 1
fi

# check that the user is on master branch
curr_branch="$(cd ${THIS_DIR} && git branch --show-current)"
if [ "$curr_branch" != "master" ] ; then
    echo "ERROR: repo is on branch ${curr_branch}" >/dev/stderr
    echo "releasing is only possible on master branch" >/dev/stderr
    exit 1
fi

# check that release is not already present in the repo
if [ "release not found" != "$(gh release view "${version}" --repo "$REPO" 2>&1 || true)" ] ; then
    echo "ERROR: desired release version ${version} is already present in the repo" >/dev/stderr
    exit 1
fi

# check that the release version is well formatted, present in CHANGELOG.md and Project.toml file.
python3 "$THIS_DIR/release_checks.py" --release "$version" --changelog "$CHANGELOG" --project "$PROJECT"

# extract release notes from CHANGELOG.md and create a release draft
echo "Submitting new release:" >/dev/stderr
python3 "$THIS_DIR/extract_release_notes.py" "$CHANGELOG" |
  gh release create \
    "${version}" \
    --repo "$REPO" \
    --title "${version}" \
    -d \
    --notes-file -

# Looks like the release appears not immediately, and if an upload is following too quickly
# it fails because it cannot find the release. So let's wait for an arbitrary amount of
# time here until the release is live.
while [ "release not found" == "$(gh release view "${version}" --repo "$REPO" 2>&1 || true)" ]; do
  echo "Waiting for release to go online"
  sleep 2
done

# Check the release once again, in case of other errors
gh release view "${version}" --repo "$REPO" >/dev/null

echo "draft release successfully submitted. Please review on github and approve."