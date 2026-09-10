#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# publish_wiki.sh -- push docs/ to the repository's GitHub wiki.
#
# GitHub wikis are a separate git repository at <repo>.wiki.git. This copies
# the Markdown in docs/ into a clone of it and pushes.
#
#   bin/publish_wiki.sh git@github.com:USER/PE_align.wiki.git
#   bin/publish_wiki.sh https://github.com/USER/PE_align.wiki.git
#
# The wiki repository must exist first: open the repository's Wiki tab on
# github.com and create any page once. GitHub does not create the wiki repo
# until a first page exists, and cloning before then fails.
# ---------------------------------------------------------------------------
set -euo pipefail

WIKI_URL=${1:?usage: publish_wiki.sh <wiki-repo-url>}
PROJECT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT

echo "cloning $WIKI_URL"
if ! git clone --quiet "$WIKI_URL" "$WORK/wiki"; then
    echo
    echo "ERROR: could not clone the wiki repository." >&2
    echo "Create the first wiki page on github.com (repository -> Wiki ->" >&2
    echo "Create the first page), then re-run this." >&2
    exit 1
fi

cp "$PROJECT_DIR"/docs/*.md "$WORK/wiki/"
cd "$WORK/wiki"

if git diff --quiet && git diff --cached --quiet && [ -z "$(git status --porcelain)" ]; then
    echo "wiki already up to date"
    exit 0
fi

git add -A
git commit -q -m "docs: sync wiki from docs/ ($(date -u +%Y-%m-%dT%H:%M:%SZ))"
git push --quiet
echo "wiki updated: ${WIKI_URL%.wiki.git}/wiki"
