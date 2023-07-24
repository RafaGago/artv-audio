#!/bin/bash

MIXMAXTRIX_VERSION="v$(cat changelog-mixmaxtrix | head -n 1 | xargs)"
TURBOPACO_VERSION="v$(cat changelog-turbopaco | head -n 1 | xargs)"

ARTIFACT_FOLDER="$1"
GIT_REV="$2"

if [[ ! -d "$ARTIFACT_FOLDER" ]]; then
    echo "Directory argument not found: '$ARTIFACT_FOLDER' " 2>&1
    exit
fi

if ! git log "$GIT_REV" -n 1 > /dev/null 2>&1; then
    echo "Git revision argument doesn't seem to be valid: '$GIT_REV' " 2>&1
    exit
fi

for SRC in $(find "$ARTIFACT_FOLDER" -type f -name "turbopaco-*${GIT_REV}*"); do
    DST="$(echo "$SRC" | sed "s|$GIT_REV|$TURBOPACO_VERSION|g")"
    mv "$SRC" "$DST"
done

for SRC in $(find "$ARTIFACT_FOLDER" -type f -name "mix-maxtrix-*${GIT_REV}*"); do
    DST="$(echo "$SRC" | sed "s|$GIT_REV|$MIXMAXTRIX_VERSION|g")"
    mv "$SRC" "$DST"
done
