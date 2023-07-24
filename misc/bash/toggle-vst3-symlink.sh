#!/bin/bash

function link() {
    local DIR="$1"
    local VST3FOLDERNAME="$2"


    if [[ -z "$DIR" ]]; then
        echo "Folder to search for \"$VST3FOLDERNAME\" needed as the first argument"
    fi

    SRC=$(find "$DIR" -type d -name "$VST3FOLDERNAME" | head -n 1)
    if [[ -z "$SRC" ]]; then
        echo "Directory \"$VST3FOLDERNAME\" wasn't found in \"$DIR\"" 1>&2
        exit 1
    fi

    if echo "$DIR" | grep -q "win"; then
        VST3DIR="$HOME/.wine/drive_c/Program\ Files/Common\ Files/VST3"
    else
        VST3DIR="$HOME/.vst3"
    fi

    set -x
    rm -f $VST3DIR/$VST3FOLDERNAME
    ln -s $PWD/$SRC $VST3DIR
    set +x
}

DIR="$1"
link "$DIR" "MixMaxTrix.vst3"
link "$DIR" "TurboPaco.vst3"
