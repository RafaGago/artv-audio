#!/bin/bash

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd)
TTF=$(realpath "$1")
cd "$SCRIPT_DIR"

if [[ ! -f "$TTF" ]]; then
    echo "$TTF is not a file" 2>&1
    exit 1
fi
cp "$TTF" . # xxd generates names including paths

TTF=$(basename "$TTF")
DSTCPP="$TTF".cpp
DSTHPP="$TTF".hpp

printf "namespace artv::font {\n\n" > $DSTCPP
xxd -i "$TTF" >> $DSTCPP
rm "$TTF"

LENLINE=$(cat $DSTCPP | tail -n 1)
head -n -1 "$DSTCPP" > temp
mv temp "$DSTCPP"
echo "} // namespace artv::font" >>  $DSTCPP

VARNAME=$(echo "$LENLINE" | cut -d ' ' -f 3 | sed 's|_len||g')
SIZE=$(echo "$LENLINE" | cut -d '=' -f 2 | sed 's|[^0-9]||g')

cat << EOF > $DSTHPP
#pragma once

namespace artv::font {
extern unsigned char ${VARNAME}[${SIZE}];
} // namespace artv::font
EOF
