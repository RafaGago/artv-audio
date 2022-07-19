#!/bin/bash

xwin_version="xwin-0.2.5"
xwin_arch="x86_64"
xwin_prefix="unknown-linux-musl"
xwin_fullname="${xwin_version}-${xwin_arch}-${xwin_prefix}"

wget "https://github.com/Jake-Shadle/xwin/releases/download/0.2.5/${xwin_fullname}.tar.gz"

mkdir -p ${xwin_version}
tar -xzvf ${xwin_fullname}.tar.gz -C ${xwin_version} --strip-components=1 "${xwin_fullname}/xwin"
rm ${xwin_fullname}.tar.gz

xwin_files="${PWD}/${xwin_version}/files"
${xwin_version}/xwin --accept-license --manifest-version 16 splat --include-debug-libs --output ${xwin_files}

## too slow...
#for file in $(find ${xwin_files} -type f); do
#    DIR="$(dirname file)"
#    FILE=$(basename file)
#    # capitalize first letter
#    FILE_CAPITAL="$(basename file | tr [:upper:] [:lower:] | sed -e 's|^.|\U&|')"
#    ln -s "$DIR/$FILE" "$DIR/$FILE_CAPITAL" 2>/dev/null || true
#done

# Manually fix some symlinks that xwin didn't handle.
ln -s "${xwin_files}/sdk/include/um/ole2.h" "${xwin_files}/sdk/include/um/Ole2.h"
ln -s "${xwin_files}/sdk/include/um/dbghelp.h" "${xwin_files}/sdk/include/um/Dbghelp.h"
ln -s "${xwin_files}/sdk/include/shared/dxgi.h" "${xwin_files}/sdk/include/shared/Dxgi.h"
ln -s "${xwin_files}/sdk/lib/um/x86_64/dbghelp.lib" "${xwin_files}/sdk/lib/um/x86_64/DbgHelp.lib"
