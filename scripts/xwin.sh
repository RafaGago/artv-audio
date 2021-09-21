#!/bin/bash

xwin_version="xwin-0.1.1"
xwin_prefix="xwin-xwin-0.1.1-x86_64-unknown-linux-musl"

mkdir -p "${xwin_version}"
curl --fail -L "https://github.com/Jake-Shadle/xwin/releases/download/${xwin_version}/${xwin_prefix}.tar.gz" | tar -xzv -C "${xwin_version}" --strip-components=1 "${xwin_prefix}/xwin"

xwin_files="${PWD}/${xwin_version}/files"
${xwin_version}/xwin --accept-license 1 --version 16 splat --output ${xwin_files}

# Fix some symlinks that xwin didn't handle.
ln -s "${xwin_files}/crt/lib/x86_64/msvcrt.lib" "${xwin_files}/crt/lib/x86_64/msvcrtd.lib"
ln -s "${xwin_files}/sdk/include/um/ole2.h" "${xwin_files}/sdk/include/um/Ole2.h"
ln -s "${xwin_files}/sdk/include/um/dbghelp.h" "${xwin_files}/sdk/include/um/Dbghelp.h"
ln -s "${xwin_files}/sdk/lib/um/x86_64/dbghelp.lib" "${xwin_files}/sdk/lib/um/x86_64/DbgHelp.lib"
