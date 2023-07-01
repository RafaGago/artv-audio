#!/bin/bash

xwin_version="xwin-0.2.14"
xwin_arch="x86_64"
xwin_prefix="unknown-linux-musl"
xwin_fullname="${xwin_version}-${xwin_arch}-${xwin_prefix}"

# Mount a case insensitive filesystem
mkdir -p ${xwin_version}
mkdir -p ${xwin_version}.ciopfs.data
umount ${xwin_version} | true
ciopfs ${xwin_version}.ciopfs.data ${xwin_version} -o nonempty | exit 1
# leave a script to remind myself that "ciopfs" has to be mounted before building
printf "#!/bin/bash\nciopfs ${xwin_version}.ciopfs.data ${xwin_version} -o nonempty\n" > ${xwin_version}.ciopfs.mount.sh
chmod +x ${xwin_version}.ciopfs.mount.sh

# Get xwin
if [[ ! -x "${xwin_version}/xwin" ]]; then
    wget "https://github.com/Jake-Shadle/xwin/releases/download/0.2.14/${xwin_fullname}.tar.gz"
    tar -xzvf ${xwin_fullname}.tar.gz -C ${xwin_version} --strip-components=1 "${xwin_fullname}/xwin"
    rm ${xwin_fullname}.tar.gz
fi

# Pass disable-symlinks, as we use a case insensitive filesystem as the destination
xwin_files="${PWD}/${xwin_version}/files"
${xwin_version}/xwin --accept-license --manifest-version 16 --cache-dir ${xwin_version}/xwin-cache  splat --disable-symlinks --include-debug-libs --output ${xwin_files}
