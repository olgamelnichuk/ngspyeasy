#!/usr/bin/env bash

$(cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

DEST=/opt/ngspyeasy

mkdir -p ${DEST}
cp -R ./ngspyeasy/* ${DEST}

echo "#!/usr/bin/env sh
python ${DEST}/ngspyeasy.py \"$@\"" > /usr/sbin/ngspyeasy
chmod a+x /usr/sbin/ngspyeasy

echo "#!/usr/bin/env sh
python ${DEST}/ngspyeasy_tool.py \"$@\"" > /usr/sbin/ngspyeasy_tool
chmod a+x /usr/sbin/ngspyeasy_tool
