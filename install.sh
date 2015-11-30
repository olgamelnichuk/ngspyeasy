#!/usr/bin/env bash

set -e

$(cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

DEST=/opt/ngspyeasy

mkdir -p ${DEST}
cp -R ./ngspyeasy/* ${DEST}

echo "#!/usr/bin/env sh
python ${DEST}/ngspyeasy.py \"\$@\"" > /usr/bin/ngspyeasy
chmod a+x /usr/bin/ngspyeasy

echo "#!/usr/bin/env sh
python ${DEST}/ngspyeasy_tool.py \"\$@\"" > /usr/bin/ngspyeasy_tool
chmod a+x /usr/bin/ngspyeasy_tool
