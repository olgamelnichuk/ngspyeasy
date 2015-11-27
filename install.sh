#!/usr/bin/env bash

pip install docker-py

$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

DIR=/opt/ngspyeasy

mkdir -p ${DIR}
cp -R ./ngspyeasy/* ${DIR}

echo "#!/usr/bin/env sh
python ${DIR}/ngspyeasy.py \"$@\"" > /usr/sbin/ngspyeasy
chmod a+x /usr/sbin/ngspyeasy

echo "#!/usr/bin/env sh
python ${DIR}/ngspyeasy_tool.py \"$@\"" > /usr/sbin/ngspyeasy_tool
chmod a+x /usr/sbin/ngspyeasy_tool
