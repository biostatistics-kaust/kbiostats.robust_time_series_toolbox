#!/bin/bash

RIOTC="node ../../../../node_modules/riot-cli/lib/index.js"

function compile_riot(){
  file=$1
  cp $file ${file%.*}.tag
  $RIOTC ${file%.*}.tag ${file%.*}.js
  rm ${file%.*}.tag
}

destfile=components.js
for file in *.html
do
  echo riot $file ${file%.*}.js
  if [[ -e ${file%.*}.js ]]; then
    if [[ $file -nt ${file%.*}.js ]]; then
      compile_riot "$file"
    else
      echo "   Skipping file: file up to date."
    fi
  else
    compile_riot "$file"
    echo "require('./${file%.*}.js');" >> $destfile
  fi
done

cp $destfile 0$destfile
cat 0$destfile | sort | uniq > $destfile
rm 0$destfile
