#!/bin/bash

usage()
{
cat << EOF
usage: $0 options

This script adds a header to all files matching the provided pattern in the given directory

OPTIONS:
   -h      Show this message
   -l      license file 
   -d      directory to search
   -p      file pattern to match
   -c      comment line prefix
EOF
}

license=
pattern=
directory=
char=
while getopts "hl:p:d:c:" OPTION
do
     case $OPTION in
         h)
             usage
             exit 1
             ;;
         l)
             license=$OPTARG
             ;;
         p)
             pattern=$OPTARG
             ;;
         d)
             directory=$OPTARG
             ;;
	 c)
	     char=$OPTARG
	     ;;
         ?)
             usage
             exit
             ;;
     esac
done
if [ ! -z "$char" ]; then
    sed -i "s/^/$char/" $license
fi

echo "Prepending ${license} to files with pattern ${pattern} in directory ${directory}"


for i in ${directory}/${pattern}
do
  if ! grep -q Copyright $i
  then
    cat ${license} $i >$i.new && mv $i.new $i #&& mv $i $i.bak && mv $i.new $i
  fi
done

if [ ! -z "$char" ]; then
    sed -i "s/^$char//" $license
fi
