#!/bin/bash

for chk in $1 $2; do
    awk '{print $(NF-2) " " $(NF-1) " " $(NF),$0}' ${chk} | sort > ${chk}.sorted
done

cmp $1.sorted $2.sorted

if [ $? -eq 1 ]; then
    diff $1.sorted $2.sorted | head -n 100
    exit 1
fi
