#!/bin/bash
# This file is part of MOM6, the Modular Ocean Model version 6.
# See the LICENSE file for licensing information.
# SPDX-License-Identifier: Apache-2.0

for chk in $1 $2; do
    awk '{print $(NF-2) " " $(NF-1) " " $(NF),$0}' ${chk} | sort > ${chk}.sorted
done

diff $1.sorted $2.sorted | grep "^[<>]" | grep "^>" > /dev/null

if [ $? -eq 0 ]; then
    diff $1.sorted $2.sorted | head -n 100
    exit 1
fi
