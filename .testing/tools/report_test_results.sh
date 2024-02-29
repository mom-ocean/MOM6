#!/bin/sh
RESULTS=${1:-${PWD}/results}

GREEN="\033[0;32m"
RESET="\033[0m"
PASS="${GREEN}PASS${RESET}"

if [ -d ${RESULTS} ]; then
  if ls ${RESULTS}/*/std.*.err &> /dev/null; then
    echo "The following tests failed to complete:"
	ls ${RESULTS}/*/std.*.out \
      | awk '{ \
        split($$0,a,"/"); \
        split(a[length(a)],t,"."); \
        v=t[2]; \
        if(length(t)>4) v=v"."t[4]; print a[length(a)-1],":",v}'
  fi

  if ls ${RESULTS}/*/ocean.stats.*.diff &> /dev/null; then
    echo "The following tests report solution regressions:"
    ls ${RESULTS}/*/ocean.stats.*.diff \
      | awk '{ \
        split($$0,a,"/"); \
        split(a[length(a)],t,"."); \
        v=t[3]; \
        if(length(t)>4) v=v"."t[4]; print a[length(a)-1],":",v}'
  fi

  if ls ${RESULTS}/*/chksum_diag.*.diff &> /dev/null; then
    echo "The following tests report diagnostic regressions:"
    ls ${RESULTS}/*/chksum_diag.*.diff \
      | awk '{ \
        split($$0,a,"/"); \
        split(a[length(a)],t,"."); \
        v=t[2]; \
        if(length(t)>4) v=v"."t[4]; print a[length(a)-1],":",v}'
  fi

  exit 1
else
  printf "${PASS}: All tests passed!\n"
fi
