#!/usr/bin/env bash

var1=$(cat ./out_core_pip_lkfunc.estimates.json | grep LogLikelihood | awk -F":" '{ print $2}' | sed 's/\"//g' | sed 's/ //g');
var2=-21.354418417018628;
status=$(echo $var1'=='$var2 | bc -l);

if [ "$status" -eq "1" ]; then
   exit 0;
else
   exit 1;
fi