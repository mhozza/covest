#!/bin/bash
for i in $@; do
    TMP=`tempfile`;
    sed '1d' $i > $TMP; mv $TMP $i;
done
