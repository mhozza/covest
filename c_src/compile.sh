#!/bin/bash
gcc --std=c99 -shared -fPIC covest_lib.c -lm -o ../lib/covestlib.so
