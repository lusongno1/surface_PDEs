#! /bin/bash
./surfpatternfm >stdout.txt 2>&1 &
tail -f stdout.txt
