 #!/usr/bin/env bash 


[ ! -z "*.raw.*.csv" ] && rm  rep*.raw.5mer.csv
[ ! -z "*delta*5mer.csv" ] && rm *delta.*.csv
[ ! -z "*sumErr*csv" ] && rm *sumErr*csv

rm *.MODEL.*
