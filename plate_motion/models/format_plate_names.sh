#!/bin/sh

cat NNR-MORVEL56-tab-delimited.text |\
	awk -F'\t' 'NR>=7 && NR<=63 && NF>1 \
	{printf("\"%-30s\", ",substr($1,4)); \
	if((NR%2==0&&NR<=31)||(NR%2==1&&NR>31)) printf(" &\n\t\t")}END{print ""}'