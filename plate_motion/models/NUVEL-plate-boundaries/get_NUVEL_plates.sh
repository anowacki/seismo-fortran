#!/bin/sh

plates="African Antarctic Arabian Australian Caribbean Cocos Eurasian Indian Juan Nazca North_Am Pacific Philippine Scotia South_Am"

for p in $plates; do 
	curl -s -f -o $p.txt http://jules.unavco.org/GMT/$p.txt
done
