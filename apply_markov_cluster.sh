#!/bin/bash

cluster_dir=/tmp/qvz_markov_clusters
rm -r $cluster_dir
mkdir $cluster_dir
cluster_file_root="$cluster_dir/cluster"
COUNTER=0
while read cluster_id ; do
    cluster_file="$cluster_file_root$cluster_id"
    COUNTER=$((COUNTER + 1))
    tail -n+$COUNTER $2 | head -n1 >> $cluster_file
    echo $cluster_file
done < "$1"
