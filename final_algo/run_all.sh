for i in data/*/; do 
  python3 parallel_merged_graph.py ${i%%/} 16 > ans/${i%%/}.log &
done
