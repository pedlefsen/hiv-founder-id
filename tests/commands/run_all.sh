for f in test_*.sh; do
  echo "$f"
  bash "$f" -H
done
