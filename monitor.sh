#!/bin/bash
while true; do
  clear
  date
  ps aux | grep 49262 | grep -v grep || { echo "PROCESSO TERMINATO"; tail -50 test_1000ast_fix_*.log; break; }
  echo "---"
  ls -lh test_1000ast_fix_*.log
  echo "---"
  tail -10 test_1000ast_fix_*.log
  sleep 30
done
