sh tmp-time-test.sh 16 LOG - -
sh tmp-time-test.sh 16 LOG_ZERO - -
sh tmp-time-test.sh 16 TABLE - -
sh tmp-time-test.sh 16 TABLE LE,LAZY -
sh tmp-time-test.sh 16 SPLIT 16 4 ALTMAP,NOSSE -
sh tmp-time-test.sh 16 SPLIT 16 4 ALTMAP,LAZY,SSE -
sh tmp-time-test.sh 16 SPLIT 16 4 ALTMAP,LAZY,NOSSE -
sh tmp-time-test.sh 16 SPLIT 16 4 ALTMAP,SSE -
sh tmp-time-test.sh 16 SPLIT 16 4 NOSSE -
sh tmp-time-test.sh 16 SPLIT 16 4 LAZY,SSE -
sh tmp-time-test.sh 16 SPLIT 16 4 LAZY,NOSSE -
sh tmp-time-test.sh 16 SPLIT 16 4 SSE -
