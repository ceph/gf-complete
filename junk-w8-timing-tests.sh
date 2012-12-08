sh tmp-time-test.sh 8 LOG - -
sh tmp-time-test.sh 8 LOG_ZERO - -
sh tmp-time-test.sh 8 TABLE - -
sh tmp-time-test.sh 8 TABLE DOUBLE -
sh tmp-time-test.sh 8 TABLE DOUBLE,LAZY -
sh tmp-time-test.sh 8 BYTWO_p - - 
sh tmp-time-test.sh 8 BYTWO_b - - 
sh tmp-time-test.sh 8 BYTWO_p SSE - 
sh tmp-time-test.sh 8 BYTWO_b SSE - 
sh tmp-time-test.sh 8 SPLIT 8 4 NOSSE -
sh tmp-time-test.sh 8 SPLIT 8 4 SSE -
sh tmp-time-test.sh 8 COMPOSITE 2 4 TABLE SINGLE,SSE - - -
sh tmp-time-test.sh 8 COMPOSITE 2 4 TABLE SINGLE,SSE - ALTMAP -
