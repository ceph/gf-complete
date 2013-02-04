if [ $# -lt 4 ]; then
  echo 'usage: sh tmp-test.sh w gf_specs (e.g. LOG - -)' >&2
  exit 1
fi

w=$1
shift
i=1024
while [ $i -le 134217728 ]; do
  iter=`echo $i | awk '{ print (134217728/$1)*1 }'`
  echo $i $iter $w $* `./gf_time $w G -1 $i $iter $* | head -n 3 | tail -n 2`
  i=`echo $i | awk '{ print $1*2 }'`
done

