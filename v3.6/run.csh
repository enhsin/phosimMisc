#!/bin/csh -f

set fil = ( u g r i z y )
set obsid = ( 2 3 4 5 6 7 )
set ckpt = ( 2 5 5 4 3 3 )

@ i = 1
while ( $i <= 6 )
    foreach vid ( 5 6 7 8 9 )
        if ( -e work/dag_9999${obsid[$i]}00${vid}.dag ) then
            @ n = `ls output/*9999${obsid[$i]}00$vid*.tar |wc -l`
            if ( $n == 378 ) then
                ./cleanup 9999${obsid[$i]}00$vid
            endif
        endif
        @ hasSub = `grep 9999${obsid[$i]}00$vid submit.dat |wc -l`
        if ( $hasSub == 0 ) then
            echo 9999${obsid[$i]}00$vid
            if ( ! -e work/dag_9999${obsid[$i]}00${vid}.dag ) then
                python phosim.py examples/flats/flat${fil[$i]}_instcat_${vid} -t 4 -g condor --checkpoint=$ckpt[$i]
            endif
            python tools/cluster_submit_v0.py dag_9999${obsid[$i]}00${vid}.dag -w /scratch/conte/e/epeng/phosim_release/work -o /scratch/conte/e/epeng/phosim_release/output
            foreach pid ( `qstat -u $user | grep trim | cut -d'.' -f1` )
                qalter $pid -l walltime=00:00:10 -l mem=1GB
            end
            exit
        endif
    end
    @ i++
end

