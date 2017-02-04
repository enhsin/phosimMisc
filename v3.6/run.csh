#!/bin/csh -f

set fil = ( u g r i z y )
set obsid = ( 2 3 4 5 6 7 )
set ckpt = ( 2 7 7 6 5 5 )

@ i = 1
while ( $i <= 6 )
    foreach vid ( 5 6 7 8 9 )
        if ( -e work/dag_9999${obsid[$i]}00${vid}.dag ) then
            @ n = `ls output/*9999${obsid[$i]}00$vid*.tar |wc -l`
            echo "9999${obsid[$i]}00$vid $n"
            if ( $n == 378 ) then
                ./cleanup 9999${obsid[$i]}00$vid
            else
                if ( `find work/ -maxdepth 1 -name dag_9999${obsid[$i]}00${vid}.dag -cmin +2160 |wc -l` == 1 ) then
                    @ nRec = `find output/ -maxdepth 1 -name "*9999${obsid[$i]}00$vid*.tar" -cmin -300 | wc -l`
                    echo "9999${obsid[$i]}00$vid $nRec"
                    if ( $nRec == 0 ) then
                        set failed = `ls work/logs/raytrace*9999${obsid[$i]}00$vid*pbs.log | xargs grep Terminated | cut -d':' -f1 | cut -d'/' -f3 | cut -d'p' -f1`
                        foreach f ( $failed )
                            echo $f
                            rm work/logs/${f}pbs.log
                            @ ckpt = `echo $f | cut -d'_' -f6`
                            @ ckpt1 = $ckpt
                            set f1 = $f
                            while ( -e work/$f1.pbs )
                                perl -pi -e 's/standby/physics/g' work/$f1.pbs
                                perl -pi -e 's/4:00/6:00/g' work/$f1.pbs
                                if ( $ckpt1 == $ckpt ) then
                                    grep -v depend work/$f1.pbs > tmp.pbs
                                    mv tmp.pbs work/$f1.pbs
                                else
                                   set pid0 = `grep depend work/$f1.pbs | cut -d':' -f2`
                                   perl -pi -e "s/$pid0/${pid}/g" work/$f1.pbs
                                endif
                                set pid = `qsub -e work/errors/${f1}.pbs.err -o work/logs/${f1}pbs.log work/$f1.pbs | cut -d'.' -f1`
                                echo $pid $f1
                                @ ckpt1++
                                set f1 = `echo $f | sed -e "s/_${ckpt}/_${ckpt1}/"`
                            end
                        end
                        #./cleanup 9999${obsid[$i]}00$vid
                    endif
                endif
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
            #@ nFree = `qlist | grep physics | awk '{print $5}'`
            #if ( $nFree > 16 ) then
            #    foreach pid ( `qstat -u $user | grep Q | cut -d'.' -f1` )
            #        qmove physics $pid
            #        @ nFree = $nFree - 4
            #        if ( $nFree < 16 ) break
            #    end
            #endif
            exit
        endif
    end
    @ i++
end

