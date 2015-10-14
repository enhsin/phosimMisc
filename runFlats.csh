#!/bin/csh -f

set myip = `hostname -i`
set machine = `hostname | cut -d'.' -f1`
set dt=`date`
set log = ~/flats_v3.5.log

#find dag finished with status 1 and 2
set filter = ( u g r i z y )
foreach f ( `find work/ -maxdepth 1 -name '*dagman.out' -cmin +20` )
    set statusnum = `tail -1 $f | awk '{print $NF}'`
    set obsid = `echo $f | cut -d'_' -f2 | cut -d'.' -f1`
    set submitHost = `grep "private command socket" $f | cut -d'<' -f2 | cut -d':' -f1 | tail -1`
    if ( $myip != $submitHost && $statusnum != 0 ) then
        continue
    endif
    if ( $statusnum != 0 ) then
       echo "$dt ${machine}: resubmit $obsid" >> $log
       set fnum = `echo $obsid | awk 'BEGIN{FS=""}{print $5-1}'`
       set opt = `echo $obsid | awk 'BEGIN{FS=""}{print $8}'`
       set oldFiles = `find work/ -maxdepth 1 -name 'lsst*'$obsid'*.gz' -cmin +7200`
       #files too old (probably corrupted)
       if ( $#oldFiles > 0 ) then
           foreach ff ( $oldFiles )
               echo "remove $ff"
               rm $ff
           end
       endif
       ./cleanlog $obsid
       python phosim.py examples/flats/flat${filter[$fnum]}_instcat_${opt} -g condor --checkpoint=110
       sleep 300
    else
        echo "$dt ${machine}: cleanup $obsid" >> $log
       ./cleanup $obsid
    endif
end

# dag with failed jobs or jobs that run too long
foreach f ( `find work/ -maxdepth 1 -name '*dagman.out' -cmin -30` )
    if ( `tail -1 $f | grep "EXITING WITH STATUS" | wc -l` > 0 ) then
        continue
    endif
    set obsid = `echo $f | cut -d'_' -f2 | cut -d'.' -f1`
    set submitHost = `grep "private command socket" $f | cut -d'<' -f2 | cut -d':' -f1 | tail -1`
    if ( $myip != $submitHost ) then
        continue
    endif
    set line = `grep -n Failed $f | tail -1 | cut -d: -f1 | awk '{print $1+2}'`
    set failJob = `head -$line $f | tail -1 | awk '{print $NF}'`
    set lastLogTime = `grep "seconds since last log event" $f | tail -1 | awk '{print int($3/60/60)}'`
    echo $obsid $failJob $lastLogTime
    if ( $failJob > 10 || $lastLogTime > 10 ) then
        echo "$dt ${machine}: remove $obsid" >> $log
        set prid = `grep exec $f | head -1 | sed -e "s/_exec./ /" | awk '{print $5}'`
        condor_rm $prid
        sleep 300
    endif
end


@ maxjobs = 3000
set filesys = `echo $machine | cut -d'-' -f1 | awk '{if($1=="hansen" || $1=="new") print "hansen"; else print "rossmann"}'`

myquota >  ~/myquota.dat
set availfile = `grep $filesys ~/myquota.dat | awk '{print $6, $7}' | awk '{ gsub(/','/, ""); print $2-$1}'`
if ( $availfile == '' ) then
    set filesys = `echo $machine | cut -d'-' -f1 | awk '{if($1=="hansen" || $1=="new") print "lustreC"; else print "lustreA"}'`
    set availfile = `grep $filesys ~/myquota.dat | awk '{print $6, $7}' | awk '{ gsub(/','/, ""); print $2-$1}'`
endif

if ( $availfile < 99 ) then
    echo "$dt ${machine}: Disk Full. left $availfile k" >> $log
    exit
endif

condor_q $user | grep running >& ~/${machine}.run
set conerr = `wc -l < ~/${machine}.run`
if ( $conerr == 0 ) then
    echo "$dt ${machine}: Condor Read Error" >> $log
        exit
    endif
endif

set idlejob = `cat ~/${machine}.run | awk '{print $7}'`
set totaljob = `cat ~/${machine}.run | awk '{print $1}'`

if ( $idlejob > 200 || $totaljob > $maxjobs ) then
    echo "$dt ${machine}: too many jobs" >> $log
    exit
endif


set doneID = `cat done.lis`
@ sub = 0
foreach fil ( u )
    #foreach fil ( u g r i z y )
    foreach i ( 0 1 2 3 4 5 6 7 8 9 )
        @ sub = 0
        set obsid = `grep obshist examples/flats/flat${fil}_instcat_${i} | awk '{print $2}'`
        foreach d ( $doneID )
            if ( $d == $obsid ) then
                @ sub = 1
                break
            endif
        end
        if ( $sub == 1 ) continue
        echo "submit $obsid"
        python phosim.py examples/flats/flat${fil}_instcat_${i} -g condor --checkpoint=110
        @ sub = 2
        echo $obsid >> done.lis
        break
    end
    if ( $sub == 2 ) break
end
