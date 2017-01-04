#!/bin/csh -f
#

set obsid = $1
set nm = `ls output/*$1*.tar |wc -l | awk '{print 378-$1}'`
echo $nm
set fil = f`echo $obsid | awk '{print int($1/1000)-2}'`
if ($nm < 25 ) then
    set files = ( `ls output/*_99992000*.tar` )
    foreach f ( $files )
        set nf = `echo $f | sed -e "s/2000/$obsid/"`
        set nf = `echo $nf | sed -e "s/f0/$fil/"`
        if ( ! -e $nf ) then 
          echo $nf
          set fid = `echo $nf | cut -d'/' -f2 | cut -d'.' -f1 | awk 'BEGIN{FS="_"} {printf("%s_%s_%s", $4, $5, $6)}'`
          set log = work/logs/raytrace_9999${obsid}_${fid}_0pbs.log
          if (-e $log ) cat $log
          echo
        endif
    end
endif
