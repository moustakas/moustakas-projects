#!/bin/tcsh
# echo "ages_mgfe_isedfit, /sdss, /isedfit, /kcorrect, /cl" | /usr/bin/nohup idl >& /moustakas/archive/projects/ages/mgfe/logs/chunk0-7.log &

# echo "ages_mgfe_isedfit, /sdss, /kcorrect, /cl" | /usr/bin/nohup idl >& /moustakas/archive/projects/ages/mgfe/logs/sdss.kcorr.log &

echo "ages_mgfe_isedfit, /write_param, /build_grids, /model_phot, /isedfit, /qaplot_sed, /kcorrect, /cl" | /usr/bin/nohup idl >& /moustakas/archive/projects/ages/mgfe/logs/ages.log &

