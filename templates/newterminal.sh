#!/bin/sh

eval "$start"
conda activate $environment
conda info|grep "active environment"

#export SWOT_HYDROLOGY_TOOLBOX=/users/alchrist/documents/tools/swot-hydrology-toolbox
#export SWOT_HYDROLOGY_TOOLBOX=/users/alchrist/documents/Tools/cnes/swot-hydrology-toolbox-release_version_01_28_2022
export SWOT_HYDROLOGY_TOOLBOX=$SWOT_HYDROLOGY_TOOLBOX
export PYTHONPATH=$SWOT_HYDROLOGY_TOOLBOX/processing/:$PYTHONPATH
export PYTHONPATH=$SWOT_HYDROLOGY_TOOLBOX/processing/src/:$PYTHONPATH
export PYTHONPATH=$SWOT_HYDROLOGY_TOOLBOX/processing/src/cnes/sas:$PYTHONPATH
export PYTHONPATH=$SWOT_HYDROLOGY_TOOLBOX/sisimp/:$PYTHONPATH
export PYTHONPATH=$RIVEROBS/src:$PYTHONPATH
export RIVEROBS=/users/alchrist/documents/tools/RiverObs

$command

## pip uninstall shapely
## pip install shapely --no-binary :all:
