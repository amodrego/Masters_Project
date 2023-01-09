#!/bin/bash
# Put the config file in the location where project3D expects to find it.

mkdir -p config
mv PhysiCell_settings.xml config/
# Create the results folder for project3D.
mkdir -p output


while [ "$1" != "" ]; do
    # All parameters appearing in python and HTCondor scripts
    case $1 in
        --odens)
            odens=$2
            shift 2
        ;;

        --cdens)
            cdens=$2
            shift 2
        ;;

        --gdens)
            gdens=$2
            shift 2
        ;;

        --ldens)
            ldens=$2
            shift 2
        ;;

        *)
            echo "command not recognized: $1, $2"
            exit
           # unknown option
        ;;
    esac
done

suffix=$(printf "o%05.2f_c%05.2f_g%05.2f" ${odens} ${cdens} ${gdens})
results_path=results_${suffix}

python data_setup.py ${odens} ${cdens} ${gdens} ${ldens}


# (change project3D to the name of the PhysiCell project you want to run)
./project3D

mkdir -p ${results_path}

cp -r *.txt ${results_path}
cp config/PhysiCell_settings.xml ${results_path}
cp -r output/. ${results_path}
rm -r output

# Clean and compress results to improve transfer speeds.
tar -cjf ${results_path}.tar.bz2 ${results_path}
mv ${results_path}.tar.bz2 results/

rm -r ${results_path}