#!/bin/bash

nevents=1000
pos_z=-1000

for particle in "mu-" "e-" #"kaon-" #"e-" "mu-" "pi-"
do
    for energy in 50 #4 6 8 10 20 30 40 50 60 70 80 90 100 125 150 175 200
    do
        for pos_x in 0.1 82.5 #0.1 0.5 1.0 2.0 5.0 10.0
        do
            for pos_y in 0.1 82.5 #0.1 0.5 1.0 2.0 5.0 10.0
            do
                echo "Running simulation for particle ${particle} with energy ${energy} GeV at position (${pos_x}, ${pos_y}, ${pos_z}) mm"
		        ./generic_condor.sh $nevents $particle $energy $pos_x $pos_y $pos_z
            done
        done
    done
done
