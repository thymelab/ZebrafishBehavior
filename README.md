The first step is to track the high-speed movies, such as the following (track.slurm file):

#!/bin/bash
SLURM info
#module load Anaconda3/5.3.1
xvfb-run -a python /data/project/thymelab/Newest_Behavior_Scripts/tracking/highspeedmovieanalysis.py -r "rois_string" -m "hsmovieTue, Feb 15, 2022_$SLURM_ARRAY_TASK_ID.avi" -e "fulltestrun_final_01_27_2020"

The following script checks if the tracking is complete on these files

../../../scripts/behavior_published2020/prepandsupportfiles/queuecheck_behavior.py

Once complete, get genotype info ready:

../../../scripts/behavior_published2020/prepandsupportfiles/splitgenotypes_flip_updated.py

Inside the genotypes file, put a * before the genotypes to compare. The above file also may need to be edited if we are not using the same naming systems with het/hom/wt.

Then, make all of the input files with the following:
../../../scripts/behavior_published2020/prepandsupportfiles/makeslurmfiles_updated.py

This script needs to be edited if we are comparing different types of samples or different genotypes than we typically use.

a jobsubmission.sh file is created, which sends out all of the comparisons based on the genotypes file.

The generated submission scripts look something like this:
#!/bin/bash
SLURM info
module load Anaconda3
cd outputfulldata_GENERATEDFILENAME
python processmotiondata.py -t "../testlog.timestamp1.Tue, Feb 15, 2022" -e "../fulltestrun_final_01_27_2020" -c "../testlog.centroid1.Tue, Feb 15, 2022" -d "../testlog.motion1.Tue, Feb 15, 2022" -m "../hsmovieTue, Feb 15, 2022_" -g "../jagn1bc2_het_vs_jagn1bc2_hom_scripted_inputgenotypeids" -s "../sectionsfile" -j "../../../scripts/behavior_published2020/prepandsupportfiles/PlotParameters"

The inputgenotypeids are generated, the centroid/motion/timestamp files are from the lab view, the hsmovie files are from the tracking above. We provide the sectionsfile, PlotParameters, and fulltestrun_final_01_27_2020. The fulltestrun file is what was originally run with the labview code.
