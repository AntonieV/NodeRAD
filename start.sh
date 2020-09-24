#!/bin/bash

echo -e ""

echo -e "https://github.com/AntonieV/NodeRAD.git               ▒ "
echo -e "by Antonie Vietor                                     ▓ "
echo -e "                                                     ░█ "
echo -e "                                                   ░░██ "
echo -e "                                    ░           ▒▓█████ "
echo -e "                                   ▒█░░           ░▒▓██ "
echo -e "                                   ░████████▓▒      ░██ "
echo -e "                                     ██▓▓▒░░      ░▒██  "
echo -e "                                      ██▒    ▒▓██████   "
echo -e "                                       ███▓    ░▒███    "
echo -e "                                         ████ ▒███░     "
echo -e "                                           ░███░        "
echo -e "                                        ░██   ░███░     "
echo -e "                                      ░███       ░██    "
echo -e "███╗   ██╗ ██████╗ ██████╗ ███████╗  ██▓        ░░▒██░  ██████╗  █████╗ ██████╗ "
echo -e "████╗  ██║██╔═══██╗██╔══██╗██╔════╝ ██▒   ▒▒▓█████████░ ██╔══██╗██╔══██╗██╔══██╗"
echo -e "██╔██╗ ██║██║   ██║██║  ██║█████╗  ▒█▓░       ░░▒▒▓▓██▒ ██████╔╝███████║██║  ██║"
echo -e "██║╚██╗██║██║   ██║██║  ██║██╔══╝  ▒██▓▒░░░        ░██▒ ██╔══██╗██╔══██║██║  ██║"
echo -e "██║ ╚████║╚██████╔╝██████╔╝███████╗ █████████▓▒▒   ▄█▀░ ██║  ██║██║  ██║██████╔╝"
echo -e "╚═╝  ╚═══╝ ╚═════╝ ╚═════╝ ╚══════╝ ░██▓▓▒▒░░░░   ███░  ╚═╝  ╚═╝╚═╝  ╚═╝╚═════╝ "
echo -e "                                      ███▒       ███    "
echo -e "                                        ████    ██▀     "
echo -e "                                           █████        "
echo -e "                                         ░█▓ ░████      "
echo -e "                                       ░██▓     ░██░    "
echo -e "                                      ██████▓▒    ░██   "
echo -e "                                     ██▓▒▒░        ░██  "
echo -e "                                    ░█░     ▒▒▓███████░ "
echo -e "                                    ██▒░       ░░▒▒▓▓█▒ "
echo -e "                                    ██████▓▒▒        ░  "
echo -e "                                    █▓▓▒░░              " 
echo -e "                                    ▓░                           NodeRAD v1.0   "
echo -e "                                    ▒                   "

echo -e ""
echo -e "Loading...";sleep 3
#TODO: init data and config directory

#TODO: if exists
source ~/miniconda3/etc/profile.d/conda.sh

#TODO: else which conda
# grep miniconda3 from /home/xxx/miniconda3/condabin/conda as $con="miniconda3", e.g. anaconda3
    #if exists
    #source ~/$con/etc/profile.d/conda.sh
    #else install conda

#TODO if conda env snakemake not exists install snakemake
#conda create -c bioconda -c conda-forge -n snakemake snakemake

conda activate snakemake

#TODO: if testing data wanted
#cd .test
#git submodule add -f https://github.com/snakemake-workflows/ngs-test-data.git
#git commit -m "add ngs-test-data submodule to get GitHub Actions tests to pass"
#cd ..

#$(nproc) eventuell als user-abfage
##TODO:do you want to execute the complete workflow (press c + Enter), make up or repeat missing steps of the workflow (press s + Enter) or execute a specific step (press r + Enter)?
##TODO:if complete workflow (input c):
###DEBUG-MODE###
#snakemake --use-conda --cores $(nproc) --directory .test --show-failed-logs --printshellcmds --verbose -R all
###USER-MODE###
snakemake --use-conda --cores $(nproc) --directory .test --show-failed-logs -R all


##TODO:else if missing steps (input s):
#snakemake --use-conda --cores $(nproc) --directory .test --show-failed-logs

##TODO:else show list of workflow rules
#please select a specific rule and press Enter: {userinput}
#snakemake --use-conda --cores $(nproc) --directory .test --show-failed-logs -R rule_{userinput}

#do you want to create a report.html file? y/n
snakemake --report --directory .test
#do you want to view report.html file? y/n
#if input = y | yes
cd .test ##TODO:add directory
#sensible-browser report.html
xdg-open report.html # TODO : intgration of multiqc.html
xdg-open results/qc/multiqc/multiqc.html
cd ..
