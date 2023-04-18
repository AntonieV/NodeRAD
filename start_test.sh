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

source ~/miniconda3/etc/profile.d/conda.sh

conda activate snakemake

snakemake --use-conda --conda-frontend mamba --cores $(nproc) --directory .test --show-failed-logs --conda-cleanup-pkgs cache -R all

snakemake --report --directory .test

cd .test

xdg-open report.html
xdg-open results/qc/multiqc/multiqc.html
cd ..
