#!/bin/bash
#$ -cwd
#$ -N RSYNC
#$ -S /bin/bash
#$ -b y
#$ -l mem_requested=4G
#$ -l h_vmem=4G
#$ -l tmp_requested=1G
#$ -l fast_dm=10
#$ -pe smp 1

#############################################################################################
#                                                                                         
# copy data from fridge to cluster storage                                                
#                                                                                         
# setting up: first, you have to setup the ssh key-based access from cluster to gridION      
#		1. generate a ssh key if you haven't done so									  
#				ssh-keygen                                                                
#		2. start a session on a node with 10Gbps access on clus2020ter that has access to gridION
#				qrsh -l fast_dm=10                                                        
#		3. copy the ssh key to grdION 														  
#				ssh-copy-id gtgadmin@129.94.120.250 -p 2020                                 
#		4. check if password-less access work											  
#				ssh  gtgadmin@129.94.120.250 -p 2020                                            
#		5. logout from gridION and then from the qrsh session                                                 
#				ctrl+d twice                                                              
#                                                                                         
#                                                                                         
# usage: qsub ./grid_dload.sge.sh /path/on/gridion /path/on/cluster/storage               
#                                                                                         
#############################################################################################

SOURCE=$1
DEST=$2

rsync -a --progress --rsh='ssh -T -o Compression=no -x -p 2020' hasindu@129.94.120.250:${SOURCE} ${DEST}


