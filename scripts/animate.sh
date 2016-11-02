#!/bin/bash

sm2 -config /tmp/csm_config.txt -file_jj /tmp/csm_journal.txt < /tmp/csm_in.log > /tmp/csm_out.log

json_extract -nth 0 < /tmp/csm_journal.txt > /tmp/csm_matching.txt

$ROS_DIR/workspace/src/MSCKF_2.0/postprocessing/sm_animate \
         -write_info 1 \
         -ref_countour_color '#00f' -sens_countour_color '#f00' \
         -ref_countour_width 0.01 -sens_countour_width 0.01 \
         -ref_rays_draw 0 -sens_rays_draw 0 \
         -ref_pose_radius 0.04 \
         -ref_horizon 65 -sens_horizon 65 \
         -ref_connect_threshold 0.3 -sens_connect_threshold 0.3 \
         -in /tmp/csm_matching.txt -out /tmp/csm_animation.pdf

evince /tmp/csm_animation.pdf
