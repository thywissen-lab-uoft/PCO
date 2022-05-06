This code controls the rotating waveplate:

To use :

(1) Use rototation_stage_home to send to waveplate to it's home position and to
listen to digital commands.
(2) Ue rotation_stage_move(steps) to make the waveplate properly aligned 
(aligning the top sharpee markers). This changes due to hysteris over time.
(3) Redefine home using rotation_stage_home
(4) Allow for analog control with rotation_stage_analog