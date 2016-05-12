::example batch file for windows; make sure gcc and opencv libaries can be found, e.g. put respective "bin" folders in PATH variable or copy (under windows) "libgcc_s_dw2-1.dll", "libstdc++-6.dll", "libopencv_calib3d310.dll" etc next to executable
::imageprocessing.exe infolder filenametemplate intstart intstop subrect_xleft subrect_ytop subrect_width subrect_height outfolder
bin\Release\imageprocessing.exe "" "img_%%09d_00-BF_EGFP_000.tif" 0 2 512 325 512 512 ""
pause