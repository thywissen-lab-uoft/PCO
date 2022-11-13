# PCO
MATLAB code for running the Pixefly QE PCO cameras. This code is not designed to run on different cameras,
but it could be feasibly altered by someone proficient in the drivers.


### Hardware Requirements
To get the camera initially working with a Windows 10 machine,
you will require the following:

- Windows 10 OS
- Pixelfly QE Ethernet Camera
- PCO 540 Board
- Magma PCI to PCIE Adapter

The Magma board is a PCIE adapter box that allows the PCO 540 Board 
which is a PCI card to interface with modern motherboards
which lack a PCI slot.

The Magma PCI to PCIE chasis box is a plug and play device.  
Use the following steps to power cycle the box.

1) Turn off both Magma and PC
2) Turn on Magma with camera plugged in
3) Turn on PC

### Initial Communication
The drivers for the pco 540 camera are provided by PCO at 

https://www.pco.de/support/interface/scientific-cameras-1/pixelfly-qe/

With the pertinent zip file

DI_540_W7_W8_W10_V201_12.zip

Use the above zip file to install the drivers.  After installation,
the pco540 board should appear under your Device Manager under
PCO cameras/pco540

### PCO Camware
PCO provides a GUI to operate the camera.  As of writing of this guide,
version above 3.17 DO NOT work with the camera.  This guide has used
version 3.17.  The software may be found at 

https://ca.pco-tech.com/software/camera-control-software/pcocamware/

SW_CAMWAREWIN64_317  

You should verifying that you can take images with the camera from here.

### Operating the code
Once all drivers have been installed, you may run the gui and the analysis programs

analysis/PCO_analysis.m
imaq_gui/PCO_gui.m
