%% Header
% Title: KBD101_test.m (modified from KDC101.m)
% Created Date: 2026-08-14
% Last modified date: 2026-09-03
% Matlab Version: R2025b
% Thorlabs DLL version: Kinesis 1.14.44
%% Notes:
%
% Modified KDC101 using the PRM1/M-Z8 stage code for KBD101 and DDR25/M

%% Start of code
% clear all; close all; clc

%% Add and Import Assemblies
devCLI = NET.addAssembly('C:\Program Files\Thorlabs\Kinesis\Thorlabs.MotionControl.DeviceManagerCLI.dll');
genCLI = NET.addAssembly('C:\Program Files\Thorlabs\Kinesis\Thorlabs.MotionControl.GenericMotorCLI.dll');
motCLI = NET.addAssembly('C:\Program Files\Thorlabs\Kinesis\Thorlabs.MotionControl.KCube.BrushlessMotorCLI.dll');

import Thorlabs.MotionControl.DeviceManagerCLI.*
import Thorlabs.MotionControl.GenericMotorCLI.*
import Thorlabs.MotionControl.GenericMotorCLI.ControlParameters.*
import Thorlabs.MotionControl.GenericMotorCLI.AdvancedMotor.*
import Thorlabs.MotionControl.GenericMotorCLI.KCubeMotor.*
import Thorlabs.MotionControl.GenericMotorCLI.Settings.*
import Thorlabs.MotionControl.KCube.BrushlessMotorCLI.*

%% Connect
%Build Device List loads the connected devices to available memory
DeviceManagerCLI.BuildDeviceList();

% Update serial number to correct device
serialNumber = '28254011'; % this is KCube closest to the door
timeout_val = 60000;

% Connect to the controller
device = KCubeBrushlessMotor.CreateKCubeBrushlessMotor(serialNumber);

device.Connect(serialNumber);
%% Move
try
    % Try/Catch statement used to disconnect correctly after an error

    device.WaitForSettingsInitialized(5000);
    device.StartPolling(250);
    
    %Pull the enumeration values from the DeviceManagerCLI
    optionTypeHandle = devCLI.AssemblyHandle.GetType('Thorlabs.MotionControl.DeviceManagerCLI.DeviceSettingsSectionBase+SettingsUseOptionType');
    optionTypeEnums = optionTypeHandle.GetEnumValues(); 
    
    %Load Settings to the controller
    motorConfiguration = device.LoadMotorConfiguration(serialNumber);
    motorConfiguration.LoadSettingsOption = optionTypeEnums.Get(1); % File Settings Option
    motorConfiguration.DeviceSettingsName = 'DDR25'; %The actuator type needs to be set here. This specifically loads an PRM1-Z8
    
    factory = KCubeMotor.KCubeBrushlessMotorSettingsFactory();
    device.SetSettings(factory.GetSettings(motorConfiguration), true, false);
    
    % Get/set velocity parameters (doesn't affect home)
    velParams   = device.GetVelocityParams();
    velParams.MaxVelocity = 360; % range is 360 to 1800 degrees/s
    device.SetVelocityParams(velParams);
    
    % Set relative move distance
    relMove = -180; % degrees
    device.SetMoveRelativeDistance(relMove);
    
    % Enable the device and start sending commands
    device.EnableDevice();
    pause(1); %wait to make sure the cube is enabled (seconds)
    
    % Home the stage
    fprintf("Homing...\n")
    device.Home(timeout_val);
    fprintf("Homed\n\n")
    
    % Move the stage to absolute value
    move1 = 145;  % degrees
    fprintf(['Moving to ' num2str(move1) ' degrees...\n'])
    device.MoveTo(145,timeout_val);
    fprintf("Moved\n")
    
    pause(5); % wait (seconds)
    
    % Move the stage again
    move2 = 210; degrees
    fprintf(['Moving to ' num2str(move2) ' degrees...\n'])
    device.MoveTo(move2,timeout_val);
    fprintf("Moved\n")
    
    pause(8); % wait (seconds)
    
    % Move the stage a relative amount
    fprintf(['Moving to ' num2str(move3) ' degrees...\n'])
    device.MoveRelative(timeout_val);
    fprintf("Moved\n")
    
    pause(13); % wait (seconds)
    
    % Home the stage
    fprintf("Homing...\n")
    device.Home(timeout_val);
    fprintf("Homed\n\n")
    
catch e
    fprintf("Error has caused the program to stop, disconnecting..\n")
    fprintf(e.identifier);
    fprintf("\n");
    fprintf(e.message);
end

%% Disconnect from controller
% comment these out if you want to test things in command line
device.StopPolling();
device.Disconnect();