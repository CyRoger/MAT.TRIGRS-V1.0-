%main program
clc
clear
T=8000;%Rainfall duration
file_path = 'C:\Users\DELL\Desktop\TRIGRSshuicheng\rainfall_tif\'; %Rainfall data storage path
t = textread('t.txt','%f');%Input time of each rainfall event
dem = GRIDobj('DEM.tif');
slope= GRIDobj('SLOPE_angle.tif');
flowdirection= GRIDobj('flow_direction.tif');
zmax = GRIDobj('Soil_depth.tif');
depthwt = GRIDobj('depthwater.tif');
Ys = GRIDobj('Weight.tif');
Yw = GRIDobj('water_weight.tif');
c = GRIDobj('Cohesion.tif');
f = GRIDobj('friction.tif');
Ks = GRIDobj('Ks.tif');
Izlt = GRIDobj('IZlT.tif');
D0 = GRIDobj('D0.tif');

[Phead,ZMAX,Fs]=TRIGRS(T,file_path,dem,slope,flowdirection,zmax,depthwt,Ys,Yw,c,f,Ks,Izlt,D0);
w1=dem.size(1);
w2=dem.size(2);
dem.Z=Phead;
GRIDobj2geotiff(dem,'C:\Users\87856\Desktop\shuicheng\shuicheng\results\wet1_Phead');%Result storage location
dem.Z=Fs;
GRIDobj2geotiff(dem,'C:\Users\87856\Desktop\shuicheng\shuicheng\results\wet1_Fs'); 
dem.Z=ZMAX;
GRIDobj2geotiff(dem,'C:\Users\87856\Desktop\shuicheng\shuicheng\results\wet1_ZMAX11');
