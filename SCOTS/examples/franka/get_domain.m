clear all
close all
clc
currentdir = pwd; 
parentdir = fileparts(currentdir);
parentdir = fileparts(parentdir);
addpath(parentdir)
path1 = fullfile(parentdir,'/mfiles');
addpath(path1)
path2 = fullfile(parentdir,'/mfiles/.');
addpath(path2)
path3 =  fullfile(path2,'/mexfiles');
addpath(path3)
path4 =  fullfile(path2,'/mexfiles/.');
addpath(path4)

controller=StaticController('scots');
% plot controller domain
dom=controller.domain;
[dom_r,dom_c]=size(dom);
for i=1:dom_r
    table(i,:)=[dom(i,:),controller.control(dom(i,:))];
end
domain_Table = array2table(table);
% T.Properties.VariableNames(1:length(v)) ="input_torque";
writetable(domain_Table,"franka.csv")
domain_data=readmatrix("franka.csv");
