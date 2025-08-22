%
% cartpole.m
%
% created on: 09.10.2015
%     author: pushpak
%
% see readme file for more information on the dcdc example
%
% you need to run ./cartpole binary first 
%
% so that the file: cartpole_controller.bdd is created
%

%function cartpole
clear set
close all


% load the symbolic set containing the controller
controller=StaticController('cartpole_sparse');

%% plot the dcdc safe set

% plot the domain of the controller
colors=get(groot,'DefaultAxesColorOrder');

% plot controller domain
dom=controller.domain;
plot(dom(:,1),dom(:,2),'.','color',0.6*ones(3,1))
hold on




