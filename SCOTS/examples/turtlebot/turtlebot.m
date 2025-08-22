%
% dcdc.m
%
% created: Jan 2017
%  author: Matthias Rungger
%
% see readme file for more information on the safety example
%
% you need to 1. have the mexfiles compiled 
%             2. run the ./dcdc binary first 
%
% so that the file: controller.scs is created
%


% addpath(parentdir);
function dcdc
clear set
close all
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

% initial state
x0=[1.5 0.5 0];
% x0=[0.2 9.4 0.4];
tau_s=0.3;
% load controller from file
controller=StaticController('scots_lab');  % tm_turtlebot_controller_with_gamma
%tm_turtlebot_controller_without_gamma
% simulate closed loop system
y=x0;
v=[];
loop=33; %43
total_time=loop*tau_s;
% plot controller domain
dom=controller.domain;
figure
plot3(dom(:,1),dom(:,2),dom(:,3),'.','color',0.6*ones(3,1))
while(y(end,1)>=2 || y(end,2)>=2)
	loop=loop-1;
  u=controller.control(y(end,:));
  
  %-------------here choose your controller input-------------%
  %in=u(end,:);
  in=u(end,:);
  %-----------------------------------------------------------%
  
  v=[v; in];
  [t x]=ode45(@unicycle_ode,[0 tau_s], y(end,:), odeset('abstol',1e-10,'reltol',1e-10),in');
    % if ( x(end,1)<=2 &&  x(end,2)<2)
    %     break
    % end
  y=[y; x(end,:)];
end


%% plot the vehicle domain
% colors
colors=get(groot,'DefaultAxesColorOrder');

% plot controller domain
dom=controller.domain;
plot3(dom(:,1),dom(:,2),dom(:,3),'.','color',0.6*ones(3,1))
hold on

% % plot initial state  and trajectory
% plot(y(:,1),y(:,2),'k.-')
% hold on
% plot(y(1,1),y(1,2),'.','color',colors(5,:),'markersize',20)

% % plot safe set
% eta_half=1.0/4e3;
% v=[1.15-eta_half 5.45-eta_half;...
%    1.55+eta_half 5.45-eta_half;...
%    1.15-eta_half 5.85+eta_half;...
%    1.55+eta_half 5.85+eta_half ];
% patch('vertices',v,'faces',[1 2 4 3],'facecolor','none','edgec',colors(2,:),'linew',1)
% hold on
% 
% box on
% axis([1.1 1.6 5.4 5.9])


%set(gcf,'paperunits','centimeters','paperposition',[0 0 16 10],'papersize',[16 10])

end

function dxdt = unicycle_ode(t,x,u)

            dxdt(1)=u(1)*cos(x(3));
			dxdt(2)=u(1)*sin(x(3));
            dxdt(3)=u(2);

end

