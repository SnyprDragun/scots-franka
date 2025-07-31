function franka
clear set
close all



%% simulation

% target set
lb1=[6 0];
ub1=lb1+0.5;
v1=[6 0; 6.5  0; 6 0.5; 6.5 .5];

lb2=[6 4];
ub2=lb2+0.5;
v2=[6 4; 6.5  4; 6 4.5; 6.5 4.5];
% initial state
x0_1=[0.4 0.4 0];
x0_2=[0.6 0.6 0];



controller1=SymbolicSet('franka_controller1.bdd','projection',[1 2 3]);
controller2=SymbolicSet('franka_controller2.bdd','projection',[1 2 3]);
target1=SymbolicSet('franka_target1.bdd');
target2=SymbolicSet('franka_target2.bdd');

y1=x0_1;
y2=x0_2;
v=[];
while(1)

  
  if ((target1.isElement(y1(end,:))) && (target2.isElement(y2(end,:))))
    break;
  end 

  u1=controller1.getInputs(y1(end,:));
  u2=controller2.getInputs(y2(end,:));

  v1=[v1; u1(1,:)];
  [t x1]=ode45(@unicycle_ode,[0 .3], y1(end,:), [],u1(1,:));

  y1=[y1; x1(end,:)];

  v2=[v2; u2(1,:)];
  [t x2]=ode45(@unicycle_ode,[0 .3], y2(end,:), [],u2(1,:));

  y2=[y2; x2(end,:)];
end



%% plot the franka domain
% colors
colors=get(groot,'DefaultAxesColorOrder');


% load the symbolic set containig the abstract state space
set=SymbolicSet('franka_ss.bdd','projection',[1 2]);
plotCells(set,'facecolor','none','edgec',[0.8 0.8 0.8],'linew',.1)
hold on

% load the symbolic set containig obstacles
set=SymbolicSet('franka_obst.bdd','projection',[1 2]);
plotCells(set,'facecolor',colors(1,:)*0.5+0.5,'edgec',colors(1,:),'linew',.1)

% plot the real obstacles and target set
plot_domain

% load the symbolic set containig target set
set=SymbolicSet('franka_target1.bdd','projection',[1 2]);
plotCells(set,'facecolor',colors(2,:)*0.5+0.5,'edgec',colors(2,:),'linew',.1)

set=SymbolicSet('franka_target2.bdd','projection',[1 2]);
plotCells(set,'facecolor',colors(2,:)*0.5+0.5,'edgec',colors(2,:),'linew',.1)

% plot initial state  and trajectory
plot(y1(:,1),y1(:,2),'k.-')
plot(y1(1,1),y1(1,1),'.','color',colors(5,:),'markersize',20)

plot(y2(:,1),y2(:,2),'k.-')
plot(y2(1,1),y2(1,1),'.','color',colors(5,:),'markersize',20)


box on
axis([-.5 10.5 -.5 10.5])


end

function dxdt = unicycle_ode(t,x,u)

  dxdt = zeros(3,1);
  c=atan(tan(u(2))/2);

  dxdt(1)=u(1)*cos(c+x(3))/cos(c);
  dxdt(2)=u(1)*sin(c+x(3))/cos(c);
  dxdt(3)=u(1)*tan(u(2));


end

function plot_domain

colors=get(groot,'DefaultAxesColorOrder');

v=[6 0; 6.5  0; 6 0.5; 6.5 .5];
patch('vertices',v,'faces',[1 2 4 3],'facec',colors(2,:),'edgec',colors(2,:));

v=[6 4; 6.5  4; 6 4.5; 6.5 4.5];
patch('vertices',v,'faces',[1 2 4 3],'facec',colors(2,:),'edgec',colors(2,:));

v=[2.2   0  ;2.4  0   ; 2.2   5    ; 2.4 5   ];
patch('vertices',v,'faces',[1 2 4 3],'facec',colors(1,:),'edgec',colors(1,:));
v=[4.6   1  ;4.8  1   ; 4.6   10   ; 4.8 10  ];
patch('vertices',v,'faces',[1 2 4 3],'facec',colors(1,:),'edgec',colors(1,:));

end

