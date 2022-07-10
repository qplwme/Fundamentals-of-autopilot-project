clear, close all
model_parameter
addpath('splines');

%% lateral direction
vx = 0.5;

%% load and init MPC_vars
load('MPC_vars.mat');

nx = ModelParams.nx;
nu = ModelParams.nu;
N = MPC_vars.N;
Ts = MPC_vars.Ts;

%% Simulation lenght and plotting
simN = 500;
%0=no plots, 1=plot predictions
plotOn = 1;
%0=real time iteration, 1=fixed number of QP iterations, 2=fixed number of damped QP iterations
QP_iter = 2;

% number of cars 
% there are two examples one with no other cars and one with 4 other cars
% (inspired by the set up shown in the paper)
n_cars = 1; % no other car
% n_cars = 5; % 4 other cars


%% add trajectory
load('traj_diy.mat');
load('MPC_vars.mat');

% scale
scale = 10;
traj = traj / scale;

track2_origin.outer = [traj(7,:); traj(8,:)];
track2_origin.inner = [traj(5,:); traj(6,:)];
track2_origin.center = [traj(1,:); traj(2,:)];

safteyScaling = 0.8;
[track_origin,track2_origin] = borderAdjustment(track2_origin,ModelParams,safteyScaling);

%% create virtual path
theta = -pi/2 + atan((track2_origin.outer(2,995) - track2_origin.inner(2,995)) / (track2_origin.outer(1,995) - track2_origin.inner(1,995)));
[line_k, b] = polyfit([track2_origin.outer(1,995) track2_origin.inner(1,995)], [track2_origin.outer(2,995) track2_origin.inner(2,995)], 1);
virture_track.outer(1,1) = track2_origin.outer(1,length(track2_origin.outer)) + 0.02*cos(theta);
for i=2:400
    virture_track.outer(1,i) = virture_track.outer(1,i-1) + 0.02*cos(theta);
end
virture_track.outer(2,1) = track2_origin.outer(2,length(track2_origin.outer)) + 0.02*sin(theta);
for i=2:400
    virture_track.outer(2,i) = virture_track.outer(2,i-1) + 0.02*sin(theta);
end

virture_track.inner = track2_origin.inner(1,length(track2_origin.inner)) + 0.02*cos(theta);
for i=2:400
    virture_track.inner(1,i) = virture_track.inner(1,i-1) + 0.02*cos(theta);
end
virture_track.inner(2,1) = track2_origin.inner(2,length(track2_origin.outer)) + 0.02*sin(theta);
for i=2:400
    virture_track.inner(2,i) = virture_track.inner(2,i-1) + 0.02*sin(theta);
end
virture_track.center = (virture_track.inner + virture_track.outer) / 2;

%% paste origin path and virtual path
track2.outer = [track2_origin.outer virture_track.outer];
track2.inner = [track2_origin.inner virture_track.inner];
track2.center = [track2_origin.center virture_track.center];

%% Fit spline to track
% TODO spline function only works with regular spaced points.
% Fix add function which given any center line and bound generates equlally
% space tracks.

[track,track2] = borderAdjustment(track2,ModelParams,safteyScaling);

[traj, borders] =splinify(track);
tl = traj.ppy.breaks(end);

% store all data in one struct
TrackMPC = struct('traj',traj,'borders',borders,'track_center',track.center,'tl',tl);

%% Define starting position
startIdx = 1; %point (in terms of track centerline array) allong the track 
% where the car starts, on the center line, aligned with the track, driving
% straight with vx0
%since the used bicycle model is not well defined for slow velocities use vx0 > 0.5
vx0 = vx;
trackWidth = norm(track.inner(:,1)-track.outer(:,1));

% find theta that coresponds to the 10th point on the centerline
[theta, ~] = findTheta([track.center(1,startIdx),track.center(2,startIdx)],track.center,traj.ppx.breaks,trackWidth,startIdx);

x0 = [track.center(1,startIdx),track.center(2,startIdx),... % point on centerline
      atan2(ppval(traj.dppy,theta),ppval(traj.dppx,theta)),... % aligned with centerline
      vx0 ,0,0,theta]'; %driving straight with vx0, and correct theta progress
    
% the find theta function performs a local search to find the projection of
% the position onto the centerline, therefore, we use the start index as an
% starting point for this local search
last_closestIdx = startIdx;

%% First initial guess
x = repmat(x0,1,N+1); % all points identical to current measurment
% first inital guess, all points on centerline aligned with centerline
% spaced as the car would drive with vx0
for i = 2:N+1
    theta_next = x(ModelParams.stateindex_theta,i-1)+Ts*vx0;
    phi_next = atan2(ppval(traj.dppy,theta_next),ppval(traj.dppx,theta_next));
    % phi_next can jump by two pi, make sure there are no jumps in the
    % initial guess
    if (x(ModelParams.stateindex_phi,i-1)-phi_next) < -pi
        phi_next = phi_next-2*pi;
    end
    if (x(ModelParams.stateindex_phi,i-1)-phi_next) > pi
        phi_next = phi_next+2*pi;
    end
    x(:,i) = [ppval(traj.ppx,theta_next),ppval(traj.ppy,theta_next),... % point on centerline
              phi_next,... % aligned with centerline
              vx0 ,0,0,theta_next]'; %driving straight with vx0, and correct theta progress
end

u = zeros(3,N); % zero inputs
uprev = zeros(3,1); % last input is zero

%% Ohter cars
Y = [];

%% Initialize logging arrays
X_log = zeros(nx*(N+1),simN);
U_log = zeros(3*N,simN);
B_log = zeros(4*N,simN);
qpTime_log = zeros(1,simN);

%% initializtion for MPCC
% solve problem 5 times without applying input
% inspiered by sequential quadratic programming (SQP)
for i = 1:5
    % formulate MPCC problem and solve it
    Iter_damping = 0.5; % 0 no damping
    [x_up, u_up, b, exitflag,info] = optimizer_mpcc(TrackMPC,MPC_vars,ModelParams, n_cars, Y, x, u, x0, uprev);
    x = Iter_damping*x + (1-Iter_damping)*x_up;
    u = Iter_damping*u + (1-Iter_damping)*u_up;

    if plotOn == 1
        % plot predictions
        PlotPrediction(x,u,b,Y,track,track2,traj,MPC_vars,ModelParams)
    end
end