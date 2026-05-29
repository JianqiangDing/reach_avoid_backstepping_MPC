% Regenerate every example's controllers (constrained "result" + "unconstrained")
% under their current, intended settings (mu, FL sampling region, input bounds), so
% that any cross-controller comparison uses controllers built consistently.
%
% Each example script exports to ../controllers/ via export_to_python with stable,
% timestamp-free names:
%   sop_bounded_control_<sys>_result.py        (constrained / ours)
%   sop_bounded_control_<sys>_unconstrained.py (vanilla, no input constraints)
%
% NOTE: each example begins with `clear`, which wipes the caller's workspace; we
% therefore use literal sequential run() calls (no loop variable that must survive).
%
% Usage:  matlab -batch "addpath('matlab'); regenerate_all_controllers"

addpath(fileparts(mfilename('fullpath')));  % the matlab/ directory
fprintf('=== Regenerating all example controllers (intended conditions) ===\n');

fprintf('\n----- example_double_integrator -----\n');
run('example_double_integrator');

fprintf('\n----- example_dubins_car -----\n');
run('example_dubins_car');

fprintf('\n----- example_manipulator -----\n');
run('example_manipulator');

fprintf('\n=== Done. Controllers written to controllers/ ===\n');
