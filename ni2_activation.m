function [data, time] = ni2_activation(varargin)

% NI2_ACTIVATION creates a timecourse of simplified activity.
%
% Use as
%   [data, time] = ni2_activation
%   [data, time] = ni2_activation('key1', value1, 'key2', value2)
%
% Input arguments:
%   key-value pairs determine the shape of the activations.
%
%   frequency, scalar value (in Hz): frequency of oscillatory burst, default = 10 Hz   
%   phase,     scalar value (in radians): phase of the oscillation at the
%                peak latency, default = 0
%   latency,   scalar value (in s): latency at which the amplitude of the
%                oscillation peaks, default = 0.5
%   length,    scalar value (in s): length of the activation time course,
%                the sampling interval is fixed to 1 ms/sampl, default = 1
%   ncycle,    scalar value: number of oscillatory cycles in the burst, default = 5
%
% Output arguments:
%   data = 1xN vector with the activation time course
%   time = 1xN vector with the time specified in seconds

frequency = ft_getopt(varargin, 'frequency', 10);
phase     = ft_getopt(varargin, 'phase',     0);
latency   = ft_getopt(varargin, 'latency',   0.5);
length    = ft_getopt(varargin, 'length',    1);
ncycle    = ft_getopt(varargin, 'ncycle',    5);

maxncycle = 10*frequency*length;
if ncycle > maxncycle
  error('the maximum number of allowable cycles is exceeded');
end

% fix the sample frequency at 1000 Hz
fsample   = 1000;

% create time axis that is 10 times as long as the requested length
time      = (-round(5*length*fsample):round(5*length*fsample))./fsample;
n         = floor(numel(time)/10);

% create phase time course of oscillation
phs       = 2.*pi.*time.*frequency+phase;

% create a hanning window of ncycles length
hannkrn    = hanning(round(fsample.*ncycle./frequency));
nkrn       = numel(hannkrn);

% create a taper that is centered at time 0
taper      = [zeros(1,floor((10*n+1-nkrn)/2)) hannkrn(:)' zeros(1,ceil((10*n+1-nkrn)/2))];

% create the data
data       = cos(phs).*taper;

% cut off the edges to the requested length
time       = time+latency;
begsmp     = nearest(time,0);
time       = time(begsmp+(1:n));
data       = data(begsmp+(1:n));
