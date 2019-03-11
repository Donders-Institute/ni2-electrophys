function [data, time] = ni2_activation_oscillation(varargin)

% NI2_ACTIVATION_OSCILLATION creates a timecourse of simplified oscillatory activity.
%
% Use as
%   [data, time] = ni2_activation_oscillation
%   [data, time] = ni2_activation_oscillation('key1', value1, 'key2', value2)
%
% Input arguments:
%   key-value pairs determine the shape of the activations.
%
%   frequency, two element vector (in Hz): frequency band of oscillatory activity   
%   length,    scalar value (in s): length of the activation time course,
%                the sampling interval is fixed to 1 ms/sampl, default = 1
%   fsample,   scalar value: sampling rate (Hz) (default 1000 Hz)
%   envelope,  1xnsmp envelope of the amplitude
%   nsignal,   scalar value: number of signals (default =1)
%
% Output arguments:
%   data = nsignalxN vector with the activation time course
%   time = 1xN vector with the time specified in seconds

frequency = ft_getopt(varargin, 'frequency', [8 12]);
length    = ft_getopt(varargin, 'length',    1);
fsample   = ft_getopt(varargin, 'fsample',   1000);
envelope  = ft_getopt(varargin, 'envelope',  ones(1,round(length*fsample)));
nsignal   = ft_getopt(varargin, 'nsignal',   1);

if size(envelope,1)==1
  envelope = repmat(envelope, [nsignal 1]);
end

% create time axis that is 10 times as long as the requested length
time      = (-round(5*length*fsample):round(5*length*fsample))./fsample;
n         = floor(numel(time)/10);

data = ft_preproc_bandpassfilter(randn(nsignal,numel(time)), fsample, frequency, [], 'firws');

% cut off the edges to the requested length
begsmp     = nearest(time,-length./2);
time       = time(begsmp+(1:n));
data       = data(:,begsmp+(1:n)).*envelope;
