%% page 4

[data, time] = ni2_activation;
figure; plot(time, data);

% verification that default values are according to the description
testdata = ni2_activation('frequency', 10, 'phase', 0, 'latency', 0.5, 'length', 1, 'ncycle', 5, 'powerup', 1, 'fsample', 1000);
figure; plot(time, data, time, testdata);

[data2, time2] = ni2_activation('latency',0.4,'frequency',5);
figure; hold on;
plot(time,data);
plot(time2,data2,'r');

%% page 5

% evaluation of different parameter settings
[testdata, testtime] = ni2_activation('frequency', 5,  'phase', pi, 'latency', 3, 'length', 5, 'ncycle', 12, 'powerup', 1, 'fsample', 1000);
figure; plot(testtime, testdata); hold on;

[testdata, testtime] = ni2_activation('frequency', 10, 'phase', 0, 'latency', 0.8, 'length', 1, 'ncycle', 10, 'powerup', 1, 'fsample', 1000);
plot(testtime, testdata);

datamix = data+data2;
figure; plot(time, datamix);

% plot the three time courses in a single figure
figure;
plot(time, [data; data2; datamix]);

datacombined = [data; data2];
figure; plot(time, datacombined);

datamix1 = data+data2;
datamix2 = sum(datacombined,1);
datamix3 = [1 1]*datacombined;

mix = [1 1];
datamix4 = mix*datacombined;

assert(isequal(datamix1, datamix4));
assert(isequal(datamix2, datamix4));
assert(isequal(datamix3, datamix4));

%% page 6

mix = [0 1; 0.1 0.9; 0.25 0.75; 0.5 0.5; 0.75 0.25; 0.9 0.1; 1 0];
datamix = mix * datacombined;
figure; plot(time, datamix+repmat((1:7)',[1 1000]));

% the effect of a negative mixing weight leads to a flip in the
% polarity of the contributing signal

%% page 8

headmodel = ni2_headmodel('type', 'spherical', 'nshell', 3);
sens = ni2_sensors('type', 'eeg');

figure; hold on;
ft_plot_headmodel(headmodel,'edgecolor','none');
ft_plot_sens(sens);

dippar1 = [0 0 6 1 0 0];
leadfield1 = ni2_leadfield(sens, headmodel, dippar1);
ni2_topoplot(sens, leadfield1); colorbar;

% the topography reflects a typical eeg dipolar potential distribution,
% with the line through the local maximum and minimum reflecting the
% right pointing orientation of the dipolar source

dippar2 = [0 0 6 2 0 0];
leadfield2 = ni2_leadfield(sens, headmodel, dippar2);
ni2_topoplot(sens, leadfield2); colorbar;

% verify that leadfield1 and leadfield2 differ by a factor of 2
figure; plot(leadfield1./leadfield2);

%% page 9

ni2_topoplot(sens, leadfield1*2);
ni2_topoplot(sens, leadfield2);

dippar3 = [0 0 5.9 1 0 0; 0 0 6.1 1 0 0];
leadfield3 = ni2_leadfield(sens, headmodel, dippar3);
topo = leadfield3*[1; 1];
ni2_topoplot(sens, topo); colorbar

% verify that the topography due to two closely spaced sources is very
% similar to a topography of a single source at the average location
testdippar1 = [0.75 0 6 0 1 0; 1.25 0 6 0 1 0];
testdippar2 = [1 0 6 0 2 0]; % this is the average location but double the dipole moment
testleadfield1 = ni2_leadfield(sens, headmodel, testdippar1);
testleadfield2 = ni2_leadfield(sens, headmodel, testdippar2);

ni2_topoplot(sens, sum(testleadfield1, 2)); colorbar
ni2_topoplot(sens, testleadfield2); colorbar
ni2_topoplot(sens, testleadfield2-sum(testleadfield1,2)); colorbar

%% page 10

% keep orientation fixed and move close to center of the head
z_coord = (0:0.5:8.5);
testdippar = zeros(numel(z_coord),6);
for k = 1:numel(z_coord)
  testdippar(k,:) = [0 0 z_coord(k) 1 0 0];
end
testleadfield = ni2_leadfield(sens, headmodel, testdippar);
for k = 1:size(testleadfield,2)
  ni2_topoplot(sens, testleadfield(:,k));
end

% change the y and z position of the source
r = 8;
angles = (-0.5:0.1:0.5).*pi;
testdippar = zeros(numel(angles),6);
for k = 1:numel(angles)
  % this results in a set of rightward pointing dipoles moving from
  % the 'back' of the head, through the vertex, to the 'front', at
  % a fixed depth
  testdippar(k,:) = [0 sin(angles(k)).*r cos(angles(k)).*r 1 0 0];
end
testleadfield = ni2_leadfield(sens, headmodel, testdippar);
for k = 1:size(testleadfield,2)
  ni2_topoplot(sens, testleadfield(:,k));
end

% taking the next three points together, investigating the effect of a
% change in orientation
testdippar = [0 0 8 1 0 0;
          0 0 8 0 1 0;
          0 0 8 0 0 1];
testleadfield = ni2_leadfield(sens, headmodel, testdippar);
for k = 1:size(testleadfield,2)
  ni2_topoplot(sens, testleadfield(:,k));
end

% rotate in the x/y plane
angles = (-1:0.1:1).*pi;
testdippar = zeros(numel(angles),6);
for k = 1:numel(angles)
  % this results in a set of rightward pointing dipoles moving from
  % the 'back' of the head, through the vertex, to the 'front', at
  % a fixed depth
  testdippar(k,:) = [0 0 8 sin(angles(k)) cos(angles(k)) 0];
end
testleadfield = ni2_leadfield(sens, headmodel, testdippar);
for k = 1:size(testleadfield,2)
  ni2_topoplot(sens, testleadfield(:,k));
end

% verify that linear mixing of the three-column leadfield is the same
% as providing the 'mixing' as dipole moment parameters
dipolemoment = randn(1,3);
testdippar = [0 0 8 dipolemoment];
testleadfield1 = ni2_leadfield(sens, headmodel, testdippar);
testleadfield2 = ni2_leadfield(sens, headmodel, testdippar(1:3));
ni2_topoplot(sens, testleadfield1);
ni2_topoplot(sens, testleadfield2*testdippar(4:6)');
